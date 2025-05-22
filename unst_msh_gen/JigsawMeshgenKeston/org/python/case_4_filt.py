
import os
import time
import math
import numpy as np
from scipy import interpolate

import jigsawpy

import os
import argparse



# DEMO-4: generate a multi-resolution mesh, via local refin-
# ement along coastlines and shallow ridges. Global grid
# resolution is 150KM, background resolution is 67KM and the
# min. adaptive resolution is 33KM.

opts = jigsawpy.jigsaw_jig_t()
topo = jigsawpy.jigsaw_msh_t()
geom = jigsawpy.jigsaw_msh_t()
mesh = jigsawpy.jigsaw_msh_t()
hmat = jigsawpy.jigsaw_msh_t()

#------------------------------------ setup files for JIGSAW

opts.geom_file = "geom.msh"
opts.jcfg_file = "opts.jig"
opts.mesh_file = "mesh.msh"
opts.hfun_file = "spac.msh"

#------------------------------------ define JIGSAW geometry

geom.mshID = "ellipsoid-mesh"
geom.radii = np.full(3, 6.371E+003, dtype=geom.REALS_t)
jigsawpy.savemsh(opts.geom_file, geom)

#------------------------------------ define spacing pattern

    #jigsawpy.loadmsh(os.path.join(
    #    src_path, "topo.msh"), topo)
jigsawpy.loadmsh("hfunD.msh", topo)
hmat.mshID = "ellipsoid-grid"
hmat.radii = geom.radii

hmat.xgrid = topo.xgrid * np.pi / 180.
hmat.ygrid = topo.ygrid * np.pi / 180.

#    hfn0 = +150.                        # global spacing
#    hfn2 = +33.                         # adapt. spacing
#    hfn3 = +67.                         # arctic spacing
#    hmat.value = np.sqrt(np.maximum(-topo.value, 0.0))
#    hmat.value = np.maximum(hmat.value, hfn2)
#    hmat.value = np.minimum(hmat.value, hfn3)

    #mask = hmat.ygrid < 40. * np.pi / 180.

  #  hmat.value[mask] = hfn0

hmat.value=topo.value

#------------------------------------ set HFUN grad.-limiter

hmat.slope = np.full(   topo.value.shape,  +0.050, dtype=hmat.REALS_t)           # |dH/dx| limits

jigsawpy.savemsh(opts.hfun_file, hmat)

jigsawpy.cmd.marche(opts, hmat)

#------------------------------------ make mesh using JIGSAW

opts.hfun_scal = "absolute"
opts.hfun_hmax = float("inf")       # null HFUN limits
opts.hfun_hmin = float(+0.00)

opts.mesh_dims = +2                 # 2-dim. simplexes

opts.optm_qlim = +9.5E-01           # tighter opt. tol
opts.optm_iter = +32
opts.optm_qtol = +1.0E-05

#   opts.optm_kern = "cvt+dqdx"

rbar = np.mean(geom.radii)          # bisect heuristic
hbar = np.mean(hmat.value)
nlev = round(math.log2(rbar / math.sin(.4 * math.pi) / hbar))

ttic = time.time()

jigsawpy.cmd.tetris(opts, nlev - 1, mesh)

ttoc = time.time()

print("CPUSEC =", (ttoc - ttic))

print("BISECT =", +nlev)

cost = jigsawpy.triscr2(            # quality metrics!
        mesh.point["coord"],
        mesh.tria3["index"])

print("TRISCR =", np.min(cost), np.mean(cost))

cost = jigsawpy.pwrscr2(
        mesh.point["coord"],
        mesh.power,
        mesh.tria3["index"])

print("PWRSCR =", np.min(cost), np.mean(cost))

tbad = jigsawpy.centre2(
        mesh.point["coord"],
        mesh.power,
        mesh.tria3["index"])

print("OBTUSE =",
          +np.count_nonzero(np.logical_not(tbad)))

ndeg = jigsawpy.trideg2(
        mesh.point["coord"],
        mesh.tria3["index"])

print("TOPOL. =",+np.count_nonzero(ndeg==+6) / ndeg.size)

#------------------------------------ save mesh for Paraview

apos = jigsawpy.R3toS2(
        geom.radii, mesh.point["coord"][:])

apos = apos * 180. / np.pi

zfun = interpolate.RectBivariateSpline(
        topo.ygrid, topo.xgrid, topo.value)

mesh.value = zfun(
        apos[:, 1], apos[:, 0], grid=False)

cell = mesh.tria3["index"]

zmsk = \
        mesh.value[cell[:, 0]] + \
        mesh.value[cell[:, 1]] + \
        mesh.value[cell[:, 2]]
zmsk = zmsk / +3.0

mesh.tria3 = mesh.tria3[zmsk < +0.]



jigsawpy.savemsh("RWPS_c4.msh",mesh)
#jigsawpy.project(mesh, proj, "inv")

#jigsawpy.savemsh("RWPS_c4Rad.msh",mesh)

mesh.point["coord"][:, :] = mesh.point["coord"][:, :]*180. / np.pi

#jigsawpy.savemsh("RWPS_c4LL.msh",mesh)

filter_ocn()


print("Saving to ../cache/case_4a.vtk")

jigsawpy.savevtk(os.path.join(
        dst_path, "case_4a.vtk"), mesh)

print("Saving to ../cache/case_4b.vtk")

jigsawpy.savevtk(os.path.join(
        dst_path, "case_4b.vtk"), hmat)





def filter_dry(mesh, mask):

#-- require dry cells > 1 dry edge, and large

    print("*filter-dry...")

    filt = np.logical_not(mask)
    tris = np.argwhere(filt).ravel()

    conn = tri_to_tri(mesh.tria3["index"][filt, :])
        
    # require dry to be adj. >=1 dry cell
    isol = np.sum(conn, axis=1) <= 1
    isol = np.ravel(isol)

    """
    # delete groups of dry if too small
    nprt, part = connected_components(
        conn, directed=False, return_labels=True)

    tris = np.argwhere(filt).ravel()
    for iprt in range(nprt):
        itri = np.argwhere(part == iprt)
        if (itri.size <= 2): mask[tris[itri]] = True
    """

    # otherwise mark isolated cell as ocn
    mask[tris[isol]] = True

    return mask


def filter_wet(mesh, mask):

#-- require wet cells > 1 wet edge, and large

    print("*filter-wet...")

    tris = np.argwhere(mask).ravel()

    conn = tri_to_tri(mesh.tria3["index"][mask, :])

    # require wet to be adj. >=1 wet cell
    isol = np.sum(conn, axis=1) <= 1
    isol = np.ravel(isol)
    
    # delete groups of wet if too small
    nprt, part = connected_components(
        conn, directed=False, return_labels=True)

    area = jigsawpy.trivol2(
        mesh.point["coord"], 
        mesh.tria3["index"][mask, :])

    for iprt in range(nprt):
        itri = np.argwhere(part == iprt)
        asum = np.sum(area[itri])
       #print(asum)
        if asum < ISOLATED: mask[tris[itri]] = False
    
    # otherwise mark isolated cell as dry
    mask[tris[isol]] = False
    
    return mask


def filter_ocn():

    args = parse_input_args()
    configurations = load_configuration(args.config)

#-- use the remapped elev. to keep ocean cells

    print("*filter-ocn...")

    elev =(mesh.value[mesh.tria3["index"][:, 0]]
         + mesh.value[mesh.tria3["index"][:, 1]]
         + mesh.value[mesh.tria3["index"][:, 2]]
         + mesh.vmids) / 4.0
    # Define the Caspian Sea region
    caspian_lat_min = 34.5
    caspian_lat_max = 50.0
    caspian_lon_min = 44.5
    caspian_lon_max = 55.5
    
    # Define the black Sea region
    blacksea_lat_min = 40 
    blacksea_lat_max = 47.25
    blacksea_lon_min = 26.15
    blacksea_lon_max = 41.5
    
   # Define the additional region1
    additional_lat_min = 39.95
    additional_lat_max = 40.6
    additional_lon_min = 26
    additional_lon_max = 26.8

   # Define the additional region2
    additional_lat_min2 = 40.3
    additional_lat_max2 = 40.6
    additional_lon_min2 = 26.8
    additional_lon_max2 = 30

   # Define the additional region3
    additional_lat_min3 = 40.6
    additional_lat_max3 = 41.25
    additional_lon_min3 = 28.9
    additional_lon_max3 = 29.1
    
    """
# Print the elevations in the additional region
    additional_elev = elev[np.logical_and.reduce((
        mesh.smids[:, 1] >= additional_lat_min,
        mesh.smids[:, 1] <= additional_lat_max,
        mesh.smids[:, 0] >= additional_lon_min,
        mesh.smids[:, 0] <= additional_lon_max
    ))]
    print("Elevations in the additional region:")
    print(additional_elev)

# Print the elevations in the additional region2
    additional_elev2 = elev[np.logical_and.reduce((
        mesh.smids[:, 1] >= additional_lat_min2,
        mesh.smids[:, 1] <= additional_lat_max2,
        mesh.smids[:, 0] >= additional_lon_min2,
        mesh.smids[:, 0] <= additional_lon_max2
    ))]
    print("Elevations in the additional region2:")
    print(additional_elev2)

# Print the elevations in the additional region3
    additional_elev3 = elev[np.logical_and.reduce((
        mesh.smids[:, 1] >= additional_lat_min3,
        mesh.smids[:, 1] <= additional_lat_max3,
        mesh.smids[:, 0] >= additional_lon_min3,
        mesh.smids[:, 0] <= additional_lon_max3
    ))]
    print("Elevations in the additional region3:")
    print(additional_elev3)

    """

    # zssh, to cull elev. against
    surf = np.zeros(elev.shape, dtype=np.float32)
    # Update the surf array to include both regions
    # Define the Caspian Sea region
    caspian_region = np.logical_and.reduce((
        mesh.smids[:, 1] >= caspian_lat_min,
        mesh.smids[:, 1] <= caspian_lat_max,
        mesh.smids[:, 0] >= caspian_lon_min,
        mesh.smids[:, 0] <= caspian_lon_max
    ))
    
    blacksea_region = np.logical_and.reduce((
        mesh.smids[:, 1] >= blacksea_lat_min,
        mesh.smids[:, 1] <= blacksea_lat_max,
        mesh.smids[:, 0] >= blacksea_lon_min,
        mesh.smids[:, 0] <= blacksea_lon_max
    ))
    
    black_sea = configurations['black_sea']
        # Activate regions based on black_sea option
    if black_sea == 1:  # Caspian and Black Sea
        surf[caspian_region] = -9999.0
        surf[blacksea_region] = -9999.0
    elif black_sea == 2:  # Only Caspian Sea
        surf[caspian_region] = -9999.0
    elif black_sea == 3:  # All except Black Sea
        surf[caspian_region] = -9999.0
        # Additional regions
        additional_region = np.logical_and.reduce((
            mesh.smids[:, 1] >= additional_lat_min,
            mesh.smids[:, 1] <= additional_lat_max,
            mesh.smids[:, 0] >= additional_lon_min,
            mesh.smids[:, 0] <= additional_lon_max
        ))
        surf[additional_region] = 200.

        additional_region2 = np.logical_and.reduce((
            mesh.smids[:, 1] >= additional_lat_min2,
            mesh.smids[:, 1] <= additional_lat_max2,
            mesh.smids[:, 0] >= additional_lon_min2,
            mesh.smids[:, 0] <= additional_lon_max2
        ))
        surf[additional_region2] = 300.

        additional_region3 = np.logical_and.reduce((
            mesh.smids[:, 1] >= additional_lat_min3,
            mesh.smids[:, 1] <= additional_lat_max3,
            mesh.smids[:, 0] >= additional_lon_min3,
            mesh.smids[:, 0] <= additional_lon_max3
        ))
        surf[additional_region3] = 400.
        
    keep = elev <= surf  # only keep tri with wet elev
   
    # iterate on dry cells until none "isolated" 
    knum = np.count_nonzero(keep)
    while (True):   
        keep = filter_dry(mesh, keep)
        if (np.count_nonzero(keep) == knum): break
        knum = np.count_nonzero(keep)
    
    mesh.tria3 = mesh.tria3[keep]

    # iterate on wet cells until none "isolated"
    keep = np.ones(mesh.tria3.size, dtype=bool)
    knum = np.count_nonzero(keep)
    while (True):
        keep = filter_wet(mesh, keep)
        if (np.count_nonzero(keep) == knum): break
        knum = np.count_nonzero(keep)
    
    mesh.tria3 = mesh.tria3[keep]

    # delete unused vertices and reindex
    ifwd = np.unique(mesh.tria3["index"].ravel())
    
    irev = np.zeros(mesh.point.size, dtype=np.int32)
    irev[ifwd] = np.arange(ifwd.size, dtype=np.int32)

    mesh.point = mesh.point[ifwd]
    mesh.value = mesh.value[ifwd]
    mesh.tria3["index"] = irev[mesh.tria3["index"]]
    return mesh



def tri_to_tri(tria):

#-- return tria-to-tria adj. as a sparse graph

    # non-unique edges in tris
    edge = np.empty((0, 2), dtype=np.int32)
    tris = np.empty((0), dtype=np.int32)
    edge = np.concatenate((edge, 
        tria[:, (0, 1)]), axis=0)
    tris = np.concatenate((tris, 
        np.arange(0, tria.shape[0])))
        
    edge = np.concatenate((edge, 
        tria[:, (1, 2)]), axis=0)
    tris = np.concatenate((tris, 
        np.arange(0, tria.shape[0])))
        
    edge = np.concatenate((edge, 
        tria[:, (2, 0)]), axis=0)
    tris = np.concatenate((tris, 
        np.arange(0, tria.shape[0])))
        
    # which edges match to which?
    edge = np.sort(edge, axis=1)
    imap = np.argsort(edge[:, 1], kind="stable")
    edge = edge[imap, :]
    tris = tris[imap]
    imap = np.argsort(edge[:, 0], kind="stable")
    edge = edge[imap, :]
    tris = tris[imap]
    
    diff = edge[1::, :] - edge[:-1:, :]

    same = np.argwhere(np.logical_and.reduce((
        diff[:, 0] == 0, 
        diff[:, 1] == 0))).ravel()
        
    # tris[same] and tris[same+1] share
    rows = np.concatenate((
        tris[same], tris[same+1]))
    cols = np.concatenate((
        tris[same+1], tris[same]))
    data = np.ones(rows.size, dtype=np.int8)
    
    # ith tri is adj. to tri in ith row
    return csr_matrix((data, (rows, cols)))


