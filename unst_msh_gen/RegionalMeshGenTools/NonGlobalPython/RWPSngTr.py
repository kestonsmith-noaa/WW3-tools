import os
import argparse
import time
import numpy as np

import jigsawpy


#def case_6_(src_path, dst_path):

# DEMO-6: generate a 2-dim. grid for the Australian coastal
# region, using scaled ocean-depth as a mesh-resolution
# heuristic. A local stereographic projection is employed.

opts = jigsawpy.jigsaw_jig_t()
topo = jigsawpy.jigsaw_msh_t()
geom = jigsawpy.jigsaw_msh_t()
mesh = jigsawpy.jigsaw_msh_t()
hmat = jigsawpy.jigsaw_msh_t()
proj = jigsawpy.jigsaw_prj_t()

#------------------------------------ setup files for JIGSAW

    #opts.geom_file = os.path.join(dst_path, "geom.msh")
    #opts.jcfg_file = os.path.join(dst_path, "aust.jig")
    #opts.mesh_file = os.path.join(dst_path, "mesh.msh")
    #opts.hfun_file = os.path.join(dst_path, "spac.msh")
    
opts.geom_file = "geom.msh"
opts.jcfg_file = "aust.jig"
opts.mesh_file = "mesh.msh"
opts.hfun_file = "spac.msh"

#------------------------------------ define JIGSAW geometry

jigsawpy.loadmsh("PSLGboundary.msh", geom)
jigsawpy.loadmsh("topo_shifted.msh", topo)

xmin = np.min(
        geom.point["coord"][:, 0])
ymin = np.min(
        geom.point["coord"][:, 1])
xmax = np.max(
        geom.point["coord"][:, 0])
ymax = np.max(
        geom.point["coord"][:, 1])

zlev = topo.value

xmsk = np.logical_and(topo.xgrid > xmin,
                          topo.xgrid < xmax)
ymsk = np.logical_and(topo.ygrid > ymin,
                          topo.ygrid < ymax)

zlev = zlev[:, xmsk]
zlev = zlev[ymsk, :]

#------------------------------------ define spacing pattern

hmat.mshID = "ellipsoid-grid"
hmat.radii = np.full(
        +3, +6371.0,
        dtype=jigsawpy.jigsaw_msh_t.REALS_t)

hmat.xgrid = \
        topo.xgrid[xmsk] * np.pi / 180.
hmat.ygrid = \
        topo.ygrid[ymsk] * np.pi / 180.

#    hmin = +1.0E+01; hmax = +1.0E+02
hmin = +10.; hmax = +100.

#    hmat.value = \
#        np.sqrt(np.maximum(-zlev, 0.)) / 0.5
hmat.value = \
        50.*np.sqrt(np.maximum(-zlev, 0.)) / 0.5

hmat.value = \
        np.maximum(hmat.value, hmin)
hmat.value = \
        np.minimum(hmat.value, hmax)

hmat.slope = np.full(
        hmat.value.shape, +0.1500,
        dtype=jigsawpy.jigsaw_msh_t.REALS_t)

#------------------------------------ do stereographic proj.

geom.point["coord"][:, :] *= np.pi / 180.

proj.prjID = 'stereographic'
proj.radii = +6.371E+003
proj.xbase = \
        +0.500 * (xmin + xmax) * np.pi / 180.
proj.ybase = \
        +0.500 * (ymin + ymax) * np.pi / 180.

jigsawpy.project(geom, proj, "fwd")
jigsawpy.project(hmat, proj, "fwd")

jigsawpy.savemsh(opts.geom_file, geom)
jigsawpy.savemsh(opts.hfun_file, hmat)

#------------------------------------ set HFUN grad.-limiter

jigsawpy.cmd.marche(opts, hmat)

#------------------------------------ make mesh using JIGSAW

opts.hfun_scal = "absolute"
opts.hfun_hmax = float("inf")       # null HFUN limits
opts.hfun_hmin = float(+0.00)

opts.mesh_dims = +2                 # 2-dim. simplexes
opts.mesh_eps1 = +1.

ttic = time.time()

jigsawpy.cmd.jigsaw(opts, mesh)

ttoc = time.time()

print("CPUSEC =", (ttoc - ttic))

cost = jigsawpy.triscr2(            # quality metrics!
        mesh.point["coord"],
        mesh.tria3["index"])

print("TRISCR =", np.min(cost), np.mean(cost))

jigsawpy.savemsh("RWPSrad.msh",mesh)

print("PWRSCR =", np.min(cost), np.mean(cost))

tbad = jigsawpy.centre2(
        mesh.point["coord"],
        mesh.power,
        mesh.tria3["index"])

print("OBTUSE =",
          +np.count_nonzero(np.logical_not(tbad)))

#------------------------------------ save mesh for Paraview
jigsawpy.savemsh("RWPS.msh",mesh)
jigsawpy.project(mesh, proj, "inv")

jigsawpy.savemsh("RWPSrad.msh",mesh)

mesh.point["coord"][:, :] = mesh.point["coord"][:, :]*180. / np.pi

jigsawpy.savemsh("RWPSlatlon_shifted.msh",mesh)

print("Saving to ../cache/case_6a.vtk")

jigsawpy.savevtk("case_6a.vtk", mesh)

print("Saving to ../cache/case_6b.vtk")

jigsawpy.savevtk("case_6b.vtk", hmat)
