% DETBNDY compute a boundary segment list for a FEM domain
%
% DETBNDY bnd=detbndy(e);
%         This function computes a boundary for the FEM domain
%         described a file containing element connectivity list (e).
%         It uses sparse matrix techniques to determine the element
%         edges on the boundary of the FEM domain.
%
% Input:  ele -  element list; 3 (.tri) or 4 (.ele) columns wide
% Output: bnd -  a 2-column list of boundary-node numbers, returned
%                to the local workspace
%
%         The output boundary list are pairs of node numbers, not 
%         coordinates, describing the edges of elements on the 
%         exterior of the domain, including islands.  The segments 
%         are not connected.
%
%         Call as: bnd=detbndy(e);
%
% Written by : Brian O. Blanton at The University of North Carolina 
%              at Chapel Hill, Mar 1995.
%
function bnd=detbndy(e)
 
 
% Check size of element list
[nelems,ncol]=size(e);

% Form (i,j) connection list from .ele element list
%
i=[e(:,1);e(:,2);e(:,3)];
j=[e(:,2);e(:,3);e(:,1)];

% Form the sparse adjacency matrix and add transpose.
%
n = max(max(i),max(j));
A = sparse(i,j,-1,n,n);
A = A + A';

% Consider only the upper part of A, since A is symmetric
% 
A=A.*triu(A);

% The boundary segments are A's with value == 1
%
B=A==1;

% Extract the row,col from B for the boundary list.
%
[ib,jb,s]=find(B);
bnd=[ib(:),jb(:)];
