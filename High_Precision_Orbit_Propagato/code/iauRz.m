%  - - - - - -
%   i a u R z
%  - - - - - -
%
%  Rotate an r-matrix about the z-axis.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  vector/matrix support function.
%
%  Given:
%     psi    double          angle (radians)
%
%  Given and returned:
%     r      double[3][3]    r-matrix, rotated
%
%  Notes:
%
%  1) Calling this function with positive psi incorporates in the
%     supplied r-matrix r an additional rotation, about the z-axis,
%     anticlockwise as seen looking towards the origin from positive z.
%
%  2) The additional rotation can be represented by this matrix:
%
%         (  + cos(psi)   + sin(psi)     0  )
%         (                                 )
%         (  - sin(psi)   + cos(psi)     0  )
%         (                                 )
%         (       0            0         1  )
%
%  This revision:  2021 May 11
%
%  SOFA release 2021-05-12
%
%  Copyright (C) 2021 IAU SOFA Board.  See notes at end.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = iauRz(psi, r)

s = sin(psi);
c = cos(psi);

a00 =   c*r(1,1) + s*r(2,1);
a01 =   c*r(1,2) + s*r(2,2);
a02 =   c*r(1,3) + s*r(2,3);
a10 = - s*r(1,1) + c*r(2,1);
a11 = - s*r(1,2) + c*r(2,2);
a12 = - s*r(1,3) + c*r(2,3);

r(1,1) = a00;
r(1,2) = a01;
r(1,3) = a02;
r(2,1) = a10;
r(2,2) = a11;
r(2,3) = a12;
