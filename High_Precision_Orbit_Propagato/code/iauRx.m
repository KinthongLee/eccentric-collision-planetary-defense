%  - - - - - -
%   i a u R x
%  - - - - - -
%
%  Rotate an r-matrix about the x-axis.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  vector/matrix support function.
%
%  Given:
%     phi    double          angle (radians)
%
%  Given and returned:
%     r      double[3][3]    r-matrix, rotated
%
%  Notes:
%
%  1) Calling this function with positive phi incorporates in the
%     supplied r-matrix r an additional rotation, about the x-axis,
%     anticlockwise as seen looking towards the origin from positive x.
%
%  2) The additional rotation can be represented by this matrix:
%
%         (  1        0            0      )
%         (                               )
%         (  0   + cos(phi)   + sin(phi)  )
%         (                               )
%         (  0   - sin(phi)   + cos(phi)  )
%
%  This revision:  2021 May 11
%
%  SOFA release 2021-05-12
%
%  Copyright (C) 2021 IAU SOFA Board.  See notes at end.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = iauRx(phi, r)

s = sin(phi);
c = cos(phi);

a10 =   c*r(2,1) + s*r(3,1);
a11 =   c*r(2,2) + s*r(3,2);
a12 =   c*r(2,3) + s*r(3,3);
a20 = - s*r(2,1) + c*r(3,1);
a21 = - s*r(2,2) + c*r(3,2);
a22 = - s*r(2,3) + c*r(3,3);

r(2,1) = a10;
r(2,2) = a11;
r(2,3) = a12;
r(3,1) = a20;
r(3,2) = a21;
r(3,3) = a22;

