%  - - - - - -
%   i a u R y
%  - - - - - -
%
%  Rotate an r-matrix about the y-axis.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  vector/matrix support function.
%
%  Given:
%     theta  double          angle (radians)
%
%  Given and returned:
%     r      double[3][3]    r-matrix, rotated
%
%  Notes:
%
%  1) Calling this function with positive theta incorporates in the
%     supplied r-matrix r an additional rotation, about the y-axis,
%     anticlockwise as seen looking towards the origin from positive y.
%
%  2) The additional rotation can be represented by this matrix:
%
%         (  + cos(theta)     0      - sin(theta)  )
%         (                                        )
%         (       0           1           0        )
%         (                                        )
%         (  + sin(theta)     0      + cos(theta)  )
%
%  This revision:  2021 May 11
%
%  SOFA release 2021-05-12
%
%  Copyright (C) 2021 IAU SOFA Board.  See notes at end.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = iauRy(theta, r)

s = sin(theta);
c = cos(theta);

a00 = c*r(1,1) - s*r(3,1);
a01 = c*r(1,2) - s*r(3,2);
a02 = c*r(1,3) - s*r(3,3);
a20 = s*r(1,1) + c*r(3,1);
a21 = s*r(1,2) + c*r(3,2);
a22 = s*r(1,3) + c*r(3,3);

r(1,1) = a00;
r(1,2) = a01;
r(1,3) = a02;
r(3,1) = a20;
r(3,2) = a21;
r(3,3) = a22;

