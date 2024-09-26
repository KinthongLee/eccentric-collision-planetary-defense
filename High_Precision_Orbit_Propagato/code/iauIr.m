%  - - - - - -
%   i a u I r
%  - - - - - -
%
%  Initialize an r-matrix to the identity matrix.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  vector/matrix support function.
%
%  Returned:
%     r       double[3][3]    r-matrix
%
%  This revision:  2021 May 11
%
%  SOFA release 2021-05-12
%
%  Copyright (C) 2021 IAU SOFA Board.  See notes at end.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = iauIr(r)

r(1,1) = 1.0;
r(1,2) = 0.0;
r(1,3) = 0.0;
r(2,1) = 0.0;
r(2,2) = 1.0;
r(2,3) = 0.0;
r(3,1) = 0.0;
r(3,2) = 0.0;
r(3,3) = 1.0;

