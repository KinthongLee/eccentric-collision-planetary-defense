%  - - - - - - - - -
%   i a u O b l 0 6
%  - - - - - - - - -
%
%  Mean obliquity of the ecliptic, IAU 2006 precession model.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  canonical model.
%
%  Given:
%     date1,date2  double   TT as a 2-part Julian Date (Note 1)
%
%  Returned (function value):
%                  double   obliquity of the ecliptic (radians, Note 2)
%
%  Notes:
%
%  1) The TT date date1+date2 is a Julian Date, apportioned in any
%     convenient way between the two arguments.  For example,
%     JD(TT)=2450123.7 could be expressed in any of these ways,
%     among others:
%
%            date1          date2
%
%         2450123.7           0.0       (JD method)
%         2451545.0       -1421.3       (J2000 method)
%         2400000.5       50123.2       (MJD method)
%         2450123.5           0.2       (date & time method)
%
%     The JD method is the most natural and convenient to use in
%     cases where the loss of several decimal digits of resolution
%     is acceptable.  The J2000 method is best matched to the way
%     the argument is handled internally and will deliver the
%     optimum resolution.  The MJD method and the date & time methods
%     are both good compromises between resolution and convenience.
%
%  2) The result is the angle between the ecliptic and mean equator of
%     date date1+date2.
%
%  Reference:
%
%     Hilton, J. et al., 2006, Celest.Mech.Dyn.Astron. 94, 351
%
%  This revision:  2021 May 11
%
%  SOFA release 2021-05-12
%
%  Copyright (C) 2021 IAU SOFA Board.  See notes at end.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eps0 = iauObl06(date1, date2)

global const

% Interval between fundamental date J2000.0 and given date (JC).
t = ((date1 - const.DJ00) + date2) / const.DJC;

% Mean obliquity.
eps0 = (84381.406     +...
       (-46.836769    +...
       ( -0.0001831   +...
       (  0.00200340  +...
       ( -0.000000576 +...
       ( -0.0000000434) * t) * t) * t) * t) * t) * const.DAS2R;

