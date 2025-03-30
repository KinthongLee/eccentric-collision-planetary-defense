%  - - - - - - - - -
%   i a u T a i t t
%  - - - - - - - - -
%
%  Time scale transformation:  International Atomic Time, TAI, to
%  Terrestrial Time, TT.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards of Fundamental Astronomy) software collection.
%
%  Status:  canonical.
%
%  Given:
%     tai1,tai2  double    TAI as a 2-part Julian Date
%
%  Returned:
%     tt1,tt2    double    TT as a 2-part Julian Date
%
%  Note:
%     tai1+tai2 is Julian Date, apportioned in any convenient way
%     between the two arguments, for example where tai1 is the Julian
%     Day Number and tai2 is the fraction of a day.  The returned
%     tt1,tt2 follow suit.
%
%  References:
%     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
%     IERS Technical Note No. 32, BKG (2004)
%
%     Explanatory Supplement to the Astronomical Almanac,
%     P. Kenneth Seidelmann (ed), University Science Books (1992)
%
%  This revision:  2011 May 14
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tt1, tt2] = iauTaitt(tai1, tai2)

constants

% TT minus TAI (days).
dtat = const.TTMTAI/const.DAYSEC;

% Result, safeguarding precision.
if (tai1 > tai2)
    tt1 = tai1;
    tt2 = tai2 + dtat;
else
    tt1 = tai1 + dtat;
    tt2 = tai2;
end

