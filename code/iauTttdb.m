%  - - - - - - - - -
%   i a u T t t d b
%  - - - - - - - - -
%
%  Time scale transformation:  Terrestrial Time, TT, to Barycentric
%  Dynamical Time, TDB.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards of Fundamental Astronomy) software collection.
%
%  Status:  canonical.
%
%  Given:
%     tt1,tt2        TT as a 2-part Julian Date
%     dtr            TDB-TT in seconds
%
%  Returned:
%     tdb1,tdb2      TDB as a 2-part Julian Date
%
%  Notes:
%  1) tt1+tt2 is Julian Date, apportioned in any convenient way between
%     the two arguments, for example where tt1 is the Julian Day Number
%     and tt2 is the fraction of a day.  The returned tdb1,tdb2 follow
%     suit.
%
%  2) The argument dtr represents the quasi-periodic component of the
%     GR transformation between TT and TCB.  It is dependent upon the
%     adopted solar-system ephemeris, and can be obtained by numerical
%     integration, by interrogating a precomputed time ephemeris or by
%     evaluating a model such as that implemented in the SOFA function
%     iauDtdb.   The quantity is dominated by an annual term of 1.7 ms
%     amplitude.
%
%  3) TDB is essentially the same as Teph, the time argument for the JPL
%     solar system ephemerides.
%
%  References:
%     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
%     IERS Technical Note No. 32, BKG (2004)
%     IAU 2006 Resolution 3
%
%  This revision:  2011 May 14
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tdb1, tdb2] = iauTttdb(tt1, tt2, dtr)

constants

% Result, safeguarding precision.
dtrd = dtr / DAYSEC;

if (tt1 > tt2)
    tdb1 = tt1;
    tdb2 = tt2 + dtrd;
else
    tdb1 = tt1 + dtrd;
    tdb2 = tt2;
end

