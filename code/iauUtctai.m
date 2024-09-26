%  - - - - - - - - - -
%   i a u U t c t a i
%  - - - - - - - - - -
%
%  Time scale transformation:  Coordinated Universal Time, UTC, to
%  International Atomic Time, TAI.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards of Fundamental Astronomy) software collection.
%
%  Status:  canonical.
%
%  Given:
%     utc1,utc2     UTC as a 2-part quasi Julian Date (Notes 1-4)
%
%  Returned:
%     tai1,tai2     TAI as a 2-part Julian Date (Note 5)
%
%  Notes:
%  1) utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
%     convenient way between the two arguments, for example where utc1
%     is the Julian Day Number and utc2 is the fraction of a day.
%
%  2) JD cannot unambiguously represent UTC during a leap second unless
%     special measures are taken.  The convention in the present
%     function is that the JD day represents UTC days whether the
%     length is 86399, 86400 or 86401 SI seconds.
%
%  3) The warning status "dubious year" flags UTCs that predate the
%     introduction of the time scale and that are too far in the future
%     to be trusted.  See iauDat  for further details.
%
%  4) The function iauDtf2d converts from calendar date and time of day
%     into 2-part Julian Date, and in the case of UTC implements the
%     leap-second-ambiguity convention described above.
%
%  5) The returned TAI1,TAI2 are such that their sum is the TAI Julian
%     Date.
%
%  Called:
%     iauJd2cal    JD to Gregorian calendar
%     iauDat       delta(AT) = TAI-UTC
%     iauCal2jd    Gregorian calendar to JD
%
%  References:
%     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
%     IERS Technical Note No. 32, BKG (2004)
%     Explanatory Supplement to the Astronomical Almanac,
%     P. Kenneth Seidelmann (ed), University Science Books (1992)
%
%  This revision:  2010 September 10
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tai1, tai2] = iauUtctai(utc1, utc2)

constants

% Put the two parts of the UTC into big-first order.
big1 = (utc1 >= utc2);
if (big1)
    u1 = utc1;
    u2 = utc2;
else
    u1 = utc2;
    u2 = utc1;
end

% Get TAI-UTC now.
[iy, im, id, fd] = iauJd2cal(u1, u2);
dats = iauDat(iy, im, id, fd);

% Get TAI-UTC tomorrow.
[iyt, imt, idt, fdt] = iauJd2cal(u1+1.5, u2-fd);
datst = iauDat(iyt, imt, idt, fdt);

% If today ends in a leap second, scale the fraction into SI days.
ddat = datst - dats;

if (abs(ddat) > 0.5)
    fd = fd + fd * ddat / const.DAYSEC;
end

% Today's calendar date to 2-part JD.
[z1, z2] = iauCal2jd(iy, im, id);

% Assemble the TAI result, preserving the UTC split and order.
a2 = z1 - u1;
a2 = a2 + z2;
a2 = a2 + fd + dats / const.DAYSEC;

if (big1)
    tai1 = u1;
    tai2 = a2;
else
    tai1 = a2;
    tai2 = u1;
end

