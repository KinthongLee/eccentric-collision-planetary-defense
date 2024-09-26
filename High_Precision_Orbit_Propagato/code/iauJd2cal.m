
function [iy, im, id, fd] = iauJd2cal(dj1, dj2)

% ----------------- find year and days of the year ---------------
jd     = dj1 + dj2;
temp   = jd - 2415019.5;
tu     = temp / 365.2425;
iy     = 1900 + floor( tu );

   % Calculate the number of leap years since 1900 up to the year before the current year
    leapyrs = 0;
    for year = 1900:(iy-1)
        if isLeapYear(year)
            leapyrs = leapyrs + 1;
        end
    end
days   = temp - ((iy-1900)*365.0 + leapyrs );
% ------------ check for case of beginning of a year -------------
if (days < 1.0)
    iy     = iy - 1;
    leapyrs = 0;
        for year = 1900:(iy-1)
            if isLeapYear(year)
                leapyrs = leapyrs + 1;
            end
        end
    days   = temp - ((iy-1900)*365.0 + leapyrs );
end
% ------------------- find remaining data  -----------------------
[im,id,hr,min,sec] = days2mdh(iy, days);

if (dj1 >= dj2)
    d1 = dj1;
    d2 = dj2;
else
    d1 = dj2;
    d2 = dj1;
end

d2 = d2 - 0.5;
f1 = mod(d1, 1);
f2 = mod(d2, 1);
fd = mod(f1 + f2, 1);

if (fd < 0.0)
    fd = fd + 1;
end

