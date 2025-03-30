function [mon,day,hr,minute,sec] = days2mdh(year, days)

% Helper function to check if the year is a leap year
function leap = isLeapYear(y)
    if mod(y, 4) == 0
        if mod(y, 100) == 0
            if mod(y, 400) == 0
                leap = true;
            else
                leap = false;
            end
        else
            leap = true;
        end
    else
        leap = false;
    end
end

% --------------- set up array of days in month  --------------
lmonth = [31 28 31 30 31 30 31 31 30 31 30 31];

if isLeapYear(year)
    lmonth(2) = 29;
end

dayofyr = floor(days);

% ----------------- find month and day of month ---------------
i = 1;
inttemp = 0;
while (dayofyr > inttemp + lmonth(i)) && (i < 12)
    inttemp = inttemp + lmonth(i);
    i = i + 1;
end

mon = i;
day = dayofyr - inttemp;

% ----------------- find hours minutes and seconds ------------
temp = (days - dayofyr) * 24.0;
hr = fix(temp);
temp = (temp - hr) * 60.0;
minute = fix(temp);
sec = (temp - minute) * 60.0;
end
