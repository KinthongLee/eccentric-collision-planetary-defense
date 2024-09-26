function numDays = daysInMonth(dt)
% Return the number of days in the given month for the given datetime object

year = dt.Year;
month = dt.Month;

if month == 2
    % Check for leap year
    if mod(year, 4) == 0 && (mod(year, 100) ~= 0 || mod(year, 400) == 0)
        numDays = 29;
    else
        numDays = 28;
    end
elseif any(month == [4, 6, 9, 11])
    numDays = 30;
else
    numDays = 31;
end