function leap = isLeapYear(year)
    if mod(year, 4) == 0
        if mod(year, 100) == 0
            if mod(year, 400) == 0
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