function [year, month, day, hour, minute, second] = fix_seconds(year, month, day, hour, minute, second)

if second > 59.999
    minute = minute + 1; % fix rounds toward zero
    second = 0;
    
    if minute >= 60
        hour = hour + 1;
        minute = 0;
        
        if hour >= 24
            day = day + 1;
            hour = 0;
            [year, month, day] = fix_date(year, month, day);
        end
    end
end

function [year, month, day] = fix_date(year, month, day)


days_in_month = [31 28 31 30 31 30 31 31 30 31 30 31];

if mod(year, 4) == 0 && (mod(year, 100) ~= 0 || mod(year, 400) == 0)
    days_in_month(2) = 29; % Leap year
end

while day > days_in_month(month)
    day = day - days_in_month(month);
    month = month + 1;


end