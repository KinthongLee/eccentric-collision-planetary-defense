function [year, month, day] = fix_date(year, month, day)


try
    dt = datetime(year, month, day);
    dt = dt + calmonths(1); % increment month
    [year, month, day] = ymd(dt);
catch
    % handle invalid date
    error('Invalid date');
end
end