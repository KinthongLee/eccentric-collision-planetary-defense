function [futureYear, futureMonth, futureDay, futureHour, futureMinute, futureSecond] = futureDate(startYear, startMonth, startDate, startHour, startMinute, startSecond, gapInSeconds)
% Calculate the future date from a start date and a gap in seconds
% Input: startYear - integer
%        startMonth - string
%        startDate - integer
%        startHour - integer
%        startMinute - integer
%        startSecond - integer
%        gapInSeconds - integer
% Output: futureYear - integer
%         futureMonth - string
%         futureDate - integer
%         futureHour - integer
%         futureMinute - integer
%         futureSecond - integer

% Convert the month string to a number
monthNum = month(datenum(startYear, startMonth, startDate));

% Convert the start date and time to a MATLAB datenum value
startDateNum = datenum(startYear, monthNum, startDate, startHour, startMinute, startSecond);

% Add the gap in seconds to the start date and time
futureDateNum = startDateNum + gapInSeconds/86400;

% Convert the future date and time back to year, month, date, hour, minute, and second
[futureYear, futureMonth, futureDay, futureHour, futureMinute, futureSecond] = datevec(futureDateNum);

% Convert the month number back to a string
futureMonth = month(futureDateNum, 'long');
end