function monthStr = monthToString(monthNum)
% Converts a number representing a month (1-12) to a string with the
% corresponding month name.

% Define an array with the month names
months = {'January', 'February', 'March', 'April', 'May', 'June', ...
    'July', 'August', 'September', 'October', 'November', 'December'};

% Check that the input is a valid month number
if monthNum < 1 | monthNum > 12
    error('Invalid month number. Must be between 1 and 12.');
end

% Return the corresponding month name as a string
monthStr = months{monthNum};
end