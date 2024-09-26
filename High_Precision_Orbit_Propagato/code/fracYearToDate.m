function [year, month, day, hour, minute, second] = fracYearToDate(fractionOfYear)
year = floor(fractionOfYear);
fraction = fractionOfYear - year;

end