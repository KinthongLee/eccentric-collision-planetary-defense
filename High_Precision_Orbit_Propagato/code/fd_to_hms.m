function [hours, minutes, seconds] = fd_to_hms(fd)
% Convert a fractional day to hours, minutes, and seconds
fd = fd - floor(fd); % Extract fractional part
hours = floor(fd * 24);
minutes = floor(mod(fd * 1440, 60));
seconds = mod(fd * 86400, 60);
end