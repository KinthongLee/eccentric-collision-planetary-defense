clc
clearvars

% This code will Read the data from PHA_table, spefically Columns A to J
% This code will calculate the Closest-approach with Earth and the relative
% error with SPICE. The result will be export as a new table temporary_result.xlsx
% If the result is correct, please manually copy and paste the result into PHA_table.xlsx.
% MAKE SURE TO DOWNLOAD THE SPECIFIC ASTEROID ".bsp" FILE AND STORE IT TO
% mice/kernel/ !! OTHERWISE IT WILL FAILED TO GET DATA OF ASTEROID FROM
% SPICE

% ---------------------------  Add Path  -----------------------------------------------
addpath('mice\')
addpath('mice\lib\')
addpath('code\')
addpath('mice\src\mice\')

% add path to kernel  
% Specify the directory where your .bsp files are located
currentDir = fileparts(which('plot_transfer_orbit.m'));
pathTokernel = fullfile(currentDir, 'mice', 'kernel');

% ---------------------------------------------------------------------------------------

close all
profile on
tic

% ---------Read data from PHA_table.xlsx-----------------------------------------
data = readtable('PHA_table');
% Extract year, month, and day from the 'Close_Approach_CA_Date' column
dateStrings = data.Close_Approach_CA_Date;
datePattern = '(\d{4})-(\w{3})-(\d{2})';
tokens = regexp(dateStrings, datePattern, 'tokens');

% Convert tokens to a matrix (each row corresponds to a match)
tokensMatrix = vertcat(tokens{:});
tokensMatrix = vertcat(tokensMatrix{:});

% Extract and convert each component
year_ini = str2double(tokensMatrix(:,1));
date_ini = str2double(tokensMatrix(:,3));

% Map the month abbreviations to month numbers and full names
monthAbbreviations = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', ...
                      'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
monthNumbers = 1:12;
monthFullNames = {'January', 'February', 'March', 'April', 'May', 'June', ...
                  'July', 'August', 'September', 'October', 'November', 'December'};

% Create containers.Map to map month abbreviations to month numbers and full names
monthNumMap = containers.Map(monthAbbreviations, monthNumbers);
monthNameMap = containers.Map(monthAbbreviations, monthFullNames);

% Convert month abbreviations to numbers and full names
month_ini_num = cell2mat(values(monthNumMap, tokensMatrix(:,2)));
month_ini_string = values(monthNameMap, tokensMatrix(:,2));
% -----------------------------------------------------------------------


%-----------------------------Load Kernels------------------------------------------------
 % List all .bsp files in the directory
bspFiles = dir(fullfile(pathTokernel, '*.bsp'));

% List all .txt files in the directory
txtFiles = dir(fullfile(pathTokernel, '*.txt'));

% Combine the .bsp and .txt files into a single array
allFiles = [bspFiles; txtFiles];

% Iterate over all files to load them
for k = 1:length(allFiles)
    % Construct the full path for each file
    filePath = fullfile(pathTokernel, allFiles(k).name);
    
    % Load the file using cspice_furnsh
    cspice_furnsh(filePath);
end
 % -----------------------------------------------------------------------------------------

% ---------------STep and total steps ---------------------------------
Step   = 3600;   % [s] integration step size
N_Step = 24*366*10; % number of integration steps 
% -----------------------------------------------------------------------------------------

% Code starts here
for p = 33 : 33
% Initial States 
    target = num2str(data.BSP_file_name(p));
    year0 = year_ini(p)-10;
    month0 = month_ini_num(p);
    day0 = date_ini(p);
    hour0 = 0;
    minute0 = 0;
    sec0 = 0;
    
    % Get starting position & velocity from SPICE:
    Mjd0_UTC = Mjday(year0, month0, day0, hour0, minute0, sec0);
    month0_string = monthToString(month0);
    str = sprintf(' %s %g , %g %g:%g:%g', month0_string, day0, year0, hour0, minute0, sec0);
    et_start = cspice_str2et(str);
    [Y0, ~] = cspice_spkezr(target, et_start, 'J2000', 'NONE', 'SUN');
    Y0 = Y0.*1000; % m/s
    Y0(4:6) = Y0(4:6);
    
    %  Propagation
    disp('Start Propagation')
     Eph_eci = Ephemeris(Y0, N_Step, Step,Mjd0_UTC);

    [year,month, day, fd] = iauJd2cal( 2400000.5, Mjd0_UTC+Eph_eci( size(Eph_eci,1) ,1)/86400);
    [hour, minute, sec] = fd_to_hms(fd);
    [year, month, day, hour, minute, sec] = fix_seconds(year, month, day, hour, minute, sec);
    month = monthToString(month);
    sec = round(sec,3);
    str_stop = sprintf(' %s %g , %g %g:%g:%g', month, day, year, hour, minute, sec);
    et_stop = cspice_str2et(str_stop);
    et = et_start : Step : et_stop;
    % Sometimes if will miss the last et_stop
    if length(et) < length(Eph_eci) 
        et(1,length(et)+1) = et_stop;
    end
    
    
    [Y_EARTH, ~] = cspice_spkezr('EARTH', et, 'J2000', 'NONE', 'SUN');
    Y_EARTH = Y_EARTH.*1000;
    
    
    % 


    
    % Calculate Standard Error of distance with Earth and Time
    diff_r = zeros(min(length(Eph_eci),length(Y_EARTH)),1);
    
    for i = 1 : min(length(Eph_eci),length(Y_EARTH))
        r = Eph_eci(i,2:4)' - Y_EARTH(1:3,i);
        diff_r(i) = norm( r );
    end
    % Get the minimum distance with Earth
    r_min = min(diff_r);
    [I,J] = find(diff_r == r_min);
    % [futureYear, futureMonth, futureDay, futureHour, futureMinute, futureSecond] = futureDate(year0,month0,day0,hour0,minute0,sec0,I*Step);
    % futureMonth = monthToString(futureMonth);
    
    % Get the minimum distance from the DE441 as a comparison
    [Y_ast, ~] = cspice_spkezr(target, et, 'J2000', 'NONE', 'SUN');
    Y_ast = Y_ast.*1000;
    diff_r_DE441 = zeros( length(Y_ast),1);
    for i = 1 : min(length(Y_ast),length(Y_EARTH))
        r = Y_ast(1:3,i) - Y_EARTH(1:3,i);
        diff_r_DE441(i) = norm( r );
    end
    r_min_DE441 = min(diff_r_DE441);

    % Relative Error
    relative_error_min_r = abs(r_min - r_min_DE441 ) / r_min_DE441 * 100;
   
    message = fprintf('HPOP for %i th asteroid: %s is Done', p,cell2mat(data.Object(p)));
    disp(message);
    % Write data to table
    data.r_min_by_matlab(p) = r_min;
    data.r_min_relative_error(p) = relative_error_min_r;
end

% Export table
writetable(data,'temporary_result.xlsx')


