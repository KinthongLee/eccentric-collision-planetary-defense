clearvars
clc
% 
% This code will Read the data from PHA_table
% This code will calculate all possible launch window by solving Lambert
% problem. 
% The result of .mat file will be store in specific location.
% The result will also be export as a new table temporary_result.xlsx,
% spefically Columns M to AO
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

% add path to save
pathTosave = [fullfile(currentDir, 'output_result', 'different_PHA','launch_window','matfile'),'\'];
% ---------------------------------------------------------------------------------------


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
 
    % Start the SPMD block for parallel execution （parfor loop）
    spmd
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
    end
    % Load kernel for non parfor loop
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




% The code starts here
% Gravitational constant of Sun
mu = 1.32712440041279e+20;   

for i = 33 : 33
    % Read BSP name of the asteroid
    target = num2str(data.BSP_file_name(i)); 
    % Mass of Asteroid
    mass_ast = 2.7e10; %kg
    % beta coefficient
    beta = 3.61;
    % Closest-Approach date with Earth
    str_min_earth = sprintf(' %s %g , %g %g:%g:%g', cell2mat(month_ini_string(i)),date_ini(i),year_ini(i),0,0,0);
    mjdate_for_min_earth = Mjday(year_ini(i),month_ini_num(i),date_ini(i),0,0,0);
    
    % Minimum distance with Earth
    min_distance_to_Earth = data.r_min_by_matlab(i);
    
    % mass of the impactor
    mass_rocket_0 = 9000; % kg
    % Specific impulse of rocket
    I_rocket = 320; % s
    % Payload capacity
    C3 = 50; %km^2 /s^2 
    
    % Specific size of the launch window
    % Begining of launch date
    % Here make the begining of launch date by 10 years before
    % closest-approach year. Modify if needed
    launch_date = Mjday(year_ini(i)-10,1,1,0,0,0); % Convert to MJD
    % Step size of each launch window
    % Modiy if needed
    launch_Step = 1; % day
    % Total Step, here is 10 years, modify if needed
    launch_NStep = 365*8+366*2;
    % minimum transfer days
    minimum_transfer_days = 0; % day
    % Last date of the launch window
    % here make the year of the closest approach as the end of thte launch
    % window, modify if needed
    last_transfer_date = Mjday(year_ini(i),1,1,0,0,0);
    % transfer days, keep it unchanged
    transfer_days = last_transfer_date - launch_date; % day
    % step size for the High Precision Orbit Propagrator
    Step   = 3600;   % seconds
    
    
    % 设置合理的发射&到达时机：
    total_launch_Step = launch_NStep / launch_Step; % 需要测试的总迭代数目
    % 内容为 way = 0代表逆时针
    possible_launch = NaT(total_launch_Step,transfer_days); 
    possible_reach = NaT(total_launch_Step,transfer_days);  
    delta_r_min = zeros(total_launch_Step,transfer_days);
    
    
    % Show the parfor progress
    parfor_progress(total_launch_Step);
    
    % Start calculate
    tic;
    parfor delta_days_launch = 1  : total_launch_Step
    
        for days_reach = 1 : transfer_days
    
            % begining date
            Mjd0_UTC = launch_date + launch_Step * delta_days_launch - 1;
    
            % last date
            Mjd1_UTC = launch_date + days_reach;
            
    
            if Mjd1_UTC - minimum_transfer_days >  Mjd0_UTC  && Mjd1_UTC < mjdate_for_min_earth 
            
                % Launch date
                [years0,months0, days0, fds0] = iauJd2cal( 2400000.5, Mjd0_UTC);
                [hours0, minutes0, secs0] = fd_to_hms(fds0);
                [years0, months0, days0, hours0, minutes0, secs0] = fix_seconds(years0, months0, days0, hours0, minutes0, secs0);
                secs0 = round(secs0,3);
                months0_string = monthToString(months0);
                str = sprintf(' %s %g , %g %g:%g:%g', months0_string, days0, years0, hours0, minutes0, secs0);
                et_start = cspice_str2et(str);
                dt0 = datetime(years0, months0, days0, hours0, minutes0, secs0);
                
                
                % Impact date
                [years1,months1, days1, fds1] = iauJd2cal( 2400000.5, Mjd1_UTC);
                [hours1, minutes1, secs1] = fd_to_hms(fds1);
                [years1, months1, days1, hours1, minutes1, secs1] = fix_seconds(years1, months1, days1, hours1, minutes1, secs1);
                secs1 = round(secs1,3);
                dt1 = datetime(years1, months1, days1, hours1, minutes1, secs1);
                months1_string = monthToString(months1);
                str = sprintf(' %s %g , %g %g:%g:%g', months1_string, days1, years1, hours1, minutes1, secs1);
                et_end = cspice_str2et(str);
                
                
                % Read the state of Earth at launch date and Asteroid at
                % impact date
                [Y0, ~] = cspice_spkezr('EARTH', et_start, 'J2000', 'NONE', 'SUN');
                r_Earth0 = Y0(1:3).*1000; % m
                v_Earth0 = Y0(4:6).*1000; % m/s
                
                [Y1, ~] = cspice_spkezr(target, et_end, 'J2000', 'NONE', 'SUN');
                r_Target1 = Y1(1:3).*1000; % m
                v_Target1 = Y1(4:6).*1000; % m/s
                
                % transfer duration, s
                delta_t = et_end - et_start;
                
               
                % Solve the Lambert problem
                [v0a,v1a,~,~,~,~]=LambSol(r_Earth0,r_Target1,delta_t,mu,0);
                delta_v_a = norm( v0a - v_Earth0 )/1000; %km/s
    
                % Total Steps for the propagator
                et_min_earth = cspice_str2et(str_min_earth);
                delta_t_to_min_point = et_min_earth - et_end;
                N_Step = ceil(delta_t_to_min_point / Step) + 86400 / Step; %有时候et的运算会少一天
                
    
                % Does the rocket capacity meet the requirements?
                if delta_v_a^2 <= C3
                    mass_rocket = mass_rocket_0 * exp(-delta_v_a*1000 / (I_rocket * 9.80665)); % kg
                    % calculate interception angle (Alpha):
                    v_r = (v1a - v_Target1);
                    % calculate delta_v, hit through Center of Geometry
                    delta_v_to_Apophis = beta* mass_rocket / (mass_rocket+mass_ast) * (v_r);
                    % Calculate the changes in CA distance through High
                    % Precision Orbit Propagator
                    delta_r_min(delta_days_launch,days_reach) = Propagation(target,delta_v_to_Apophis,min_distance_to_Earth,dt1,Step,N_Step);
                    % Write the possible launch and impact date
                    possible_launch(delta_days_launch,days_reach) = dt0; % 逆时针
                    possible_reach(delta_days_launch,days_reach) = dt1;
                end
    
            end
        end
        parfor_progress;
    end
    toc;
    message = sprintf('calculate for %s is Done', cell2mat(data.Object(i)));
    disp(message);
    filename = [pathTosave, cell2mat(data.Object(i)), '.mat']; 

    % save the .mat file
    save(filename, 'possible_launch', 'possible_reach', 'delta_r_min');
    
    % Find the best Lambert transfer
    [largestElement, index] = max(delta_r_min(:));
    [row, col] = ind2sub(size(delta_r_min), index);
    

    % Write the best lambert transfer to table
    data.Best_Launch_Year(i) = year(possible_launch(row, col));
    data.Best_Launch_Month(i) = month(possible_launch(row, col));
    data.Best_Launch_Date(i) = day(possible_launch(row, col));
    data.Best_Impact_Year(i) = year(possible_reach(row, col));
    data.Best_Impact_Month(i) = month(possible_reach(row, col));
    data.Best_Impact_Date(i) = day(possible_reach(row, col));
    data.delta_r_min(i) = delta_r_min(row, col);

    % Write the data to table
    if data.delta_r_min(i) ~= 0
         Mjd0_UTC = Mjday(year(possible_launch(row, col)),month(possible_launch(row, col)),day(possible_launch(row, col))) ; 
         Mjd1_UTC =  Mjday(year(possible_reach(row, col)),month(possible_reach(row, col)),day(possible_reach(row, col)));
         % Launch date
        [years0,months0, days0, fds0] = iauJd2cal( 2400000.5, Mjd0_UTC);
        [hours0, minutes0, secs0] = fd_to_hms(fds0);
        [years0, months0, days0, hours0, minutes0, secs0] = fix_seconds(years0, months0, days0, hours0, minutes0, secs0);
        secs0 = round(secs0,3);
        months0_string = monthToString(months0);
        str = sprintf(' %s %g , %g %g:%g:%g', months0_string, days0, years0, hours0, minutes0, secs0);
        et_start = cspice_str2et(str);
        
        % Impact date
        [years1,months1, days1, fds1] = iauJd2cal( 2400000.5, Mjd1_UTC);
        [hours1, minutes1, secs1] = fd_to_hms(fds1);
        [years1, months1, days1, hours1, minutes1, secs1] = fix_seconds(years1, months1, days1, hours1, minutes1, secs1);
        secs1 = round(secs1,3);
        months1_string = monthToString(months1);
        str = sprintf(' %s %g , %g %g:%g:%g', months1_string, days1, years1, hours1, minutes1, secs1);
        et_end = cspice_str2et(str);
        
        
        % Read the state of Earth and Asteroid at Launch and Impact date
        [Y0, ~] = cspice_spkezr('EARTH', et_start, 'J2000', 'NONE', 'SUN');
        r_Earth0 = Y0(1:3).*1000; % m
        v_Earth0 = Y0(4:6).*1000; % m/s
        
        [Y1, ~] = cspice_spkezr(target, et_end, 'J2000', 'NONE', 'SUN');
        r_Target1 = Y1(1:3).*1000; % m
        v_Target1 = Y1(4:6).*1000; % m/s
        
        % transfer duration, s
        delta_t = et_end - et_start;
        
        % Solve the Lambert prob
        [v0a,v1a,~,~,~,~]=LambSol(r_Earth0,r_Target1,delta_t,mu,0);
        Delta_v_a = norm(v0a - v_Earth0)./1000;
    
        mass_rocket = mass_rocket_0 * exp(-Delta_v_a*1000 / (I_rocket * 9.80665));
    
        [a,e,i_orbit,omega,w,f] = input_r_v(r_Target1,v_Target1,mu);

        
        % Calculate the true anomaly at the closest-approach
        et_min_earth = cspice_str2et(str_min_earth);
        [Y_min_Earth, ~] = cspice_spkezr(target, et_min_earth, 'J2000', 'NONE', 'SUN');
        r_min_Earth = Y_min_Earth(1:3).*1000; % m
        v_min_Earth = Y_min_Earth(4:6).*1000; % m/s
    
        [~,~,~,~,~,f_MOID] = input_r_v(r_min_Earth,v_min_Earth,mu);

        
        delta_t_to_min_Earth = et_min_earth - et_end;
       
        data.f_MOID(i) = f_MOID;
        data.delta_t(i) = delta_t_to_min_Earth;
        data.delta_t_by_T(i) = delta_t_to_min_Earth / ( 2*pi*sqrt( a^(3) / mu ) );
        data.x_ast_impact(i) = r_Target1(1);
        data.y_ast_impact(i) = r_Target1(2);
        data.z_ast_impact(i) = r_Target1(3);
        data.v_imp_X(i) = v1a(1);
        data.v_imp_Y(i) = v1a(2);
        data.v_imp_Z(i) = v1a(3);
        data.v_ast_X(i) = v_Target1(1);
        data.v_ast_Y(i) = v_Target1(2);
        data.v_ast_Z(i) = v_Target1(3);
        data.Relative_U(i) = norm(v1a - v_Target1)./1000;
        data.Alpha(i) = acosd(dot(v1a,v_Target1) / norm(v1a) / norm(v_Target1));
        data.Delta_v_a(i) = Delta_v_a;
        data.Imp_Mass(i) = mass_rocket;
        data.a(i) = a;
        data.e(i) = e;
        data.i(i) = i_orbit;
        data.omega(i) = omega;
        data.w(i) = w;
        data.f(i) = f;
        
    end
end

    writetable(data,'temporary_result.xlsx')