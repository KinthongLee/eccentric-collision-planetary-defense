% clearvars

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
pathTosave = [fullfile(currentDir, 'output_result', 'different_PHA','launch_window','transfer_orbit_png'),'\'];
% ---------------------------------------------------------------------------------------


% --------------------------Read Data from Table-----------------------------------------
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
% ----------------------------------------------------------------------------------------


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





 % The code starts here
 % The gravitational parameter of the Sun
    mu = 1.32712440041279e+20;
 % Step Size: seconds
    Step = 86400; % s

 % Plot all the results from No.1 to No.32
 for p = 1 : 32
    if p ~= 9 && p ~= 14  % Exclude No.9 & No.14, due to unavailable Lambert Transfer
        target = num2str(data.BSP_file_name(p));
        name = cell2mat(data.Object(p));

        % Launch Date
        Mjd0_UTC = Mjday(data.Best_Launch_Year(p),data.Best_Launch_Month(p),data.Best_Launch_Date(p)) ; 
        [years0,months0, days0, fds0] = iauJd2cal( 2400000.5, Mjd0_UTC);
        [hours0, minutes0, secs0] = fd_to_hms(fds0);
        [years0, months0, days0, hours0, minutes0, secs0] = fix_seconds(years0, months0, days0, hours0, minutes0, secs0);
        secs0 = round(secs0,3);
        months0_string = monthToString(months0);
        str = sprintf(' %s %g , %g %g:%g:%g', months0_string, days0, years0, hours0, minutes0, secs0);
        et_start = cspice_str2et(str);
        
        % Impact date
        Mjd1_UTC =  Mjday(data.Best_Impact_Year(p),data.Best_Impact_Month(p),data.Best_Impact_Date(p));
        [years1,months1, days1, fds1] = iauJd2cal( 2400000.5, Mjd1_UTC);
        [hours1, minutes1, secs1] = fd_to_hms(fds1);
        [years1, months1, days1, hours1, minutes1, secs1] = fix_seconds(years1, months1, days1, hours1, minutes1, secs1);
        secs1 = round(secs1,3);
        months1_string = monthToString(months1);
        str = sprintf(' %s %g , %g %g:%g:%g', months1_string, days1, years1, hours1, minutes1, secs1);
        et_end = cspice_str2et(str);
        
        
        % Read Earth's position and velocity at Launch date
        [Y0, ~] = cspice_spkezr('EARTH', et_start, 'J2000', 'NONE', 'SUN');
        r_Earth0 = Y0(1:3).*1000; % m
        v_Earth0 = Y0(4:6).*1000; % m/s
        % Read PHA's position and velocity at Impact date
        [Y1, ~] = cspice_spkezr(target, et_end, 'J2000', 'NONE', 'SUN');
        r_Target1 = Y1(1:3).*1000; % m
        v_Target1 = Y1(4:6).*1000; % m/s
        
        % transfer duration: seconds
        delta_t = et_end - et_start; % s
        

        % Solving Lambert Problem：
        [v0a,v1a,~,~,~,~]=LambSol(r_Earth0,r_Target1,delta_t,mu,0);
        
        % Prepare to calculate orbit by HPOP
        transfer_orbit_a(1:3,1) = r_Earth0(1:3);
        transfer_orbit_a(4:6,1) = v0a;
        
        N_Step = delta_t / Step;

        % Transfer orbit (Lambert Transfer)
        transfer_orbit1 = Ephemeris_transfer(transfer_orbit_a, N_Step, Step,Mjd0_UTC);
    
        % A period orbit of Earth
        et_Earth = et_start : Step : et_start + 86400*366;
        [EARTH, ~] = cspice_spkezr('EARTH', et_Earth, 'J2000', 'NONE', 'SUN');
        EARTH = EARTH.*1000; % m

        % A period orbit of PHA
        et_Target = et_start - sqrt( (4*pi^2)/mu*data.a(p)^3 )*1.2 : Step : et_start  ;
        [TARGET, ~] = cspice_spkezr(target, et_Target, 'J2000', 'NONE', 'SUN');
        TARGET = TARGET.*1000; % m
        
        % Convert to astronomical units (AU)
        EARTH(1:3,:) = EARTH(1:3,:)./(1.496E11); % AU
        TARGET(1:3,:) = TARGET(1:3,:)./(1.496E11); % AU
        transfer_orbit1(:,2:4) = transfer_orbit1(:,2:4)./(1.496E11); % AU
    
    
        % Plot the graph
        gca = figure(1);
        % plot Earth's orbit
        plot3(EARTH(1,:),EARTH(2,:),EARTH(3,:),'blue','LineWidth', 4); m1 = 'Earth Orbit';
        hold on
        % plot PHA's orbit
        plot3(TARGET(1,:),TARGET(2,:),TARGET(3,:),'green','LineWidth', 4);m2 = 'Target Orbit';
        hold on
        % plot Transfer trajactory
        plot3(transfer_orbit1(:,2),transfer_orbit1(:,3),transfer_orbit1(:,4),'red','LineWidth', 4);
        m3 = 'Transfer Orbit';
        % lg = legend(m1,m2,m3);
        xlabel('x (AU)')
        ylabel('y (AU)')
        zlabel('z (AU)' )

        view(-45,45)
        t = title([num2str(p),'. ',name]);
        % To better fit the paper, the images have been restricted to a small size, but they can be adjusted as needed.
        set(gcf, 'Units', 'inches');
        set(gcf, 'Position', [1, 1, 3, 3])
        % Export graphic to file
        exportgraphics(gca, [pathTosave,num2str(p),'_transfer_orbit','.png']);
        close           
    end

 end

