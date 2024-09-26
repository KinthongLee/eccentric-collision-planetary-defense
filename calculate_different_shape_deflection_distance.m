clearvars
clc

% This code will Read the data from PHA_table

% This code will calculate all the distribution of deflection distance by
% both COG & BIP strategies through different types of 3D model of the
% asteroid

% The result of .mat file will be store in specific location.

% Unlikely other code, this code WILL NOT generate a new table, only save
% .mat file to futher code to plot the candle distribution of different
% shape

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
pathTosave = [fullfile(currentDir, 'output_result', 'different_shape','Apophis_orbit','matfile','different_shape'),'\'];
% ---------------------------------------------------------------------------------------



% ---------Read data and shape from PHA_table.xlsx-----------------------------------------
% Load shape
shape_data = readtable('PHA_shape');
% Load date
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
 
    % Start the SPMD block for parallel execution
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
 % -----------------------------------------------------------------------------------------





% Code starts here

% Number of Monte-Carlo samples modify if needed
 sample_size = 100000;



% Properties, modify if needed 
mass_Apophis = 2.7e10;
mass_rocket_0 = 9000; % kg
I_rocket = 320; % 比冲
% Use Apophis impact moment
p = 1;
target = num2str(data.BSP_file_name(p));
v_imp = [data.v_imp_X(p);data.v_imp_Y(p);data.v_imp_Z(p)];
v_ast = [data.v_ast_X(p);data.v_ast_Y(p);data.v_ast_Z(p)];
delta_v_a = data.Delta_v_a(p); % v_INF in the article
CA_distance = data.r_min_by_matlab(p); % m
dt1 = datetime(data.Best_Impact_Year(p),data.Best_Impact_Month(p),data.Best_Impact_Date(p),0,0,0);% Impact date
str1 = sprintf(' %s %g , %g %g:%g:%g', monthToString(data.Best_Impact_Month(p)),data.Best_Impact_Date(p),data.Best_Impact_Year(p),0,0,0); % Impact date
str2 = sprintf(' %s %g , %g %g:%g:%g', cell2mat(month_ini_string(p)),date_ini(p),year_ini(p),0,0,0); % Closest-approach
v_r = v_imp - v_ast;
r_ast = [data.x_ast_impact(p);data.y_ast_impact(p);data.z_ast_impact(p)];
mass_rocket = mass_rocket_0 * exp(-delta_v_a*1000 / (I_rocket * 9.80665)); % kg
et_impact = cspice_str2et(str1);
et_min_earth = cspice_str2et(str2);
delta_t_to_min_point = et_min_earth - et_impact;
Step   = 3600;   % [s] integration step size
N_Step = ceil(delta_t_to_min_point / Step) + 86400*2/3600;

et_impact = cspice_str2et(str1);
et_min_earth = cspice_str2et(str2);
delta_t_to_min_point = et_min_earth - et_impact;
Step   = 3600;   % [s] integration step size

% pp is the row from the table PHA_shape
for pp = 2 : 2
model = [cell2mat(shape_data.Model(pp)),'.obj'];
[vertices, faces] = readObj(model);

% Show progress in the parfor loop
    parfor_progress(sample_size);

    % starts the loop
    tic;
    norm_delta_v = zeros(sample_size,1);
    norm_delta_v_COM = zeros(sample_size,1);
    delta_r_min_BIP = zeros(sample_size,1);
    delta_r_min_COG = zeros(sample_size,1);
    delta_v = zeros(sample_size,3);
    delta_v_COM = zeros(sample_size,3);
    angles = zeros(sample_size,3);
parfor i = 1 : sample_size 
    beta = 3.61 + unifrnd(-0.25,0.19); % generate beta coefficien from a skewness distribution as mentioned in article
    [dv_BIP,dv_COG,ang] = get_delta_v_from_momentum(v_imp,v_ast,mass_rocket,vertices,faces,beta);
    delta_v(i,:) = dv_BIP;
    delta_v_COM(i,:) = dv_COG';
    angles(i,:) = ang;
    norm_delta_v(i,1) = norm(dv_BIP);
    norm_delta_v_COM(i,1) = norm(dv_COG);
    delta_r_min_BIP(i,1) = Propagation(target,dv_BIP',CA_distance,dt1,Step,N_Step);
    delta_r_min_COG(i,1) = Propagation(target,dv_COG,CA_distance,dt1,Step,N_Step);
    parfor_progress;
end
toc;
    vt = v_ast / norm(v_ast);
    vh = cross(r_ast,vt) / norm(cross(r_ast,vt) );
    vn = cross(vh,vt) / norm(cross(vh,vt));
    delta_v_BIP_avg = [mean(delta_v(:,1));mean(delta_v(:,2));mean(delta_v(:,3))];
    delta_v_COG_avg = [mean(delta_v_COM(:,1));mean(delta_v_COM(:,2));mean(delta_v_COM(:,3))];
    % project delta v to t,n,h direction
    delta_v_BIP_final = [dot(delta_v_BIP_avg,vt),dot(delta_v_BIP_avg,vn),dot(delta_v_BIP_avg,vh)];
    delta_v_COG_final = [dot(delta_v_COG_avg,vt),dot(delta_v_COG_avg,vn),dot(delta_v_COG_avg,vh)];

 message = sprintf('calculate for %s is Done', cell2mat(shape_data.Model(pp)));
    disp(message);


    filename = [pathTosave,cell2mat(shape_data.Model(pp)), '_distribution.mat'];  % Append .mat extension
    % 保存.mat文件
    save(filename, 'delta_r_min_COG','delta_r_min_BIP','norm_delta_v_COM','norm_delta_v','angles','delta_v','delta_v_COM','delta_v_BIP_final','delta_v_COG_final');


end
