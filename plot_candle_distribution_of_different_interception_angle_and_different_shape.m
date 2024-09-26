% ---------------------------  Add Path  -----------------------------------------------
currentDir = fileparts(which('plot_launch_window.m'));

% add path to save
pathTosave = [fullfile(currentDir, 'output_result', 'different_shape','Apophis_orbit','pngfile'),'\'];

% add path of the .mat file stored
% different degree
pathOfmatfile_degree = [fullfile(currentDir, 'output_result', 'different_shape','Apophis_orbit','matfile','different_degree'),'\'];
% different shape
pathOfmatfile_shape = [fullfile(currentDir, 'output_result', 'different_shape','Apophis_orbit','matfile','different_shape'),'\'];
% ---------------------------------------------------------------------------------------


% --------------------------Read Data from Table-----------------------------------------
data = readtable('NEO_different_shape');
shape_data = readtable('NEO_different_shape');
% ---------------------------------------------------------------------------------------




% % Code starts here
boxplot_degree = zeros(100000,52);
boxplot_shape = zeros(100000,24);
% ---------------------------different interception angle----------------------------------------------- 

% Read Data
% Read negative, 0 means dot, for example 3205 means 32.5
angle_n = [40,3705,35,3205,30,2705,25,2205,20,1705,15,1205,10,705,5,205]; % Modify if  necessary, make sure the specific matfile ady ran and stored
for i = 1:length(angle_n)
    name = ['Apophis_Model_n_',num2str(angle_n(i))];
    filename = [name,'_distribution.mat']; % Change to your desired file name
    % Full path for the file
    fullFilePath = fullfile(pathOfmatfile_degree, filename);

    % Check if the file exists
    if isfile(fullFilePath)
        % File exists, read it using readtable
        load(fullFilePath);
        disp('File read successfully.');
    else
        % File does not exist
        error([ filename, ' not found in ',pathOfmatfile_degree]);
    end
    boxplot_degree(:,2*i-1) = delta_r_min_COG./1000;
    boxplot_degree(:,2*i) = delta_r_min_BIP./1000;
end


% read zero degree
name = 'Apophis_Model_0';
    filename = [name,'_distribution.mat']; % Change to your desired file name
    % Full path for the file
    fullFilePath = fullfile(pathOfmatfile_degree, filename);

    % Check if the file exists
    if isfile(fullFilePath)
        % File exists, read it using readtable
        load(fullFilePath);
        disp('File read successfully.');
    else
        % File does not exist
         error([ filename, ' not found in ',pathOfmatfile_degree]);
    end
    boxplot_degree(:,2*(i+1)-1) = delta_r_min_COG./1000;
    boxplot_degree(:,2*(i+1)) = delta_r_min_BIP./1000;




% read positive, 0 means dot, for example 3205 means 32.5
angle_p = [205,5,705,10,1205,15,1705,20,2205,25]; % Modify if  necessary, make sure the specific matfile ady ran and stored
for i = 1:length(angle_p)
    name = ['Apophis_Model_p_',num2str(angle_p(i))];
    filename = [name,'_distribution.mat']; % Change to your desired file name
    % Full path for the file
    fullFilePath = fullfile(pathOfmatfile_degree, filename);

    % Check if the file exists
    if isfile(fullFilePath)
        % File exists, read it using readtable
        load(fullFilePath);
        disp('File read successfully.');
    else
        % File does not exist
         error([ filename, ' not found in ',pathOfmatfile_degree]);
    end
    boxplot_degree(:,(length(angle_n)+1)*2 + 2*i-1) = delta_r_min_COG./1000;
    boxplot_degree(:,(length(angle_n)+1)*2 + 2*i) = delta_r_min_BIP./1000;

end


% plot candles

xpositions = [-40.5, -39.5,-38,-37, -35.5, -34.5,-33,-32 -30.5, -29.5,-28,-27 -25.5,-24.5,-23,-22,-20.5,-19.5,-18,-17,-15.5,-14.5,-13,-12,-10.5,-9.5,...
    -8,-7,-5.5,-4.5,-3,-2,-0.5,0.5,2,3,4.5,5.5,7,8,9.5,10.5,12,13,14.5,15.5,17,18,19.5,20.5,22,23,24.5,25.5];

% Shift the 0 point of the x-axis to the original angle of Apophis interception mission (5.22 deg).
xpositions = xpositions + 5.22; 
boxplot(boxplot_degree,'Positions',xpositions)

% Remove all x-axis ticks
set(gca, 'XTick', []);
gca1 = figure(1);
% Set specific x-axis ticks and labels
xticks = [-40,-37.5,-35,-32.5,-30,-27.5,-25,-22.5,-20,-17.5,-15,-12.5,-10,-7.5,-5,-2.5,0,2.5,5,7.5,10,12.5,15,17.5,20,22.5,25] + 5.22;  % Define specific positions where you want ticks
xticklabel = strings(length(xticks),1);
for i = 1 : length(xticks)
    xticklabel(i) = num2str(xticks(i));
end
set(gca, 'XTick', xticks, 'XTickLabel', xticklabel);
ylabel('Deflection Distance (km)')
xlabel('Interception Angle (deg)')
title('Distribution of COG and BIP Strategies Pair (Sequential) under Different Interception Angles')

hold on;

% Highlight the original angle of Apophis interception mission (5.22 deg).
h1 = fill([4.22, 6.22, 6.22, 4.22], ...
     [-8000, -8000, 4050, 4050], [0.9, 0.9, 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5); 

lgd = legend(h1, {'Initial Angle'}, 'Location', 'BestOutside');
ylim([-8000 4050])
set(gcf, 'Position', [100, 100, 1600, 700]);
set(lgd, 'Position', [0.62,0.3,0.1,0.05]);
set(gca, 'FontSize', 14)

% Export graph
exportgraphics(gca1, [pathTosave,'Candles_distribution_of_different_interception_angle','.png']);



% Plot graph of Gain(%) BIP to COG strategy
% Calculate data ready to be used
mean_COG_degree = zeros(1,length(boxplot_degree(1,:))/2);
std_COG_degree = zeros(1,length(boxplot_degree(1,:))/2);
std_BIP_degree = zeros(1,length(boxplot_degree(1,:))/2);
mean_BIP_degree = zeros(1,length(boxplot_degree(1,:))/2);
Gain_degree = zeros(1,length(boxplot_degree(1,:))/2);
num_less_zero_COG = zeros(1,length(boxplot_degree(1,:))/2);
num_less_zero_BIP = zeros(1,length(boxplot_degree(1,:))/2);


for i = 1 : length(boxplot_degree(1,:))/2
    mean_COG_degree(i) = mean(boxplot_degree(:,2*i-1));
    mean_BIP_degree(i) = mean(boxplot_degree(:,2*i));
    std_COG_degree(i) = std(boxplot_degree(:,2*i-1));
    std_BIP_degree(i) = std(boxplot_degree(:,2*i));
    Gain_degree(i) = ( mean_BIP_degree(i) - mean_COG_degree(i) ) / abs(mean_COG_degree(i)) *100 ;
    num_less_zero_COG(1,i) = sum(boxplot_degree(:,2*i-1) < 0);
    num_less_zero_BIP(1,i) = sum(boxplot_degree(:,2*i) < 0);
end

xticks = [-40,-37.5,-35,-32.5,-30,-27.5,-25,-22.5,-20,-17.5,-15,-12.5,-10,-7.5,-5,-2.5,0,2.5,5,7.5,10,12.5,15,17.5,20,22.5,25] + 5.22;  % Define specific positions where you want ticks


% Plot scatter graph of GAIN(%)
gca2 = figure(2);
yyaxis right
plot1 = scatter(xticks(5:23),Gain_degree(5:23),'black'); m1 = 'GAIN(%)';
ylabel('GAIN(%)')
hold on
set(gca, 'YColor', 'k')  % Set the right y-axis to black
yyaxis left
plot2 = scatter(xticks(5:23),mean_BIP_degree(5:23),'blue'); m2 = 'BIP_{mean}';
ylabel('Deflection Distance (km)')
hold on
plot3 = scatter(xticks(5:23),mean_COG_degree(5:23),'red'); m3 = 'COG_{mean}';
hold on
set(gca, 'YColor', 'k')  % Set the right y-axis to black
xlabel('Interception Angle (deg)')
title('Deflection Distance and GAIN(%) of different Interception Angle')
% poly fit
p1 = polyfit(xticks(5:15),mean_BIP_degree(5:15),1);
p2 = polyfit(xticks(16:23),mean_BIP_degree(16:23),1);
p3 = polyfit(xticks(5:23),mean_COG_degree(5:23),2);
syms  x
y1 = p1(1)*x + p1(2);
y2 = p2(1)*x + p2(2);
y3 = p3(1)*x^2 + p3(2)*x + p3(3);
y4 = ((y2 - y3) / y3)*100;
y5 = ((y1 - y3) / y3)*100;
fplot(x,y1,[-25,-2.5],'blue:')
hold on
fplot(x,y2,[2.5,20],'blue:')
hold on
fplot(x,y3,[-25,20],'red:')
hold on
h1 = fill([4.22, 6.22, 6.22, 4.22], ...
     [0, 0, 3000, 3000], [0.9, 0.9, 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5); 
m4 = 'Initial Angle';
hold on

yyaxis right
fplot(x,y4,[0,20.28],'black:')
fplot(x,y5,[-24.78,0],'black:')
ylim([0 4500])
lgd5 = legend([plot1,plot2,plot3,h1],m1,m2,m3,m4);
set(lgd5, 'Position', [0.62,0.8,0.1,0.05]);
yyaxis left
equation_str1 = sprintf('y = %.2fx + %.2f', p1(1), p1(2));
equation_str2 = sprintf('y = %.2fx + %.2f', p2(1), p2(2));
equation_str3 = sprintf('y = %.2fx^2 %.2fx + %.2f', p3(1), p3(2), p3(3));
txt1 = text(-20, 2800, equation_str1, 'FontSize', 12, 'Color', 'blue'); % Adjust the position [x, y] as needed
txt2 = text(0, 2200, equation_str2, 'FontSize', 12, 'Color', 'blue'); % Adjust the position [x, y] as needed
txt3 = text(-13, 650, equation_str3, 'FontSize', 12, 'Color', 'red'); % Adjust the position [x, y] as needed

% Export graph

exportgraphics(gca2, [pathTosave,'Gain_rate_from_BIP_to_COG_from_different_interception_angle','.png']);

% ------------------------------------------------------------------------------------------------------------






% -------------------------different shape-----------------------------------------------------------------------------
for i = 1:12
    name = [cell2mat(shape_data.Model(i))];
    filename = [name,'_Model_distribution.mat']; % Change to your desired file name
    % Full path for the file
    fullFilePath = fullfile(pathOfmatfile_shape, filename);

    % Check if the file exists
    if isfile(fullFilePath)
        % File exists, read it using readtable
        load(fullFilePath);
        disp('File read successfully.');
    else
        % File does not exist
        error([ filename, ' not found in ',pathOfmatfile_shape]);
    end

    boxplot_shape(:,2*i-1) = delta_r_min_COG./1000; % km
    boxplot_shape(:,2*i) = delta_r_min_BIP./1000; % km


end


% plot graph
% Both COG & BIP strategies
gca3 = figure(3);
ax1 = boxplot(boxplot_shape);
ylabel('Deflection Distance (km)')
title('Candle distribution of Both BIP and COG ')
xPositions = 1:24; % Modify if necessary

% Set custom XTick positions and labels, modify if necessary
set(gca, 'XTick', [1.5, 3.5, 5.5, 7.5, 9.5, 11.5, 13.5, 15.5, 17.5, 19.5, 21.5, 23.5], ...
    'XTickLabel', {'Apophis', 'Iris Hanus', 'Buda', 'Deira', 'Lutetia', 'Vesta', 'Itokawa', 'Bennu', ...
    'Golevka', 'Castalia', 'Betulia', 'Bacchus'});

% Hold the current plot to overlay the background colors
hold on;

% Highlight the first 8 pairs for Light-Curve (positions 1-8)
h1 = fill([xPositions(1)-1, xPositions(8)+0.5, xPositions(8)+0.5, xPositions(1)-1], ...
     [-450,-450,1550,1550], [0.9, 0.9, 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5); % Light blue

% Highlight the middle 8 pairs for Satellite-Orbservation (positions 9-16)
h2 = fill([xPositions(9)-0.5, xPositions(16)+0.5, xPositions(16)+0.5, xPositions(9)-0.5], ...
     [-450,-450,1550,1550], [0.9, 1, 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5); % Light green

% Highlight the last 8 pairs for Rader-based (positions 17-24)
h3 = fill([xPositions(17)-0.5, xPositions(24)+0.5, xPositions(24)+0.5, xPositions(17)-0.5], ...
     [-450,-450,1550,1550], [1, 0.9, 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5); % Light red

% Add a legend to describe each background color
lgd = legend([h1 h2 h3], {'Light-curve', 'Satellite Orbserve', 'Rader-Based'}, 'Location', 'BestOutside');

set(gcf, 'Position', [100, 500, 1400, 500]);
set(lgd, 'Position', [0.3,0.162,0,0]);

% Release the hold on the plot
hold off;


% plot only COG
gca4 = figure(4);
boxplot(boxplot_shape(:,1:2:24))
ylabel('Deflection Distance (km)')
title('Candle distribution of COG ')
xPositions = 1:12;

% Set custom XTick positions and labels, modify if necessary
set(gca, 'XTick', 1:12 , ...
    'XTickLabel', {'Apophis', 'Iris Hanus', 'Buda', 'Deira', 'Lutetia', 'Vesta', 'Itokawa', 'Bennu', ...
    'Golevka', 'Castalia', 'Betulia', 'Bacchus'});

% Hold the current plot to overlay the background colors
hold on;

% Highlight the first 8 pairs for Light-Curve (positions 1-8)
h1 = fill([0,4.5,4.5,0], ...
     [-450,-450,1550,1550], [0.9, 0.9, 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5); % Light blue

% Highlight the middle 8 pairs for Satellite-Orbservation (positions 9-16)
h2 = fill([xPositions(5)-0.5, xPositions(8)+0.5, xPositions(8)+0.5, xPositions(5)-0.5], ...
     [-450,-450,1550,1550], [0.9, 1, 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5); % Light green

% Highlight the last 8 pairs for Rader-based (positions 17-24)
h3 = fill([xPositions(9)-0.5, xPositions(12)+0.5, xPositions(12)+0.5, xPositions(9)-0.5], ...
     [-450,-450,1550,1550], [1, 0.9, 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5); % Light red

% Add a legend to describe each background color
lgd2 = legend([h1 h2 h3], {'Light-curve', 'Satellite Orbserve', 'Rader-Based'}, 'Location', 'BestOutside');

set(gcf, 'Position', [100, 50, 700, 500]);
set(lgd2, 'Position', [0.235,0.45,0,0]);
% Release the hold on the plot
hold off;



% plot only BIP
gca5 = figure(5);
boxplot(boxplot_shape(:,2:2:24))
ylabel('Deflection Distance (km)')
title('Candle distribution of BIP ')
xPositions = 1:12;

% Set custom XTick positions and labels
set(gca, 'XTick', 1:12 , ...
    'XTickLabel', {'Apophis', 'Iris Hanus', 'Buda', 'Deira', 'Lutetia', 'Vesta', 'Itokawa', 'Bennu', ...
    'Golevka', 'Castalia', 'Betulia', 'Bacchus'});

% Hold the current plot to overlay the background colors
hold on;

% Highlight the first 8 pairs for Light-Curve (positions 1-8)
h1 = fill([0,4.5,4.5,0], ...
     [-450,-450,1550,1550], [0.9, 0.9, 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5); % Light blue

% Highlight the middle 8 pairs for Satellite-Orbservation (positions 9-16)
h2 = fill([xPositions(5)-0.5, xPositions(8)+0.5, xPositions(8)+0.5, xPositions(5)-0.5], ...
     [-450,-450,1550,1550], [0.9, 1, 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5); % Light green

% Highlight the last 8 pairs for Rader-based (positions 17-24)
h3 = fill([xPositions(9)-0.5, xPositions(12)+0.5, xPositions(12)+0.5, xPositions(9)-0.5], ...
     [-450,-450,1550,1550], [1, 0.9, 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5); % Light red

% Add a legend to describe each background color
lgd = legend([h1 h2 h3], {'Light-curve', 'Satellite Orbserve', 'Rader-Based'}, 'Location', 'BestOutside');
set(gcf, 'Position', [800, 50, 700, 500]);
set(lgd, 'Position', [0.3,0.162,0,0]);
% Release the hold on the plot
hold off;

% Export graph
exportgraphics(gca3, [pathTosave,'Candles_distribution_of_different_shape_both_COG_and_BIP','.png']);
exportgraphics(gca4, [pathTosave,'Candles_distribution_of_different_shape_both_COG','.png']);
exportgraphics(gca5, [pathTosave,'Candles_distribution_of_different_shape_both_BIP','.png']);
%---------------------------------------------------------------------------------------------------------------






