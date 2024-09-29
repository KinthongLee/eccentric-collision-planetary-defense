clearvars
clc
% ---------------------------  Add Path  -----------------------------------------------
currentDir = fileparts(which('plot_launch_window.m'));

% add path to save
pathTosave = [fullfile(currentDir, 'output_result', 'different_PHA','distribution','100k_result','pngfile'),'\'];

% add path of the .mat file stored
pathOfmatfile = [fullfile(currentDir, 'output_result', 'different_PHA','distribution','100k_result','matfile'),'\'];

% ---------------------------------------------------------------------------------------


% --------------------------Read Data from Table-----------------------------------------
data = readtable('PHA_table');
% ---------------------------------------------------------------------------------------


 % Plot all the results from No.1 to No.32
for i = 32:32
    if i ~= 9 && i ~= 14 % Exclude No.9 & No.14, due to unavailable Lambert Transfer
        name = cell2mat(data.Object(i));
        filename = [name,'_distribution.mat']; % Change to your desired file name
        % Full path for the file
        fullFilePath = fullfile(pathOfmatfile, filename);
        
        % Check if the file exists
        if isfile(fullFilePath)
            % File exists, read it using readtable
            load(fullFilePath);
            disp('File read successfully.');
        else
            % File does not exist
            error([ filename, ' not found in ',pathOfmatfile]);
        end
    
        % Deflection distance
        delta_r_min_BIP = delta_r_min_BIP./1000; % BIP, km
        delta_r_min_COG = delta_r_min_COG./1000; % COG (Center of Mass), km
        % Fix the scale of x-axis to the 0.01% to 99.99% of the total
        % samples(include both BIP & COG strategy)
        minimum = min( prctile(delta_r_min_BIP, 0.01),prctile(delta_r_min_COG, 0.01) );
        maximum = max( prctile(delta_r_min_BIP, 99.99),prctile(delta_r_min_COG, 99.99) );

        % Plot graph
        % BIP strategy
        gca = figure(1);
        ax1 = subplot(2,1,1);
        histogram(delta_r_min_BIP);
        title([num2str(i),'.',name])
        ylabel('BIP')
        xlim([minimum maximum])
        set(ax1, 'XTickLabel', []); % Remove x-axis labels
 
        % COG strategyy
        ax2 =  subplot(2,1,2);
        histogram(delta_r_min_COG);
        xlabel('Deflection Distance (km)')
        ylabel('COG')
        xlim([minimum maximum])

        % Adjust positions
        ax1_pos = get(ax1, 'Position'); % Get the position of the first subplot
        ax2_pos = get(ax2, 'Position'); % Get the position of the second subplot

        % Reduce the vertical distance between the plots
        ax2_pos(2) = ax1_pos(2) - ax2_pos(4) - 0.04; % Move up the second subplot closer to the first

        % Set the new positions
        set(ax1, 'Position', ax1_pos);
        set(ax2, 'Position', ax2_pos);

        set(gcf, 'Units', 'inches');
        set(gcf, 'Position', [1, 1, 3.2, 3])

        exportgraphics(gca, [pathTosave,num2str(i),'_distribution','.png']);
        close 
        


    end


end