% ---------------------------  Add Path  -----------------------------------------------
currentDir = fileparts(which('plot_launch_window.m'));

% add path to save
pathTosave = [fullfile(currentDir, 'output_result', 'different_PHA','launch_window','pngfile'),'\'];

% add path of the .mat file stored
pathOfmatfile = [fullfile(currentDir, 'output_result', 'different_PHA','launch_window','matfile'),'\'];

% ---------------------------------------------------------------------------------------


% --------------------------Read Data from Table-----------------------------------------
data = readtable('PHA_table');
% ---------------------------------------------------------------------------------------



% Code starts here
% Plot the graph
for i = 1 : 32
    name = cell2mat(data.Object(i));
    filename = [name,'.mat']; % Change to your desired file name
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

    % Find the best launch window and mark it with a red circle
    [largestElement, index] = max(delta_r_min(:));
    [row, col] = ind2sub(size(delta_r_min), index);
    fig = figure(1);
    fig.Units = 'inches';
    % To better fit the paper, the images have been restricted to a small size, but they can be adjusted as needed.
    fig.Position = [1, 1, 3, 3]; % [left, bottom, width, height]
    scatter3(possible_launch(row, col),possible_reach(row, col),delta_r_min(row, col)./1000,'red');
    m1 = 'Best Launch Window';
    title([num2str(i),'.',name]) 

    hold on
    
    % Plot all the possible launch window
    surf(possible_launch,possible_reach,delta_r_min./1000,'EdgeColor','none')
    xlabel('launch date (UTC)')
    ylabel('reach date (UTC)')
    zlabel('Deflection Distance (km)')
    view(0,0)
    % Export the png file to specific file
    exportgraphics(fig, [pathTosave,num2str(i),'_launch_window.png']);
    close

end