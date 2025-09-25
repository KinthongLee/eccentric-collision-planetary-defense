% Replace with your main entry point file
mainFile = 'calculate_launch_window_and_best_transfer.m';

% Get a list of all required files (functions, classes, etc.)
[files, products] = matlab.codetools.requiredFilesAndProducts(mainFile);

for k = 1:numel(files)
    disp(files{k})
end