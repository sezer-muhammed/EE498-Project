% Set the folder path where the figures will be saved
folderPath = 'ode45_2';

% Get handles of all open figures
figHandles = findall(0, 'Type', 'figure');

% Loop through each figure and save as PNG
for i = 1:numel(figHandles)
    % Get the current figure handle
    figHandle = figHandles(i);
    
    % Set the filename for saving (using the figure number)
    filename = sprintf('Figure_%d.png', i);
    
    % Set the full file path
    filePath = fullfile(folderPath, filename);
    
    % Save the figure as PNG
    saveas(figHandle, filePath, 'png');
end
