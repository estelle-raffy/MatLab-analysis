clear; clc; close all;

%% Load the Excel Files 
files = {
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\Hardware\THESIS_TESTS\LOCAL.xlsx',
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\Hardware\THESIS_TESTS\NEIGHBOUR.xlsx',
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\Hardware\THESIS_TESTS\SELFISH.xlsx',
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\Hardware\THESIS_TESTS\GLOBAL.xlsx', 
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\Hardware\THESIS_TESTS\GLOBAL_ONLY.xlsx',
};

% Define colours (one per file)
colors = lines(length(files));

% Initialize variables
labels = cell(length(files), 1);
meanErrors = zeros(length(files), 1);

%% === Plot 1: Trajectories ===
figure(1);
hold on;
title('Object Trajectories (X vs Y)');
xlabel('X Position');
ylabel('Y Position');
grid on;

for i = 1:length(files)
    % Get just the file name (no path or extension)
    [~, fname, ~] = fileparts(files{i});
    
    % Normalize and clean name
    cleanName = upper(strrep(fname, '_', '')); % e.g., GLOBALONLY
    switch cleanName
        case 'LOCAL'
            label = 'LOCAL';
        case 'NEIGHBOUR'
            label = 'NEIGHBOUR';
        case 'SELFISH'
            label = 'SELFISH';
        case 'GLOBAL'
            label = 'GLOBAL';
        case 'GLOBALONLY'
            label = 'GLOBAL_ONLY';
        otherwise
            label = 'UNKNOWN';
    end
    labels{i} = label;

    % Load table and strip column name whitespace
    data = readtable(files{i});
    data.Properties.VariableNames = strtrim(data.Properties.VariableNames);

    % Plot X vs Y Position
    if all(ismember({'X_Position', 'Y_Position'}, data.Properties.VariableNames))
        x = data.X_Position;
        y = data.Y_Position;
        plot(x, y, 'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName', label);
    else
        warning('Missing X/Y position in: %s', files{i});
    end

    % Get mean Error_Distance for bar chart
    if ismember('Error_Distance', data.Properties.VariableNames)
        meanErrors(i) = mean(data.Error_Distance, 'omitnan');
    else
        warning('Missing Error_Distance in: %s', files{i});
        meanErrors(i) = NaN;
    end
end

legend('Location', 'best');
hold off;

%% === Plot 2: Mean Error Distance Bar Chart ===
figure(2);
bar(meanErrors, 'FaceColor', 'flat');
colormap(colors); % apply same colors as Plot 1
title('Mean Error Distance per Strategy');
ylabel('Mean Error Distance');
set(gca, 'XTickLabel', labels, 'XTick', 1:length(labels));
xtickangle(45);
grid on;


%% === Plot 3: Total Movement Distance per Strategy === 
% total distance traveled from the (X, Y) data over time for each file. 
% often called the path length = sum of Euclidean distances between each consecutive point

totalDistances = zeros(length(files), 1);

for i = 1:length(files)
    data = readtable(files{i});
    data.Properties.VariableNames = strtrim(data.Properties.VariableNames);
    
    if all(ismember({'X_Position', 'Y_Position'}, data.Properties.VariableNames))
        x = data.X_Position;
        y = data.Y_Position;
        
        % Compute pairwise Euclidean distances and sum
        dx = diff(x);
        dy = diff(y);
        distances = sqrt(dx.^2 + dy.^2);
        totalDistances(i) = sum(distances, 'omitnan');
    else
        warning('Missing X/Y position data in: %s', files{i});
        totalDistances(i) = NaN;
    end
end

% Plot the total distances
figure(3);
bar(totalDistances, 'FaceColor', 'flat');
colormap(colors);
title('Total Distance Traveled per Strategy');
ylabel('Total Distance Traveled (pixels)');
set(gca, 'XTickLabel', labels, 'XTick', 1:length(labels));
xtickangle(45);
grid on;

