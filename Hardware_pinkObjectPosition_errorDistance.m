clear; clc; close all;

%% === PARAMETERS ===
numSamplePoints = 100;    % Resample trajectories to this length
distThreshold = 0.1;      % <-- Change this to control max allowed distance in clusters

%% Load the Excel Files 
files = {
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\Hardware\THESIS_TESTS\Preliminary_tests\LOCAL.xlsx',
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\Hardware\THESIS_TESTS\Preliminary_tests\NEIGHBOUR.xlsx',
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\Hardware\THESIS_TESTS\Preliminary_tests\SELFISH.xlsx',
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\Hardware\THESIS_TESTS\Preliminary_tests\GLOBAL.xlsx', 
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\Hardware\THESIS_TESTS\Preliminary_tests\GLOBAL_ONLY.xlsx',
};

colors = lines(length(files));
labels = cell(length(files), 1);
trajectories = cell(length(files), 1);

%% === Load and Resample Trajectories ===
for i = 1:length(files)
    [~, fname, ~] = fileparts(files{i});
    cleanName = upper(strrep(fname, '_', ''));
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

    data = readtable(files{i});
    data.Properties.VariableNames = strtrim(data.Properties.VariableNames);
    
    if all(ismember({'X_Position', 'Y_Position'}, data.Properties.VariableNames))
        x = data.X_Position;
        y = data.Y_Position;
        
        % Resample to same length using linear interpolation
        origLen = length(x);
        tOrig = linspace(0,1,origLen);
        tNew = linspace(0,1,numSamplePoints);
        x_resampled = interp1(tOrig, x, tNew);
        y_resampled = interp1(tOrig, y, tNew);
        
        trajectories{i} = [x_resampled(:), y_resampled(:)];
    else
        warning('Missing X/Y position in: %s', files{i});
        trajectories{i} = nan(numSamplePoints, 2);
    end
end

%% === Plot 1: Raw Trajectories ===
figure(1);
hold on;
title('Object Trajectories (Resampled)');
xlabel('X Position');
ylabel('Y Position');
grid on;

for i = 1:length(trajectories)
    plot(trajectories{i}(:,1), trajectories{i}(:,2), 'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName', labels{i});
end
legend('Location', 'best');
hold off;

%% === Calculate Procrustes Distance Matrix ===
n = length(trajectories);
distMatrix = zeros(n);

for i = 1:n
    for j = i+1:n
        [d, ~, ~] = procrustes(trajectories{i}, trajectories{j});  % <-- FIXED HERE
        distMatrix(i,j) = d;
        distMatrix(j,i) = d;
    end
end

%% === Plot 2: Heatmap of Pairwise Distances ===  dark = more different.
figure(2);
imagesc(distMatrix);
colorbar;
title('Procrustes Distance Matrix (Trajectories)');
set(gca, 'XTick', 1:n, 'XTickLabel', labels, 'XTickLabelRotation', 45);
set(gca, 'YTick', 1:n, 'YTickLabel', labels);
axis square;

%% === Plot 3: Dendrogram + Cluster Mean Shapes Subplots === (groups trajectories by similarity)
% Left subplot: the hierarchical dendrogram showing clustering hierarchy and labels.
%Right subplot: overlaid trajectories colored by cluster, with each cluster5s mean shape plotted boldly.

Y = squareform(distMatrix);
Z = linkage(Y, 'average');

% Cluster using distance threshold instead of fixed number of clusters
clusterIDs = cluster(Z, 'cutoff', distThreshold, 'criterion', 'distance');

uniqueClusters = unique(clusterIDs);
numClusters = length(uniqueClusters);

clusterColors = lines(numClusters);

figure('Position',[100 100 900 400]);

% Subplot 3a: Dendrogram
subplot(1,2,1);
dendrogram(Z, 'Labels', labels);
title('Trajectory Clustering (Dendrogram)');
set(gca,'XTickLabelRotation',45);
hold on;
% Draw a horizontal line to show cut-off
yl = ylim;
plot(xlim, [distThreshold distThreshold], 'r--', 'LineWidth', 1.5);
hold off;

% Subplot 3b: Mean Shapes + Trajectories by Cluster
subplot(1,2,2);
hold on;
title('Cluster Mean Shapes with Trajectories');
xlabel('X Position');
ylabel('Y Position');
grid on;

for c = 1:numClusters
    idx = find(clusterIDs == uniqueClusters(c));
    ref = trajectories{idx(1)};
    aligned = zeros(size(ref,1), 2, length(idx));
    
    for k = 1:length(idx)
        [~, Zmat] = procrustes(ref, trajectories{idx(k)});
        aligned(:,:,k) = Zmat;
        plot(Zmat(:,1), Zmat(:,2), 'Color', [clusterColors(c,:) 0.3]);
    end
    
    meanShape = mean(aligned, 3);
    plot(meanShape(:,1), meanShape(:,2), 'Color', clusterColors(c,:), 'LineWidth', 2.5, 'DisplayName', ['Cluster ' num2str(uniqueClusters(c))]);
end

legend show;
hold off;


%% === Plot 4: Final Error Distance Bar Chart ===   taking the final error (last data) for strategies and comparing
finalErrors = zeros(length(files), 1);  % preallocate

for i = 1:length(files)
    data = readtable(files{i});
    data.Properties.VariableNames = strtrim(data.Properties.VariableNames);

    if ismember('Error_Distance', data.Properties.VariableNames)
        ed = data.Error_Distance;
        ed = ed(~isnan(ed));  % remove NaNs if any
        if ~isempty(ed)
            finalErrors(i) = ed(end); % last error distance
        else
            finalErrors(i) = NaN;
            warning('No valid Error_Distance data in file: %s', files{i});
        end
    else
        finalErrors(i) = NaN;
        warning('Missing Error_Distance in: %s', files{i});
    end
end

figure(4);
bar(finalErrors, 'FaceColor', 'flat');
colormap(colors); % apply same colors as other plots
title('Final Error Distance per Strategy');
ylabel('Error Distance at Last Data Point');
set(gca, 'XTickLabel', labels, 'XTick', 1:length(labels));
xtickangle(45);
grid on;


%% === Plot 5: Error Evolution Over Time ===
figure(5);
hold on;
title('Error Evolution Over Time');
xlabel('Sample Number');
ylabel('Error (Error_Python)');
grid on;

for i = 1:length(files)
    data = readtable(files{i});
    data.Properties.VariableNames = strtrim(data.Properties.VariableNames);
    
    if all(ismember({'Elapsed_Time_s', 'Error_Distance'}, data.Properties.VariableNames))
        x = data.Elapsed_Time_s;
        y = data.Error_Distance;
        plot(x, y, 'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName', labels{i});
    else
        warning('Missing Sample_Number or Error_Python in file: %s', files{i});
    end
end

legend('Location', 'best');
hold off;


%% === Plot 6: Total Movement Distance per Strategy === 
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
figure(6);
bar(totalDistances, 'FaceColor', 'flat');
colormap(colors);
title('Total Distance Traveled per Strategy');
ylabel('Total Distance Traveled (pixels)');
set(gca, 'XTickLabel', labels, 'XTick', 1:length(labels));
xtickangle(45);
grid on;



%% === Plot 7: Productivity vs Cumulative Distance over Time === 
figure(7);
numFiles = length(files);
tiledlayout(numFiles, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

for i = 1:numFiles
    nexttile;
    
    data = readtable(files{i});
    data.Properties.VariableNames = strtrim(data.Properties.VariableNames);

    if all(ismember({'X_Position', 'Y_Position', 'Error_Distance'}, data.Properties.VariableNames))
        x = data.X_Position;
        y = data.Y_Position;
        errorDist = data.Error_Distance;
        
        validIdx = ~isnan(errorDist);
        x = x(validIdx);
        y = y(validIdx);
        errorDist = errorDist(validIdx);
        
        dx = diff(x);
        dy = diff(y);
        stepDistances = sqrt(dx.^2 + dy.^2);
        cumDist = [0; cumsum(stepDistances)];
        
        deltaError = errorDist(1:end-1) - errorDist(2:end);
        deltaDist = stepDistances;
        
        productivity = deltaError ./ deltaDist;
        productivity(deltaDist == 0) = 0;
        
        plot(cumDist(2:end), productivity, 'Color', colors(i,:), 'LineWidth', 1.5);
        grid on;
        title(labels{i}, 'Interpreter', 'none');
        
        % Only bottom subplot has xlabel
        if i == numFiles
            xlabel('Cumulative Distance Traveled');
        else
            set(gca, 'XTickLabel', []);
        end
        
        % Only first subplot has ylabel
        if i == 1
            ylabel('Productivity (Error Improvement per Unit Distance)');
        else
            set(gca, 'YTickLabel', []);
        end
    else
        warning('Missing required columns in file: %s', files{i});
    end
end
