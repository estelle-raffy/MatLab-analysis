clear; clc; close all;

%% === PARAMETERS ===
numSamplePoints = 200; % Resample trajectories to this length, ONLY USING error_distance from Python not Arduino! 
distThreshold = 0.1;

% === Set Folder Where Excel Files Are Stored ===
baseFolder = 'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\Hardware\THESIS_TESTS\5 repeats Preliminary\300x300_1xELONGATE_length'; 

%% === Load Excel File Paths ===
% Each strategy has 5 repetitions
files = {
    fullfile(baseFolder, {'LOCAL1.xlsx', 'LOCAL2.xlsx', 'LOCAL3.xlsx', 'LOCAL4.xlsx', 'LOCAL5.xlsx'});
    fullfile(baseFolder, {'NEIGHBOUR1.xlsx', 'NEIGHBOUR2.xlsx', 'NEIGHBOUR3.xlsx', 'NEIGHBOUR4.xlsx', 'NEIGHBOUR5.xlsx'});
    fullfile(baseFolder, {'SELFISH1.xlsx', 'SELFISH2.xlsx', 'SELFISH3.xlsx', 'SELFISH4.xlsx', 'SELFISH5.xlsx'});
    fullfile(baseFolder, {'GLOBAL1.xlsx', 'GLOBAL2.xlsx', 'GLOBAL3.xlsx', 'GLOBAL4.xlsx', 'GLOBAL5.xlsx'});
    fullfile(baseFolder, {'GLOBAL_ONLY1.xlsx', 'GLOBAL_ONLY2.xlsx', 'GLOBAL_ONLY3.xlsx', 'GLOBAL_ONLY4.xlsx', 'GLOBAL_ONLY5.xlsx'});
};

strategyNames = {'LOCAL', 'NEIGHBOUR', 'SELFISH', 'GLOBAL', 'GLOBAL_ONLY'};
colors = lines(length(files));

meanTrajectories = cell(length(files),1);
allTrajectories = cell(length(files),1);  % Cell of cells
labels = strategyNames;

%% === Load, Resample, and Average Trajectories per Strategy ===
for i = 1:length(files)
    group = files{i};
    trajGroup = cell(length(group),1);
    
    for j = 1:length(group)
        data = readtable(group{j});
        data.Properties.VariableNames = strtrim(data.Properties.VariableNames);
        
        if all(ismember({'X_Position', 'Y_Position'}, data.Properties.VariableNames))
            x = data.X_Position;
            y = data.Y_Position;
            
            origLen = length(x);
            tOrig = linspace(0,1,origLen);
            tNew = linspace(0,1,numSamplePoints);
            x_resampled = interp1(tOrig, x, tNew);
            y_resampled = interp1(tOrig, y, tNew);
            
            trajGroup{j} = [x_resampled(:), y_resampled(:)];
        else
            trajGroup{j} = nan(numSamplePoints, 2);
            warning('Missing XY data in: %s', group{j});
        end
    end
    
    % Store all and mean
    allTrajectories{i} = trajGroup;
    trajMat = cat(3, trajGroup{:});
    meanTrajectories{i} = mean(trajMat, 3, 'omitnan');
end

%% === Plot 1: Raw Trajectories (All Repetitions + Mean) ===
figure(1); hold on;
title('Pink endcap (Motor 1) Trajectories (Resampled + mean in bold)');
xlabel('X Position'); ylabel('Y Position'); grid on;

for i = 1:length(allTrajectories)
    for j = 1:length(allTrajectories{i})
        plot(allTrajectories{i}{j}(:,1), allTrajectories{i}{j}(:,2), ...
             'Color', [colors(i,:) 0.3], 'HandleVisibility', 'off');
    end
    plot(meanTrajectories{i}(:,1), meanTrajectories{i}(:,2), ...
         'Color', colors(i,:), 'LineWidth', 2.0, 'DisplayName', labels{i});
end
legend('Location', 'best'); hold off;

%% === Plot 2: Procrustes Distance Heatmap (Mean Shapes) ===
n = length(meanTrajectories);
distMatrix = zeros(n);
for i = 1:n
    for j = i+1:n
        [d, ~, ~] = procrustes(meanTrajectories{i}, meanTrajectories{j});
        distMatrix(i,j) = d;
        distMatrix(j,i) = d;
    end
end

figure(2);
imagesc(distMatrix);
colorbar;
title('Procrustes Distance Matrix (Mean Trajectories)');
set(gca, 'XTick', 1:n, 'XTickLabel', labels, 'XTickLabelRotation', 45);
set(gca, 'YTick', 1:n, 'YTickLabel', labels);
axis square;

%% === Plot 3: Dendrogram + Clustering (Mean Shapes) ===
Y = squareform(distMatrix);
Z = linkage(Y, 'average');
clusterIDs = cluster(Z, 'cutoff', distThreshold, 'criterion', 'distance');

uniqueClusters = unique(clusterIDs);
numClusters = length(uniqueClusters);
clusterColors = lines(numClusters);

figure('Position',[100 100 900 400]);

% Dendrogram
subplot(1,2,1);
dendrogram(Z, 'Labels', labels);
title('Trajectory Clustering (Dendrogram)');
set(gca,'XTickLabelRotation',45);
hold on;
plot(xlim, [distThreshold distThreshold], 'r--', 'LineWidth', 1.5);
hold off;

% Cluster Plot
subplot(1,2,2); hold on;
title('Cluster Mean Trajectories Shapes');
xlabel('X'); ylabel('Y'); grid on;

for c = 1:numClusters
    idx = find(clusterIDs == uniqueClusters(c));
    ref = meanTrajectories{idx(1)};
    aligned = zeros(size(ref,1), 2, length(idx));
    
    for k = 1:length(idx)
        [~, Zmat] = procrustes(ref, meanTrajectories{idx(k)});
        aligned(:,:,k) = Zmat;
        plot(Zmat(:,1), Zmat(:,2), 'Color', [clusterColors(c,:) 0.3], 'HandleVisibility', 'off');
    end
    
    meanShape = mean(aligned,3);
    plot(meanShape(:,1), meanShape(:,2), 'Color', clusterColors(c,:), 'LineWidth', 2.5, ...
         'DisplayName', ['Cluster ' num2str(uniqueClusters(c))]);
end
legend show; hold off;

%% === Plot 4: Final Error Distance (boxplot showing variability across 5 runs) === 

numStrategies = length(files);
numRepeats = 5;
finalErrors = nan(numStrategies, numRepeats);  % store all runs

for i = 1:numStrategies
    for j = 1:numRepeats
        data = readtable(files{i}{j});
        data.Properties.VariableNames = strtrim(data.Properties.VariableNames);
        ed = data.Error_Distance; % error_distance not global error because nice to see time at rest before tensegrity moves!
        ed = ed(~isnan(ed));
        finalErrors(i, j) = ed(end);
    end
end

% Reshape data for boxplot
errorDataReshaped = [];
groupLabels = {};

for i = 1:numStrategies
    errorDataReshaped = [errorDataReshaped; finalErrors(i, :)'];
    groupLabels = [groupLabels; repmat(labels(i), numRepeats, 1)];
end

figure(4);
boxplot(errorDataReshaped, groupLabels, 'LabelOrientation', 'inline');
title('Final Error Distance per Strategy');
ylabel('Error Distance');
xtickangle(45);
grid on;

%% === Plot 5: Error Over Time (boxplot showing variability across 5 runs, Aligned by Time) ===
figure(5); clf; hold on;
title('Error Over Time (Mean Across Reps)');
xlabel('Elapsed Time (s)');
ylabel('Error (Error_Distance)');
grid on;

% Define common time base (e.g. 200 points between 0 and max time)
numPoints = 200;

for i = 1:length(files)
    curves = zeros(5, numPoints);
    tMax = 0;

    for j = 1:5
        data = readtable(files{i}{j});
        data.Properties.VariableNames = strtrim(data.Properties.VariableNames);
        t = data.Elapsed_Time_s;
        y = data.Error_Distance;

        valid = ~isnan(t) & ~isnan(y);
        t = t(valid);
        y = y(valid);

        if isempty(t) || isempty(y)
            continue;
        end

        tMax = max(tMax, max(t));
        tNorm = linspace(min(t), max(t), numPoints);
        yInterp = interp1(t, y, tNorm, 'linear', 'extrap');
        curves(j,:) = yInterp;
    end

    % Compute mean and std
    meanCurve = mean(curves, 1, 'omitnan');
    stdCurve = std(curves, 0, 1, 'omitnan');
    tPlot = linspace(0, tMax, numPoints);

    % Shaded region (±std), suppress from legend
    hFill = fill([tPlot fliplr(tPlot)], ...
                 [meanCurve + stdCurve, fliplr(meanCurve - stdCurve)], ...
                 colors(i,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    hFill.Annotation.LegendInformation.IconDisplayStyle = 'off';

    % Mean line (shown in legend)
    plot(tPlot, meanCurve, 'Color', colors(i,:), 'LineWidth', 1.5, ...
         'DisplayName', labels{i});
end

legend('Location', 'best');
hold off;


%% === Plot 6: Total Movement Distance (boxplot showing variability across 5 runs) ===

numStrategies = length(files);
numRepeats = 5;
totalDistances = nan(numStrategies, numRepeats);

for i = 1:numStrategies
    for j = 1:numRepeats
        data = readtable(files{i}{j});
        data.Properties.VariableNames = strtrim(data.Properties.VariableNames);
        x = data.X_Position;
        y = data.Y_Position;
        dx = diff(x); dy = diff(y);
        totalDistances(i, j) = sum(sqrt(dx.^2 + dy.^2), 'omitnan');
    end
end

% Reshape for boxplot
distancesReshaped = [];
groupLabels = {};

for i = 1:numStrategies
    distancesReshaped = [distancesReshaped; totalDistances(i, :)'];
    groupLabels = [groupLabels; repmat(labels(i), numRepeats, 1)];
end

figure(6);
boxplot(distancesReshaped, groupLabels, 'LabelOrientation', 'inline');
title('Total Distance Traveled');
ylabel('Distance');
xtickangle(45); grid on;


%% === Plot 7: Productivity vs Cumulative Distance ===
figure(7);
tiledlayout(length(files), 1, 'Padding', 'compact', 'TileSpacing', 'compact');

for i = 1:length(files)
    prodMat = [];
    distMat = [];
    for j = 1:5
        data = readtable(files{i}{j});
        data.Properties.VariableNames = strtrim(data.Properties.VariableNames);
        x = data.X_Position;
        y = data.Y_Position;
        err = data.Error_Distance;
        
        valid = ~isnan(err);
        x = x(valid); y = y(valid); err = err(valid);
        
        dx = diff(x); dy = diff(y);
        stepDist = sqrt(dx.^2 + dy.^2);
        cumDist = [0; cumsum(stepDist)];
        deltaError = err(1:end-1) - err(2:end);
        productivity = deltaError ./ stepDist;
        productivity(stepDist == 0) = 0;

        % Interpolate to common length for averaging
        tNew = linspace(0, 1, 100);
        prodMat(j,:) = interp1(linspace(0,1,length(productivity)), productivity, tNew, 'linear', 'extrap');
        distMat(j,:) = interp1(linspace(0,1,length(cumDist(2:end))), cumDist(2:end), tNew, 'linear', 'extrap');
    end
    
    nexttile;
    plot(mean(distMat,1), mean(prodMat,1), 'Color', colors(i,:), 'LineWidth', 1.5);
    % Compute and display mean productivity value
    meanProdValue = mean(prodMat(:), 'omitnan');
    text(0.98, 0.95, sprintf('Mean: %.3f', meanProdValue), ...
    'Units', 'normalized', 'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', 'FontSize', 9, 'Color', [0.2 0.2 0.2]);

    grid on; title(labels{i});
    
    if i == 1
        ylabel('Productivity');
    else
        set(gca,'YTickLabel',[]);
    end
    if i == length(files)
        xlabel('Cumulative Distance');
    else
        set(gca,'XTickLabel',[]);
    end
end

%% === Plot 8: Productivity Summary Boxplot (Normalized by Distance) ===

% Preallocate matrices:
numStrategies = length(files);
numRepeats = 5;
productiveSum = nan(numStrategies, numRepeats);    % Sum of negative productivity (productive phases)
unproductiveSum = nan(numStrategies, numRepeats);  % Sum of positive productivity (unproductive phases)

for i = 1:numStrategies
    for j = 1:numRepeats
        data = readtable(files{i}{j});
        data.Properties.VariableNames = strtrim(data.Properties.VariableNames);
        
        x = data.X_Position;
        y = data.Y_Position;
        err = data.Error_Distance;
        
        valid = ~isnan(err);
        x = x(valid); y = y(valid); err = err(valid);
        
        dx = diff(x);
        dy = diff(y);
        stepDist = sqrt(dx.^2 + dy.^2);
        cumDist = [0; cumsum(stepDist)];
        
        deltaError = err(1:end-1) - err(2:end);
        productivity = deltaError ./ stepDist;
        productivity(stepDist == 0) = 0;
        
        % Normalize sums by total distance to compare runs
        % especially if they have different lengths or step counts.
        totalDist = sum(stepDist);
        
        % Sum productive and unproductive phases separately
        productiveSum(i,j) = sum(-productivity(productivity < 0)) / totalDist;    % Negative productivity summed as positive number
        unproductiveSum(i,j) = sum(productivity(productivity > 0)) / totalDist;   % Positive productivity summed as is
    end
end

% Prepare data for boxplot:
% Combine productive and unproductive data in one array
allData = [productiveSum(:); unproductiveSum(:)];

% Group labels: For each strategy, two groups: 'Productive' and 'Unproductive'
groupLabels = {};
types = {'Productive', 'Unproductive'};
for t = 1:2
    for i = 1:numStrategies
        groupLabels = [groupLabels; repmat({[labels{i} ' - ' types{t}]}, numRepeats, 1)];
    end
end

figure(8);
boxplot(allData, groupLabels, 'LabelOrientation', 'inline');
title('Normalized Productivity Summary per Strategy');
ylabel('Sum of Productivity / Total Distance');
xtickangle(45);
grid on;

% ===
% Notes:
% - Normalised productivity sums by total traveled distance per run to make comparisons (python has inconsistent # of data vs Arduino).
% - "Productive" sum aggregates all phases where error decreased per unit distance.
% - "Unproductive" sum aggregates all phases where error increased per unit distance.
% - This helps interpret how much overall productive vs unproductive effort each strategy made.
% - Boxplots show variability across the 5 repeats for each strategy.


%% === Plot 9: Final Motor Lengths (M0, M1, M2) ===

% Preallocate matrix: [strategy x repetition x motor]
numStrategies = length(files);
numRepeats = 5;
numMotors = 3;
motorData = nan(numStrategies, numRepeats, numMotors);

for i = 1:numStrategies
    for j = 1:numRepeats
        data = readtable(files{i}{j});
        data.Properties.VariableNames = strtrim(data.Properties.VariableNames);
        
        if all(ismember({'M0', 'M1', 'M2'}, data.Properties.VariableNames))
            motorData(i,j,:) = [data.M0(1), data.M1(1), data.M2(1)];
        else
            warning('Missing M0/M1/M2 in file: %s', files{i}{j});
        end
    end
end

% Reshape data for boxplot: [row = sample, col = motor+strategy]
motorDataReshaped = [];
groupLabels = {};
motorLabels = {'M0', 'M1', 'M2'};

for m = 1:numMotors
    for i = 1:numStrategies
        vals = squeeze(motorData(i,:,m))';
        motorDataReshaped = [motorDataReshaped; vals(:)];
        groupLabels = [groupLabels; repmat({[labels{i} '-' motorLabels{m}]}, numRepeats, 1)];
    end
end

% Create boxplot
figure(9);
boxplot(motorDataReshaped, groupLabels, 'LabelOrientation', 'inline');
title('Final Motor Lengths per Strategy (M0, M1, M2)');
ylabel('Final Length (mm or unit)');
xtickangle(45);
grid on;

