clear; clc; close all;

%% === PARAMETERS ===
numSamplePoints = 200;
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

%% === Plot 4: Final Error Distance (Mean Across 5 Runs) ===
finalErrors = zeros(length(files),1);
for i = 1:length(files)
    vals = zeros(5,1);
    for j = 1:5
        data = readtable(files{i}{j});
        data.Properties.VariableNames = strtrim(data.Properties.VariableNames);
        ed = data.Error_Distance;
        ed = ed(~isnan(ed));
        vals(j) = ed(end);
    end
    finalErrors(i) = mean(vals);
end

figure(4);
bar(finalErrors, 'FaceColor', 'flat'); colormap(colors);
title('Mean Final Error Distance per Strategy');
ylabel('Error Distance'); set(gca, 'XTickLabel', labels);
xtickangle(45); grid on;

%% === Plot 5: Error Over Time (Mean Curve) ===
figure(5); hold on;
title('Error Over Time (Mean Across Reps)'); xlabel('Time'); ylabel('Error'); grid on;

for i = 1:length(files)
    curves = [];
    for j = 1:5
        data = readtable(files{i}{j});
        data.Properties.VariableNames = strtrim(data.Properties.VariableNames);
        x = data.Elapsed_Time_s;
        y = data.Error_Distance;
        curves(j,1:length(y)) = y(:)';
    end
    meanCurve = mean(curves, 1, 'omitnan');
    t = 1:length(meanCurve);
    plot(t, meanCurve, 'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName', labels{i});
end
legend('Location', 'best'); hold off;

%% === Plot 6: Total Movement Distance (Mean Across 5) ===
totalDistances = zeros(length(files),1);
for i = 1:length(files)
    dists = zeros(5,1);
    for j = 1:5
        data = readtable(files{i}{j});
        data.Properties.VariableNames = strtrim(data.Properties.VariableNames);
        x = data.X_Position;
        y = data.Y_Position;
        dx = diff(x); dy = diff(y);
        dists(j) = sum(sqrt(dx.^2 + dy.^2), 'omitnan');
    end
    totalDistances(i) = mean(dists);
end

figure(6);
bar(totalDistances, 'FaceColor', 'flat'); colormap(colors);
title('Total Distance Traveled (Mean)'); ylabel('Distance');
set(gca, 'XTickLabel', labels); xtickangle(45); grid on;

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
