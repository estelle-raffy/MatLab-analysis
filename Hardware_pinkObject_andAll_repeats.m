clear; clc; close all;

%% === PARAMETERS ===
numSamplePoints = 200; % Resample trajectories to this length, ONLY USING error_distance from Python not Arduino! 
distThreshold = 0.1;

% === Set Folder Where Excel Files Are Stored ===
baseFolder = 'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\Hardware\THESIS_TESTS\M3 Success Region\REPEATS_top'; 

%% === Load Excel File Paths ===
% Each strategy has 5 repetitions
files = {
    fullfile(baseFolder, {'16B_1_LOCAL.xlsx', '16B_2_LOCAL.xlsx', '16B_3_LOCAL.xlsx', '16B_4_LOCAL.xlsx', '16B_5_LOCAL.xlsx'});
    fullfile(baseFolder, {'16B_1_NEIGHBOUR.xlsx', '16B_2_NEIGHBOUR.xlsx', '16B_3_NEIGHBOUR.xlsx', '16B_4_NEIGHBOUR.xlsx', '16B_5_NEIGHBOUR.xlsx'});
    fullfile(baseFolder, {'16B_1_SELFISH.xlsx', '16B_2_SELFISH.xlsx', '16B_3_SELFISH.xlsx', '16B_4_SELFISH.xlsx', '16B_5_SELFISH.xlsx'});
    fullfile(baseFolder, {'16B_1_GLOBAL.xlsx', '16B_2_GLOBAL.xlsx', '16B_3_GLOBAL.xlsx', '16B_4_GLOBAL.xlsx', '16B_5_GLOBAL.xlsx'});
    fullfile(baseFolder, {'16_1_GLOBAL_ONLY.xlsx', '16_2_GLOBAL_ONLY.xlsx', '16_3_GLOBAL_ONLY.xlsx', '16_4_GLOBAL_ONLY.xlsx', '16_5_GLOBAL_ONLY.xlsx'});
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


%{

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

%}

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
ylabel('Sum of Productivity / Total Distance'); % sum of normalized productivity values (unitless ratio).
% How much error improvement (productive) or deterioration (unproductive) you got per unit of total distance traveled in that run.
% productive sum = 0.3 → You improved error by a total amount equivalent to 0.3 times your total distance traveled.
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUPPORT PLOTS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;

for i = 1:length(files)
    subplot(3, 2, i); hold on;

    cleanName = strrep(strategyNames{i}, '_', ' ');

    trialM0 = {};
    trialM1 = {};
    trialM2 = {};
    maxT = 0;

    for j = 1:length(files{i})
        data = readtable(files{i}{j});
        data.Properties.VariableNames = matlab.lang.makeValidName(strtrim(data.Properties.VariableNames));

        data = sortrows(data, {'Sample_Number', 'Motor'});
        sampleNums = unique(data.Sample_Number);
        nSteps = length(sampleNums);
        maxT = max(maxT, nSteps);

        M0 = nan(nSteps, 1);
        M1 = nan(nSteps, 1);
        M2 = nan(nSteps, 1);

        for k = 1:nSteps
            block = data(data.Sample_Number == sampleNums(k), :);
            for m = 1:height(block)
                motorLabel = strtrim(block.Motor{m});
                val = block.Final_Actuation(m);
                if contains(motorLabel, '0')
                    M0(k) = val;
                elseif contains(motorLabel, '1')
                    M1(k) = val;
                elseif contains(motorLabel, '2')
                    M2(k) = val;
                end
            end
        end

        trialM0{end+1} = M0;
        trialM1{end+1} = M1;
        trialM2{end+1} = M2;
    end

    % Pad with NaNs
    all_M0 = nan(length(trialM0), maxT);
    all_M1 = nan(length(trialM1), maxT);
    all_M2 = nan(length(trialM2), maxT);

    for j = 1:length(trialM0)
        pad = @(v) [v; nan(maxT - length(v), 1)];
        all_M0(j, :) = pad(trialM0{j});
        all_M1(j, :) = pad(trialM1{j});
        all_M2(j, :) = pad(trialM2{j});
    end

    t = 1:maxT;

    % Calculate mean and std deviation
    M0_mean = mean(all_M0, 1, 'omitnan');
    M1_mean = mean(all_M1, 1, 'omitnan');
    M2_mean = mean(all_M2, 1, 'omitnan');

    M0_std = std(all_M0, 0, 1, 'omitnan');
    M1_std = std(all_M1, 0, 1, 'omitnan');
    M2_std = std(all_M2, 0, 1, 'omitnan');

    % Plot shaded areas (mean ± std)
    fill([t fliplr(t)], [M0_mean - M0_std fliplr(M0_mean + M0_std)], ...
        [1 0 0], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
    fill([t fliplr(t)], [M1_mean - M1_std fliplr(M1_mean + M1_std)], ...
        [0 1 0], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
    fill([t fliplr(t)], [M2_mean - M2_std fliplr(M2_mean + M2_std)], ...
        [0 0 1], 'FaceAlpha', 0.4, 'EdgeColor', 'none');

    % Plot mean lines on top
    plot(t, M0_mean, 'r', 'LineWidth', 2, 'DisplayName', 'M0 Mean');
    plot(t, M1_mean, 'g', 'LineWidth', 2, 'DisplayName', 'M1 Mean');
    plot(t, M2_mean, 'b', 'LineWidth', 2, 'DisplayName', 'M2 Mean');

    title(['Actuation Over Time - ', cleanName], 'FontSize', 16);
    xlabel('Time Steps', 'FontSize', 14);
    ylabel('Final Actuation', 'FontSize', 14);
    legend({'M0 Mean', 'M1 Mean', 'M2 Mean'}, 'Location', 'northeastoutside', 'FontSize', 12);
    set(gca, 'FontSize', 14);
    hold off;
end


%% === Plot 12b: Distribution of Equilibrium Durations (all 3 motors = 0) ===
figure;

disp('Files per strategy being analyzed in Plot 12b:');
for i = 1:length(files)
    fprintf('Strategy: %s\n', strategyNames{i});
    disp(files{i});
end

for i = 1:length(files)
    allDurations = [];
    
    for j = 1:length(files{i})
        fileName = files{i}{j};  % store for printing
        data = readtable(fileName);
        data.Properties.VariableNames = matlab.lang.makeValidName(strtrim(data.Properties.VariableNames));

        % Confirm needed columns
        if ~all(ismember({'Sample_Number', 'Motor', 'Final_Actuation'}, data.Properties.VariableNames))
            warning('Missing required columns in: %s', fileName);
            continue;
        end

        % Group by Sample_Number
        [uniqueSamples, ~, sampleIdx] = unique(data.Sample_Number);
        isEquilibrium = false(length(uniqueSamples), 1);

        for k = 1:length(uniqueSamples)
            rows = (sampleIdx == k);
            group = data(rows, :);

            % Equilibrium if all three actuations are zero
            if height(group) == 3 && all(group.Final_Actuation == 0)
                isEquilibrium(k) = true;
            end
        end

        % Find durations
        padded = [0; isEquilibrium; 0];
        change = diff(padded);
        starts = find(change == 1);
        stops = find(change == -1);
        durations = stops - starts;

        if strcmp(strategyNames{i}, 'GLOBAL_ONLY') && any(durations > 0)
            fprintf('🔍 File: %s\n', fileName);
            fprintf('  Detected durations: %s\n', mat2str(durations'));
        end

        allDurations = [allDurations; durations];
    end

    subplot(3,2,i);
    histogram(allDurations, 'BinWidth', 1);
    title(strrep(['Equilibrium Durations - ', strategyNames{i}], '_', ' '), 'FontSize', 16);
    xlabel('Duration (timesteps)', 'FontSize', 14);
    ylabel('Frequency', 'FontSize', 14);
    set(gca, 'FontSize', 14);
end


%%%
%% LOCAL & NEIGHBOUR SUPPORT - productivity underover-exploring phases 
%%%

%%% === Productivity Support Plot: Actuation vs Error Over Time ===
%%% === Visualizing Productive vs Unproductive Phases ===

targetStrategy = 'LOCAL';  % Change to 'LOCAL' etc. as needed

figure;
plotCount = 1;

for i = 1:length(files)
    if strcmp(strategyNames{i}, targetStrategy)
        for j = 1:length(files{i})
            data = readtable(files{i}{j});
            data.Properties.VariableNames = matlab.lang.makeValidName(strtrim(data.Properties.VariableNames));
            
            % Ensure the data is ordered by Sample_Number and Motor
            data = sortrows(data, {'Sample_Number', 'Motor'});

            % Keep only one row per Sample_Number (e.g., take M0 always)
            % Assumes motors are in consistent order (M0, M1, M2)
            % If needed, filter by Motor == 'M0'
            uniqueSamples = unique(data.Sample_Number);
            idx = find(strcmp(strtrim(data.Motor), 'M0'));  % or use 'Motor' == '0'
            actuation = data.Final_Actuation(idx);
            error = data.Error_python(idx);
            samples = data.Sample_Number(idx);

            % Calculate productivity: Δerror / Δsample
            dErr = diff(error);
            dT = diff(samples);  % should be 1 usually
            productivity = -dErr ./ dT;  % Positive if error is decreasing
            productivity(dT == 0) = 0;

            % Label as productive or unproductive
            isProductive = productivity > 0;
            isUnproductive = productivity < 0;

            % Prepare plot
            subplot(length(files{i}), 1, plotCount); hold on;

            % Plot error
            yyaxis right;
            plot(samples(2:end), error(2:end), 'b-', 'LineWidth', 1.2);
            ylabel('Error');

            % Plot actuation
            yyaxis left;
            plot(samples(2:end), actuation(2:end), 'k-', 'LineWidth', 1.2);
            ylabel('Actuation');

            % Highlight productive/unproductive phases
            % Slice relevant segments
            samples_sub = samples(2:end);
            actuation_sub = actuation(2:end);
            
            % Now use logical indexing on these
            scatter(samples_sub(isProductive), actuation_sub(isProductive), 20, 'g', 'filled', 'DisplayName', 'Productive');
            scatter(samples_sub(isUnproductive), actuation_sub(isUnproductive), 20, 'r', 'filled', 'DisplayName', 'Unproductive');

            title(sprintf('%s Trial %d: Productivity over Time', targetStrategy, j), 'Interpreter', 'none');
            xlabel('Sample Number');
            legend({'Actuation', 'Error', 'Productive', 'Unproductive'});
            set(gca, 'FontSize', 12);
            grid on;

            plotCount = plotCount + 1;
        end
    end
end



%%%
%% SELFISH SUPPORT - Erratic behaviour: local frustration and Wb regulation 
%%%

figure;
for i = 1:length(files)
    if strcmp(strategyNames{i}, 'SELFISH')
        for j = 1:length(files{i})
            data = readtable(files{i}{j});
            
            % Filter for M0 and M2 (where local frustration and weight are meaningful)
            validRows = ismember(strtrim(data.Motor), {'M0', 'M2'});
            sample = data.Sample_Number(validRows);
            frustration = data.Local_Frustration(validRows);
            weight = data.Weight(validRows);

            subplot(length(files{i}),2,2*j-1);
            plot(sample, frustration, 'r-', 'LineWidth', 1.2);
            title(sprintf('SELFISH Trial %d - Local Frustration (M0+M2)', j));
            xlabel('Sample Number');
            ylabel('Local Frustration');
            set(gca,'FontSize',12);

            subplot(length(files{i}),2,2*j);
            plot(sample, weight, 'b-', 'LineWidth', 1.2);
            title(sprintf('SELFISH Trial %d - Weight Regulation (M0+M2)', j));
            xlabel('Sample Number');
            ylabel('Weight');
            set(gca,'FontSize',12);
        end
    end
end

%%%
%% GLOBAL SUPPORT - better regulation? : local frustration & Wb regulation VS global frustration & Global Wb regulation
%%%

figure; hold on;
for i = 1:length(files)
    if strcmp(strategyNames{i}, 'GLOBAL') || strcmp(strategyNames{i}, 'SELFISH')
        allWeight = [];
        allFrustration = [];
        allSamples = [];

        for j = 1:length(files{i})
            data = readtable(files{i}{j});

            if strcmp(strategyNames{i}, 'GLOBAL')
                % GLOBAL: Only M1 rows matter
                validRows = ismember(strtrim(data.Motor), {'M1'});
                weight = data.Global_Neighbour_Weight(validRows);
                frustration = data.Global_Frustration(validRows);
            else
                % SELFISH: Only M0 and M2 rows matter
                validRows = ismember(strtrim(data.Motor), {'M0', 'M2'});
                weight = data.Weight(validRows);
                frustration = data.Local_Frustration(validRows);
            end

            samples = data.Sample_Number(validRows);

            allWeight = [allWeight; weight(:)];
            allFrustration = [allFrustration; frustration(:)];
            allSamples = [allSamples; samples(:)];
        end

        % Aggregate by unique sample number
        [tUnique, ~, idx] = unique(allSamples);
        meanWeight = accumarray(idx, allWeight, [], @mean, NaN);
        meanFrustration = accumarray(idx, allFrustration, [], @mean, NaN);

        if strcmp(strategyNames{i}, 'GLOBAL')
            plot(tUnique, meanWeight, 'b-', 'LineWidth', 2, 'DisplayName', 'GLOBAL Weight');
            plot(tUnique, meanFrustration, 'c--', 'LineWidth', 2, 'DisplayName', 'GLOBAL Frustration');
        else
            plot(tUnique, meanWeight, 'r-', 'LineWidth', 2, 'DisplayName', 'SELFISH Weight');
            plot(tUnique, meanFrustration, 'm--', 'LineWidth', 2, 'DisplayName', 'SELFISH Frustration');
        end
    end
end

title('Weight and Frustration Comparison: GLOBAL vs SELFISH');
xlabel('Sample Number');
ylabel('Mean Value');
legend('Location', 'best');
set(gca,'FontSize',14);
