clear; clc;

% Folder with your Excel files
dataFolder = 'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons';
files = dir(fullfile(dataFolder, '*.xlsx'));


%% DEBUG SELFISH 
disp('Checking filenames for SELFISH...');
for i = 1:length(files)
    if contains(upper(files(i).name), 'SELFISH')
        disp(['Found SELFISH file: ', files(i).name]);
    end
end

disp('Checking filenames for HOMEO...');
for i = 1:length(files)
    if contains(upper(files(i).name), 'HOMEO')
        disp(['Found HOMEO file: ', files(i).name]);
    end
end

% Define strategy groups
strategyGroups = {'LOCAL', 'NEIGHBOUR','SELFISH', 'GLOBAL_ONLY', 'GLOBAL', 'HOMEO'};
groupNames = {'Local', 'Neighbour', 'Selfish', 'M3 only', 'Global Mod', 'Homeo'};
numGroups = numel(strategyGroups);
strategyColors = lines(numGroups);

% Define spring connections
spring_pairs = [1 3; 3 2; 1 2; 2 5; 3 5; 5 4; 2 4];

theta = linspace(0, 2*pi, 20); % for balls
radius_body = 0.05;

% Expanded and compressed reference shapes
expanded_shape = [0 0; 1 0; 2 0; 3 0; 4 0];
compressed_shape = [0 0; 0.5 1; 1 0; 1.5 1; 2 0];

% ===============================================
%% ======= PLOT 1/5: Final Shapes Grouped by Strategy (Updated, Label Underneath)
% ===============================================
figure('Name','1/4: Final Shapes Grouped by Strategy');
hold on;
axis equal;
xlabel('X');
ylabel('Y');
title('1/4: Diversity of Solutions by Strategy');
grid on;
spacing_x = 5; % horizontal spacing between groups
spacing_y = 3; % vertical spacing between solutions in a group
max_per_column = 8; % maximum solutions per column
% Store a count for each group to place them nicely
groupSolutionCount = zeros(1, numGroups);
theta = linspace(0, 2*pi, 20); % for body circles
radius_body = 0.05;
% Set up consistent experiment colors (same as Plot 3/4)
expColors = lines(100);
uniqueExpNames = {}; % to track and index experiment names
% Loop through each file (experiment)
for i = 1:length(files)
    T = readtable(fullfile(files(i).folder, files(i).name), 'Sheet', 'Tabelle1');
    
    fname = upper(files(i).name);

    % Match and extract experiment name (e.g., EXP1)
    expMatch = regexp(fname, '(EXP[_\s]?0*\d+)', 'tokens', 'ignorecase');
    if ~isempty(expMatch)
        expName = upper(strrep(expMatch{1}{1}, '_', '')); % remove underscore
    else
        expName = 'Unknown';
    end
    % Get experiment color
    if ~ismember(expName, uniqueExpNames)
        uniqueExpNames{end+1} = expName;
    end
    expIdx = find(strcmp(uniqueExpNames, expName));
    color = expColors(expIdx, :);
    % If there is no success, plot a cross instead of the shapes
    if T.success_log(end) == 0
        groupIdx = [];
        for k = 1:numGroups
            if contains(fname, strategyGroups{k})
                groupIdx = k;
                break;
            end
        end
        if isempty(groupIdx)
            warning('Unknown group in file: %s', files(i).name);
            continue;
        end
        
        % Compute offset
        idx = groupSolutionCount(groupIdx);
        dx = (groupIdx - 1) * spacing_x;
        dy = -mod(idx, max_per_column) * spacing_y;
        dx = dx + floor(idx / max_per_column) * (spacing_x/2);
        groupSolutionCount(groupIdx) = groupSolutionCount(groupIdx) + 1;
        % Plot a cross
        label_x = dx + 1;
        label_y = dy;
        plot(label_x, label_y, 'x', 'MarkerSize', 8, 'LineWidth', 2, 'Color', [0.7, 0.7, 0.7]);
        % Add experiment label
        text(label_x, label_y - 0.3, expName, 'FontSize', 8, 'Color', [0.7, 0.7, 0.7], ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontWeight', 'bold');
        continue;
    end
    % Final body positions (for successful experiments)
    bx = [T.body1_pos_x(end), T.body2_pos_x(end), T.body3_pos_x(end), ...
          T.body4_pos_x(end), T.body5_pos_x(end)];
    by = [T.body1_pos_y(end), T.body2_pos_y(end), T.body3_pos_y(end), ...
          T.body4_pos_y(end), T.body5_pos_y(end)];
    % Normalization
    bx = bx - bx(1);
    by = by - by(1);
    bx = bx - min(bx);
    bx = bx / max(bx) * 2;
    by = by - min(by);
    by = by / max(by) * 1;
    % Determine group
    groupIdx = [];
    for k = 1:numGroups
        if contains(fname, strategyGroups{k})
            groupIdx = k;
            break;
        end
    end
    if isempty(groupIdx)
        warning('Unknown group in file: %s', files(i).name);
        continue;
    end
    % Compute offset
    idx = groupSolutionCount(groupIdx);
    dx = (groupIdx - 1) * spacing_x;
    dy = -mod(idx, max_per_column) * spacing_y;
    dx = dx + floor(idx / max_per_column) * (spacing_x/2);
    groupSolutionCount(groupIdx) = groupSolutionCount(groupIdx) + 1;
    % Apply offset
    bx = bx + dx;
    by = by + dy;
    % Plot springs
    for j = 1:size(spring_pairs, 1)
        i1 = spring_pairs(j,1);
        i2 = spring_pairs(j,2);
        plot([bx(i1), bx(i2)], [by(i1), by(i2)], '-', 'Color', color, 'LineWidth', 1.5);
    end
    % Plot bodies (balls)
    for j = 1:5
        fill(bx(j) + radius_body*cos(theta), by(j) + radius_body*sin(theta), ...
            color, 'FaceAlpha', 0.6, 'EdgeColor', 'k', 'LineWidth', 0.5);
    end
    % === NEW: Add EXP# Label Underneath Shape ===
    label_x = mean(bx);
    label_y = min(by) - 0.3;
    text(label_x, label_y, expName, 'FontSize', 8, 'Color', color, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontWeight', 'bold');
end
% Add group names as labels
for k = 1:numGroups
    xpos = (k-1) * spacing_x + 1;
    ypos = 2;
    text(xpos, ypos, groupNames{k}, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 12);
end
% === Ensure everything fits in the figure ===
xlim([-1, numGroups * spacing_x]);
ylim([-max_per_column * spacing_y - 1, 3]);


% ===============================================
%% ======= PLOT 2/5: Time Efficiency of Successful Solutions
% ===============================================
figure('Name','2/4: Time Efficiency of Successful Solutions');
hold on;
title('2/4: Time to Achieve Stable Success (Grouped by Strategy)');
ylabel('Time of Final Success Onset (s)');
xlabel('Strategy Group');
grid on;

% Initialize storage per group
groupSuccessTimes = cell(1, numGroups);
groupExpNames = cell(1, numGroups);
groupExpColors = cell(1, numGroups); % Store colors per experiment

% Prepare color mapping (same as Plot 1/4)
expColors = lines(100);
uniqueExpNames = {}; % make sure this matches Plot 1/4

% Loop through all files again
for i = 1:length(files)
    T = readtable(fullfile(files(i).folder, files(i).name), 'Sheet', 'Tabelle1');
    fname = upper(files(i).name);

    % Match experiment name
    expMatch = regexp(fname, '(EXP[_\s]?0*\d+)', 'tokens', 'ignorecase');
    if ~isempty(expMatch)
        expName = upper(strrep(expMatch{1}{1}, '_', ''));
    else
        expName = sprintf('File%d', i);
    end

    % Get experiment color (same method as Plot 1/4)
    if ~ismember(expName, uniqueExpNames)
        uniqueExpNames{end+1} = expName;
    end
    expIdx = find(strcmp(uniqueExpNames, expName));
    color = expColors(expIdx, :);

    % Match strategy group
    groupIdx = [];
    for k = 1:numGroups
        if contains(fname, strategyGroups{k})
            groupIdx = k;
            break;
        end
    end
    if isempty(groupIdx)
        warning('Unknown group in file: %s', files(i).name);
        continue;
    end

    % Only analyze successful runs
    if T.success_log(end) == 1
        s_log = T.success_log;
        t = T.current_time;

        % Find start of final continuous success period
        idx_end = length(s_log);
        idx_start = idx_end;
        while idx_start > 1 && s_log(idx_start - 1) == 1
            idx_start = idx_start - 1;
        end

        onset_time = t(idx_start);

        % Store info for this experiment
        groupSuccessTimes{groupIdx}(end+1) = onset_time;
        groupExpNames{groupIdx}{end+1} = expName;
        groupExpColors{groupIdx}(end+1, :) = color;
    end
end

% Plot grouped bar chart
groupSpacing = 5;
barWidth = 0.5;
x_ticks = [];
x_labels = {};

legendEntries = {};      % For legend text
legendColors = [];       % For corresponding colors
plottedExps = {};        % To avoid duplicates in legend

for k = 1:numGroups
    times = groupSuccessTimes{k};
    labels = groupExpNames{k};
    colors = groupExpColors{k};
    n = numel(times);

    for j = 1:n
        x = (k - 1) * groupSpacing + j;

        % Plot bar
        b = bar(x, times(j), barWidth, 'FaceColor', colors(j,:), 'EdgeColor', 'k');

        % Add to legend if not already included
        expName = labels{j};
        if ~ismember(expName, plottedExps)
            legendEntries{end+1} = expName;
            legendColors(end+1, :) = colors(j,:);
            plottedExps{end+1} = expName;
        end
    end

    % For group label under the bars
    mid_x = (k - 1) * groupSpacing + n / 2 + 0.5;
    x_ticks(end+1) = mid_x;
    x_labels{end+1} = groupNames{k};
end

xticks(x_ticks);
xticklabels(x_labels);
ylim([0, 65]); % Assuming max t = ~60s
xlim([0, x_ticks(end) + groupSpacing]);

% Add legend using invisible bars
hold on;
for i = 1:numel(legendEntries)
    hLegend(i) = bar(nan, nan, 'FaceColor', legendColors(i,:), 'EdgeColor', 'k');
end
legend(hLegend, legendEntries, 'Location', 'eastoutside');


% ===============================================
%% ======= PLOT 3/5: Area Dynamics Over Time (With M3 Dominance)
% ===============================================

figure('Name','3/5: Module Dynamics Over Time');
sgtitle('3/5: Evolution of Module Areas for Each Strategy');

for g = 1:numGroups
    subplot(3,3,g);
    hold on; grid on;
    title(groupNames{g});
    xlabel('Time');
    ylabel('Area Size');

    cmap = lines(3); % blue, orange, yellow for M1, M2, M3

    hM1 = []; hM2 = []; hM3 = []; % handles for legend

    for i = 1:length(files)
        T = readtable(fullfile(files(i).folder, files(i).name), 'Sheet', 'Tabelle1');

        fname = upper(files(i).name);
        if T.success_log(end) == 0, continue; end
        if ~contains(fname, strategyGroups{g}), continue; end

        A1 = T.area_M1;
        A2 = T.area_M2;

        % Get actuation magnitudes
        a1 = abs(T.M1_actuation_final(end));
        a2 = abs(T.M2_actuation_final(end));
        a3 = abs(T.M3_actuation_final(end)); % NEW for M3

        % Determine dominant module
        [~, dominantIdx] = max([a1, a2, a3]);
        color = cmap(dominantIdx, :);

        t = 1:length(A1); % assuming A1 and A2 same length

        % Plot both area curves, color by dominant
        h = plot(t, A1, '-', 'Color', color, 'LineWidth', 1);
        plot(t, A2, '-', 'Color', color, 'LineWidth', 1);

        % Save handle for legend
        switch dominantIdx
            case 1
                if isempty(hM1), hM1 = h; end
            case 2
                if isempty(hM2), hM2 = h; end
            case 3
                if isempty(hM3), hM3 = h; end
        end
    end

    % Add legend based on which handles were used
    legendEntries = {};
    legendHandles = [];

    if ~isempty(hM1)
        legendEntries{end+1} = 'M1 Dominant';
        legendHandles(end+1) = hM1;
    end
    if ~isempty(hM2)
        legendEntries{end+1} = 'M2 Dominant';
        legendHandles(end+1) = hM2;
    end
    if ~isempty(hM3)
        legendEntries{end+1} = 'M3 Dominant';
        legendHandles(end+1) = hM3;
    end

    if ~isempty(legendHandles)
        legend(legendHandles, legendEntries, 'Location', 'best');
    end
end



% ===============================================
%% ======= PLOT 4/5: Global Error vs Local Frustration (Per Experiment)
% ===============================================

% Set up color map for up to 100 experiments
expColors = lines(100);

% Track all unique experiment names
uniqueExpNames = {};

% Prepare per-experiment error and frustration values
groupedError = cell(1, numGroups);
groupedFrustration = cell(1, numGroups);
groupedExpNames = cell(1, numGroups);  % Per group

fprintf('\n=== DEBUG OUTPUT FOR PLOT 3/4 ===\n');

for i = 1:length(files)
    fname = upper(files(i).name);
    fprintf('Processing file: %s\n', fname);

    T = readtable(fullfile(files(i).folder, files(i).name), 'Sheet', 'Tabelle1');
    if T.success_log(end) == 0
        fprintf('  -> success_log(end) = 0\n');
        fprintf('  -> Skipping (no final success)\n');
        continue;
    end

    % Determine group
    groupIdx = [];
    for g = 1:numGroups
        if contains(fname, upper(strategyGroups{g}))
            groupIdx = g;
            fprintf('  -> Matched group: %s\n', groupNames{g});
            break;
        end
    end
    if isempty(groupIdx), continue; end

    % Error and frustration
    globalErr = mean(T.global_error);
    frustration = mean((T.M1_local_frustration + T.M2_local_frustration) / 2);

    groupedError{groupIdx}(end+1) = globalErr;
    groupedFrustration{groupIdx}(end+1) = frustration;

    % Match and store experiment name
    expMatch = regexp(fname, '(EXP[_\s]?0*\d+)', 'tokens', 'ignorecase');
    if ~isempty(expMatch)
        expName = upper(strrep(expMatch{1}{1}, '_', '')); % remove underscores
    else
        expName = 'Unknown';
    end
    groupedExpNames{groupIdx}{end+1} = expName;

    % Track unique experiment names globally
    if ~ismember(expName, uniqueExpNames)
        uniqueExpNames{end+1} = expName;
    end

    fprintf('  -> Matched EXP name: %s\n', expName);
end

% Summary
fprintf('\nSummary per group:\n');
for g = 1:numGroups
    fprintf('Strategy: %s | Num Experiments: %d\n', groupNames{g}, length(groupedError{g}));
end

% ===== Plotting
figure;
hold on;

numExps = length(uniqueExpNames);
data = nan(numGroups, numExps);          % rows = groups, cols = unique exps
frustrationPoints = nan(numGroups, numExps);

% Fill matrices with data matching unique EXP names
for g = 1:numGroups
    for j = 1:length(groupedExpNames{g})
        expName = groupedExpNames{g}{j};
        expIdx = find(strcmp(uniqueExpNames, expName));
        data(g, expIdx) = groupedError{g}(j);
        frustrationPoints(g, expIdx) = groupedFrustration{g}(j);
    end
end

x = 1:numGroups;

% Bar plot for global error
yyaxis left;
bh = bar(x, data, 0.5, 'stacked', 'BarWidth', 0.5);
for i = 1:numExps
    bh(i).FaceColor = expColors(i, :);
end
ylabel('Average Global Error');
set(gca, 'YColor', 'k');

% Overlay frustration lines
yyaxis right;
set(gca, 'YColor', [0.4 0.4 0.4]);

% Plot the frustration lines AFTER the bar plot
for i = 1:numExps
    % Plot each frustration dots after bars 
    plot(x, frustrationPoints(:, i), 'o', 'Color', expColors(i,:), ...
         'LineWidth', 1.5, 'MarkerFaceColor', expColors(i,:), ...
         'MarkerEdgeColor', 'k', 'MarkerSize', 6, 'LineJoin', 'round');
end

ylabel('Average Local Frustration');
set(gca, 'XTick', 1:numGroups, 'XTickLabel', groupNames);
xlabel('Control Strategy');
title('3/4: Error vs Frustration Per Experiment');
grid on;

% Legend
legend(uniqueExpNames, 'Location', 'northeastoutside');


% ===============================================
%% ======= PLOT 5/5: Integrated Information (Φ)
% ===============================================
% Set up color map and initialize
expColors = lines(100);  % Color map for up to 100 experiments
uniqueExpNames = {};  % List to store unique experiment names
phiMatrix = nan(numGroups, 100);  % Rows = strategy groups, columns = experiments
fprintf('\n=== DEBUG OUTPUT FOR PLOT 4/4 ===\n');
for i = 1:length(files)
    fname = upper(files(i).name);  % Get the experiment file name
    fprintf('Processing file: %s\n', fname);
    T = readtable(fullfile(files(i).folder, files(i).name), 'Sheet', 'Tabelle1');
    if T.success_log(end) == 0
        fprintf('  -> success_log(end) = 0\n');
        fprintf('  -> Skipping (no final success)\n');
        continue;
    end
    % Get experiment name (e.g., EXP1, EXP2)
    expMatch = regexp(fname, '(EXP[_\s]?0*\d+)', 'tokens', 'ignorecase');
    if ~isempty(expMatch)
        expName = upper(strrep(expMatch{1}{1}, '_', ''));  % Remove underscores
    else
        expName = 'Unknown';
    end
    if ~ismember(expName, uniqueExpNames)
        uniqueExpNames{end+1} = expName;
    end
    expIdx = find(strcmp(uniqueExpNames, expName));  % Get experiment index
    % Identify the strategy group for this experiment
    groupIdx = [];
    for k = 1:numGroups
        if contains(fname, upper(strategyGroups{k}))
            groupIdx = k;
            fprintf('  -> Matched group: %s\n', groupNames{k});
            break;
        end
    end
    if isempty(groupIdx), continue; end
    % Get the data columns
    M1Error = T.M1_new_local_err;
    M2Error = T.M2_new_local_err;
    M3Error = T.global_error;
    M1_neigh_diff = T.M1_neigh_diff;
    M2_neigh_diff = T.M2_neigh_diff;
    local_Wb_M1 = T.local_Wb_M1;
    local_Wb_M2 = T.local_Wb_M2;
    final_Wb_M1 = T.final_Wb_M1;
    final_Wb_M2 = T.final_Wb_M2;
    % ====== UPDATED Φ LOGIC BASED ON STRATEGY ======
    if groupIdx == 1  % LOCAL (M1, M2, M3)
        signal = [M1Error, M2Error, M3Error];
        phi = temporal_mutual_information(signal);
    elseif groupIdx == 2  % NEIGHBOUR
        signal = [M1Error, M2Error, M1_neigh_diff, M2_neigh_diff, M3Error];
        phi = temporal_mutual_information(signal);
    elseif groupIdx == 3  % SELFISH
        signal = [M1Error, M2Error, M3Error];
        if any(local_Wb_M1 > 0)
            signal = [signal, M1_neigh_diff];
        end
        if any(local_Wb_M2 > 0)
            signal = [signal, M2_neigh_diff];
        end
        phi = temporal_mutual_information(signal);
    elseif groupIdx == 4  % GLOBAL ONLY
        phi = temporal_mutual_information(M3Error);
    elseif groupIdx == 5  % GLOBAL
        signal = [M1Error, M2Error, M3Error];
        if any(final_Wb_M1 > 0)
            signal = [signal, M1_neigh_diff];
        end
        if any(final_Wb_M2 > 0)
            signal = [signal, M2_neigh_diff];
        end
        phi = temporal_mutual_information(signal);
      elseif groupIdx == 6  % HOMEO
        signal = [M1Error, M2Error, M3Error];
        if any(final_Wb_M1 > 0)
            signal = [signal, M1_neigh_diff];
        end
        if any(final_Wb_M2 > 0)
            signal = [signal, M2_neigh_diff];
        end
        phi = temporal_mutual_information(signal);
    end
    % Store Φ in the phiMatrix
    phiMatrix(groupIdx, expIdx) = phi;
    fprintf('  -> Stored Φ = %.4f in group %d, EXP = %s (column %d)\n', ...
        phi, groupIdx, expName, expIdx);
end
% Trim unused experiment columns
numExps = length(uniqueExpNames);
phiMatrix = phiMatrix(:, 1:numExps);
% ===== Plotting
figure;
hold on;
x = 1:numGroups;  % x-axis positions for groups
bh = bar(x, phiMatrix, 0.5, 'stacked');  % Create a stacked bar chart
for i = 1:numExps
    bh(i).FaceColor = expColors(i, :);  % Assign color for each experiment
end
ylabel('Integrated Information (Φ)');
set(gca, 'XTick', 1:numGroups, 'XTickLabel', groupNames);  % Label the x-axis with group names
xlabel('Control Strategy');
title('4/4: Integrated Information (Φ)');
grid on;
% Add the legend with experiment names
legend(uniqueExpNames, 'Location', 'northeastoutside');
function phi = temporal_mutual_information(signal)
    % Remove rows with any NaNs
    signal = signal(~any(isnan(signal), 2), :);
    % Not enough data
    if size(signal, 1) < 2
        phi = NaN;
        return;
    end
    % Z-score normalization
    signal = (signal - mean(signal)) ./ std(signal);
    % Create time-lagged vectors
    x_t = signal(1:end-1, :);
    x_t1 = signal(2:end, :);
    % Correlation matrix
    R = corr(x_t, x_t1, 'rows', 'complete');
    % Check if correlation matrix is valid
    if any(isnan(R(:)))
        phi = 0;
        return;
    end
    % ✅ Correct: elementwise square
    R_squared = R .* R;
    % I - R^2
    I_minus_R2 = eye(size(R)) - R_squared;
    % Compute determinant safely
    detVal = det(I_minus_R2);
    if detVal <= 0 || ~isfinite(detVal)
        phi = 0;
    else
        phi = 0.5 * log(1 / detVal);
    end
end


