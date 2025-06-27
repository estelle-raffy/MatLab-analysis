clear; clc; close all;

% Folder with your Excel files
dataFolder = 'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons_NOGRAV\Wa3Wb3Wc3\target_air_500\NoBall';
files = dir(fullfile(dataFolder, '*.xlsx'));

%{
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
%}

% Define strategy groups
strategyGroups = {'LOCAL','NEIGHBOUR','SELFISH', 'GLOBAL_ONLY', 'GLOBAL', 'HOMEO'};
groupNames = {'Local', 'Neighbour', 'Selfish', 'M3 only', 'Global Mod', 'Homeo'};
%strategyGroups = {'LOCAL', 'NEIGHBOUR_ONLY','NEIGHBOUR','SELFISH', 'GLOBAL_ONLY', 'GLOBAL', 'HOMEO'};
%groupNames = {'Local', 'Neigh only', 'Neighbour', 'Selfish', 'M3 only', 'Global Mod', 'Homeo'};
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
%% ======= PLOT 1.a/5: Final Shapes Grouped by Strategy (Updated, Label Underneath)
% ===============================================
figure('Name','1.a/5: Final Shapes Grouped by Strategy');
hold on;
axis equal;
xlabel('X');
ylabel('Y');
title('1.a/5: Diversity of Solutions by Strategy');
grid on;
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);

spacing_x = 5;
spacing_y = 3;
max_per_column = 8;
groupSolutionCount = zeros(1, numGroups);
theta = linspace(0, 2*pi, 20);
radius_body = 0.05;
expColors = lines(100);
uniqueExpNames = {};

for i = 1:length(files)
    T = readtable(fullfile(files(i).folder, files(i).name), 'Sheet', 'Tabelle1');
    
    fname = upper(files(i).name);

    expMatch = regexp(fname, '(EXP[_\s]?0*\d+)', 'tokens', 'ignorecase');
    if ~isempty(expMatch)
        expName = upper(strrep(expMatch{1}{1}, '_', ''));
    else
        expName = 'Unknown';
    end
    if ~ismember(expName, uniqueExpNames)
        uniqueExpNames{end+1} = expName;
    end
    expIdx = find(strcmp(uniqueExpNames, expName));
    color = expColors(expIdx, :);

    % ==== New Logic: Find the longest ≥5s of all-zero actuation ====
    m1 = T.M1_actuation_final;
    m2 = T.M2_actuation_final;
    m3 = T.M3_actuation_final;
    t  = T.current_time;

    isZeroAll = (m1 == 0) & (m2 == 0) & (m3 == 0);
    d = diff([0; isZeroAll; 0]);
    starts = find(d == 1);
    ends = find(d == -1) - 1;

    longestDuration = 0;
    zeroEndIndex = NaN;

    for z = 1:length(starts)
        duration = t(ends(z)) - t(starts(z));
        if duration >= 5 && duration > longestDuration
            longestDuration = duration;
            zeroEndIndex = ends(z); % Use end of longest qualifying zero-streak
        end
    end

    foundStreak = ~isnan(zeroEndIndex);

    if ~foundStreak
        % ==== Plot grey cross and label if no valid streak ====
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

        idx = groupSolutionCount(groupIdx);
        dx = (groupIdx - 1) * spacing_x;
        dy = -mod(idx, max_per_column) * spacing_y;
        dx = dx + floor(idx / max_per_column) * (spacing_x/2);
        groupSolutionCount(groupIdx) = groupSolutionCount(groupIdx) + 1;

        label_x = dx + 1;
        label_y = dy;
        plot(label_x, label_y, 'x', 'MarkerSize', 8, 'LineWidth', 2, 'Color', [0.7, 0.7, 0.7]);
        text(label_x, label_y - 0.3, expName, 'FontSize', 8, 'Color', [0.7, 0.7, 0.7], ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontWeight', 'bold');
        continue;
    end

    % ==== Extract body positions at the end of the longest zero-streak ====
    bx = [T.body1_pos_x(zeroEndIndex), T.body2_pos_x(zeroEndIndex), T.body3_pos_x(zeroEndIndex), ...
          T.body4_pos_x(zeroEndIndex), T.body5_pos_x(zeroEndIndex)];
    by = [T.body1_pos_y(zeroEndIndex), T.body2_pos_y(zeroEndIndex), T.body3_pos_y(zeroEndIndex), ...
          T.body4_pos_y(zeroEndIndex), T.body5_pos_y(zeroEndIndex)];

    % === Improved Normalization ===
    bx = bx - bx(1);
    bx = bx - min(bx);
    range_bx = max(bx);
    if range_bx > 0
        bx = bx / range_bx * 2;
    end
    
    by = by - by(1);
    by = by - min(by);
    range_by = max(by);
    if range_by > 0
        by = by / range_by * 1;
    end

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

    idx = groupSolutionCount(groupIdx);
    dx = (groupIdx - 1) * spacing_x;
    dy = -mod(idx, max_per_column) * spacing_y;
    dx = dx + floor(idx / max_per_column) * (spacing_x/2);
    groupSolutionCount(groupIdx) = groupSolutionCount(groupIdx) + 1;

    bx = bx + dx;
    by = by + dy;

    for j = 1:size(spring_pairs, 1)
        i1 = spring_pairs(j,1);
        i2 = spring_pairs(j,2);
        plot([bx(i1), bx(i2)], [by(i1), by(i2)], '-', 'Color', color, 'LineWidth', 1.5);
    end

    for j = 1:5
        fill(bx(j) + radius_body*cos(theta), by(j) + radius_body*sin(theta), ...
            color, 'FaceAlpha', 0.6, 'EdgeColor', 'k', 'LineWidth', 0.5);
    end

    label_x = mean(bx);
    label_y = min(by) - 0.3;
    text(label_x, label_y, expName, 'FontSize', 8, 'Color', color, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontWeight', 'bold');
end

for k = 1:numGroups
    xpos = (k-1) * spacing_x + 1;
    ypos = 2;
    text(xpos, ypos, groupNames{k}, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 12);
end

xlim([-1, numGroups * spacing_x]);
ylim([-max_per_column * spacing_y - 1, 3]);



% ===============================================
%% ======= PLOT 1.b/5: All Final Shapes Grouped by Strategy (Grey if failed)
% ===============================================
figure('Name','1.b/5: All Final Shapes by Strategy (Success = Colour, Failed = Grey)');
hold on;
axis equal;
xlabel('X');
ylabel('Y');
title('1.b/5: All Final Solutions by Strategy (Success = Colour, Failed = Grey)');
grid on;
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);

spacing_x = 5;
spacing_y = 3;
max_per_column = 8;
groupSolutionCount = zeros(1, numGroups);
theta = linspace(0, 2*pi, 20);
radius_body = 0.05;
expColors = lines(100);
uniqueExpNames = {};

for i = 1:length(files)
    T = readtable(fullfile(files(i).folder, files(i).name), 'Sheet', 'Tabelle1');
    fname = upper(files(i).name);

    % Match experiment name
    expMatch = regexp(fname, '(EXP[_\s]?0*\d+)', 'tokens', 'ignorecase');
    if ~isempty(expMatch)
        expName = upper(strrep(expMatch{1}{1}, '_', ''));
    else
        expName = 'Unknown';
    end
    if ~ismember(expName, uniqueExpNames)
        uniqueExpNames{end+1} = expName;
    end
    expIdx = find(strcmp(uniqueExpNames, expName));
    color = expColors(expIdx, :);

    % ==== Updated Logic: Find longest ≥5s of all-zero actuation ====
    m1 = T.M1_actuation_final;
    m2 = T.M2_actuation_final;
    m3 = T.M3_actuation_final;
    t  = T.current_time;

    isZeroAll = (m1 == 0) & (m2 == 0) & (m3 == 0);
    d = diff([0; isZeroAll; 0]);
    starts = find(d == 1);
    ends = find(d == -1) - 1;

    longestDuration = 0;
    zeroEndIndex = NaN;

    for z = 1:length(starts)
        duration = t(ends(z)) - t(starts(z));
        if duration >= 5 && duration > longestDuration
            longestDuration = duration;
            zeroEndIndex = ends(z);
        end
    end

    foundStreak = ~isnan(zeroEndIndex);

    % Determine strategy group
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

    if foundStreak
        % ==== Use end of the longest valid zero-streak ====
        bx = [T.body1_pos_x(zeroEndIndex), T.body2_pos_x(zeroEndIndex), T.body3_pos_x(zeroEndIndex), ...
              T.body4_pos_x(zeroEndIndex), T.body5_pos_x(zeroEndIndex)];
        by = [T.body1_pos_y(zeroEndIndex), T.body2_pos_y(zeroEndIndex), T.body3_pos_y(zeroEndIndex), ...
              T.body4_pos_y(zeroEndIndex), T.body5_pos_y(zeroEndIndex)];

        % ==== Convert back to Python coordinates ====
        bodyX_current_pos_x = bx + 300;
        bodyX_current_pos_y = (-1) * by + 550;
        
        % ==== Display in MATLAB command window ====
        fprintf('--- SUCCESS: %s (Group: %s) ---\n', expName, groupNames{groupIdx});
        for b = 1:5
            fprintf('Body %d: (%.2f, %.2f)\n', b, bodyX_current_pos_x(b), bodyX_current_pos_y(b));
        end

    else
        % ==== Failed: fallback to final timestep, use grey color ====
        bx = [T.body1_pos_x(end), T.body2_pos_x(end), T.body3_pos_x(end), ...
              T.body4_pos_x(end), T.body5_pos_x(end)];
        by = [T.body1_pos_y(end), T.body2_pos_y(end), T.body3_pos_y(end), ...
              T.body4_pos_y(end), T.body5_pos_y(end)];
        color = [0.7, 0.7, 0.7];
    end

    % === Improved Normalization ===
    bx = bx - bx(1);
    bx = bx - min(bx);
    range_bx = max(bx);
    if range_bx > 0
        bx = bx / range_bx * 2;
    end

    by = by - by(1);
    by = by - min(by);
    range_by = max(by);
    if range_by > 0
        by = by / range_by * 1;
    end

    % Apply offset
    bx = bx + dx;
    by = by + dy;

    % Plot springs
    for j = 1:size(spring_pairs, 1)
        i1 = spring_pairs(j,1);
        i2 = spring_pairs(j,2);
        plot([bx(i1), bx(i2)], [by(i1), by(i2)], '-', 'Color', color, 'LineWidth', 1.5);
    end

    % Plot bodies
    for j = 1:5
        fill(bx(j) + radius_body*cos(theta), by(j) + radius_body*sin(theta), ...
            color, 'FaceAlpha', 0.6, 'EdgeColor', 'k', 'LineWidth', 0.5);
    end

    % Label
    label_x = mean(bx);
    label_y = min(by) - 0.3;
    text(label_x, label_y, expName, 'FontSize', 8, 'Color', color, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontWeight', 'bold');
end

% Group names
for k = 1:numGroups
    xpos = (k-1) * spacing_x + 1;
    ypos = 2;
    text(xpos, ypos, groupNames{k}, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 12);
end

xlim([-1, numGroups * spacing_x]);
ylim([-max_per_column * spacing_y - 1, 3]);


% ==================================================================
%% ======= PLOT 2/5: Stability (Subplots per EXP)
% ==================================================================

expList = {'EXP1', 'EXP2', 'EXP3', 'EXP4'};
nExps = numel(expList);

% === Setup data containers ===
expData = struct();
for e = 1:nExps
    expData(e).starts = cell(1, numGroups);
    expData(e).durations = cell(1, numGroups);
    expData(e).colors = cell(1, numGroups);
end

expColors = lines(100);
uniqueExpNames = {};
strategyColors = lines(numGroups);

% === Loop through files ===
for i = 1:length(files)
    T = readtable(fullfile(files(i).folder, files(i).name), 'Sheet', 'Tabelle1');
    fname = upper(files(i).name);

    % Strategy group detection
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

    % Experiment name detection
    expMatch = regexp(fname, '(EXP[_\s]?0*\d+)', 'tokens', 'ignorecase');
    if ~isempty(expMatch)
        expName = upper(strrep(expMatch{1}{1}, '_', ''));
    else
        expName = sprintf('FILE%d', i);
    end

    expIdx = find(strcmp(expList, expName));
    if isempty(expIdx)
        warning('File %s does not match any known experiment.', files(i).name);
        continue;
    end

    if ~ismember(expName, uniqueExpNames)
        uniqueExpNames{end+1} = expName;
    end
    expColorIdx = find(strcmp(uniqueExpNames, expName));
    expColor = expColors(expColorIdx, :);

    % Check for simultaneous zero actuation
    m1 = T.M1_actuation_final;
    m2 = T.M2_actuation_final;
    m3 = T.M3_actuation_final;
    t = T.current_time;

    isZeroAll = (m1 == 0) & (m2 == 0) & (m3 == 0);
    d = diff([0; isZeroAll; 0]);
    starts = find(d == 1);
    ends = find(d == -1) - 1;

    if isempty(starts)
        continue;
    end

    durations = t(ends) - t(starts);
    onsetTimes = t(starts);

    for idx = 1:length(onsetTimes)
        expData(expIdx).starts{groupIdx}(end+1) = onsetTimes(idx);
        expData(expIdx).durations{groupIdx}(end+1) = durations(idx);
        expData(expIdx).colors{groupIdx}(end+1, :) = expColor;
    end
end

% === Plotting ===
figure('Name','Actuation Silence Periods (All Zeros, by EXP)');
groupSpacing = 5;

for e = 1:nExps
    subplot(2, 2, e);
    hold on;
    title(['Periods Where M1, M2, M3 = 0 — ' expList{e}]);
    ylabel('Start Time (s)');
    xlabel('Strategy Group');
    grid on;

    x_ticks = [];
    x_labels = {};

    for k = 1:numGroups
        starts = expData(e).starts{k};
        durations = expData(e).durations{k};
        colors = expData(e).colors{k};

        n = numel(starts);
        if n == 0
            x_pos = (k - 1) * groupSpacing + 1;
            x_ticks(end+1) = x_pos;
            x_labels{end+1} = groupNames{k};
            continue;
        end

        x = (k - 1) * groupSpacing + 1;
        x_ticks(end+1) = x;
        x_labels{end+1} = groupNames{k};

        for j = 1:n
            startTime = starts(j);
            duration = durations(j);
            endTime = startTime + duration;
            colorToUse = colors(j,:);

            % Vertical bar (error-bar style)
            line([x, x], [startTime, endTime], ...
                'Color', colorToUse, 'LineWidth', 1.5);

            % Top & bottom caps
            capWidth = 0.2;
            line([x - capWidth, x + capWidth], [startTime, startTime], ...
                'Color', colorToUse, 'LineWidth', 1.2);
            line([x - capWidth, x + capWidth], [endTime, endTime], ...
                'Color', colorToUse, 'LineWidth', 1.2);

            % No center dot (was removed)
        end
    end

    % Format axis
    xticks(sort(unique(x_ticks)));
    xticklabels(x_labels);
    ylim([0, 65]);
    xlim([0, max(x_ticks) + groupSpacing]);
end

% ========================================================================
%% ======= PLOT 3/5: Successful Shapes Comparisons Within/Between Files
% ========================================================================

% === Input file ===
file = 'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons\Wa1Wb1Wc0.33\NoBall\EXP2_NEIGHBOUR.xlsx'; % <-- Update this
T = readtable(file);

% === Extract strategy name from filename ===
[~, filename, ~] = fileparts(file);
tokens = regexp(filename, 'EXP\d+_(\w+)', 'tokens');
if ~isempty(tokens)
    strategyName = tokens{1}{1};
else
    strategyName = 'Unknown Strategy';
end

% === Define spring connections and drawing parameters ===
spring_pairs = [1 3; 3 2; 1 2; 2 5; 3 5; 5 4; 2 4];
theta = linspace(0, 2*pi, 20);
radius_body = 0.05;

% === Extract actuation data ===
m1 = T.M1_actuation_final;
m2 = T.M2_actuation_final;
m3 = T.M3_actuation_final;
t  = T.current_time;

isZeroAll = (m1 == 0) & (m2 == 0) & (m3 == 0);
d = diff([0; isZeroAll; 0]);
starts = find(d == 1);
ends = find(d == -1) - 1;

% === Setup figure ===
figure('Name', ['Intra-file shape diversity — ' strategyName]);
hold on; axis equal; grid on;
title(['All Final Positions (≥5s 0-actuation) — Strategy: ' strategyName]);
xlabel('X'); ylabel('Y');
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
colors = lines(length(starts));

% === Plot each qualifying shape snapshot ===
for z = 1:length(starts)
    duration = t(ends(z)) - t(starts(z));
    if duration < 5, continue; end

    idx = ends(z);
    bx = [T.body1_pos_x(idx), T.body2_pos_x(idx), T.body3_pos_x(idx), ...
          T.body4_pos_x(idx), T.body5_pos_x(idx)];
    by = [T.body1_pos_y(idx), T.body2_pos_y(idx), T.body3_pos_y(idx), ...
          T.body4_pos_y(idx), T.body5_pos_y(idx)];

    % Normalize positions
    bx = bx - bx(1); bx = bx - min(bx); if max(bx) > 0, bx = bx / max(bx) * 2; end
    by = by - by(1); by = by - min(by); if max(by) > 0, by = by / max(by); end

    % Draw springs
    for j = 1:size(spring_pairs,1)
        plot([bx(spring_pairs(j,1)), bx(spring_pairs(j,2))], ...
             [by(spring_pairs(j,1)), by(spring_pairs(j,2))], ...
             '-', 'Color', colors(z,:), 'LineWidth', 1.2);
    end

    % Draw bodies
    for j = 1:5
        fill(bx(j) + radius_body*cos(theta), ...
             by(j) + radius_body*sin(theta), ...
             colors(z,:), 'FaceAlpha', 0.5, ...
             'EdgeColor', 'k', 'LineWidth', 0.3);
    end
end

% === Add legend for strategy name ===
%legend(strategyName, 'Location', 'northeastoutside');


%{
%% MAY NEED COMMENTING OFF IF WANT ACTUATION SIGNAL TO WORK
% === Plot 3.b: Inter-File Comparison of Best Shapes ===

files = {
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons\Wa1Wb1Wc0.33\EXP1_NEIGHBOUR.xlsx',
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons\Wa3Wb3Wc3\EXP2_NEIGHBOUR.xlsx',
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons\Wa1Wb1Wc3\EXP2_NEIGHBOUR.xlsx',
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons_NOGRAV\Wa1Wb1Wc3\target_air_500\EXP3_NEIGHBOUR.xlsx',
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons_NOGRAV\Wa3Wb3Wc3\target_air_500\EXP3_NEIGHBOUR.xlsx',
};

spring_pairs = [1 3; 3 2; 1 2; 2 5; 3 5; 5 4; 2 4];
theta = linspace(0, 2*pi, 20);
radius_body = 0.05;

figure('Name','Inter-file shape comparison');
hold on; axis equal; grid on;
title('Best Final Position Across Files');
xlabel('X'); ylabel('Y');
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
colors = lines(length(files));

legend_entries = {};
plot_handles = [];

for f = 1:length(files)
    T = readtable(files{f});
    m1 = T.M1_actuation_final;
    m2 = T.M2_actuation_final;
    m3 = T.M3_actuation_final;
    t  = T.current_time;

    isZeroAll = (m1 == 0) & (m2 == 0) & (m3 == 0);
    d = diff([0; isZeroAll; 0]);
    starts = find(d == 1);
    ends = find(d == -1) - 1;

    longest = 0; idx = NaN;
    for z = 1:length(starts)
        duration = t(ends(z)) - t(starts(z));
        if duration >= 5 && duration > longest
            longest = duration;
            idx = ends(z);
        end
    end
    if isnan(idx)
        continue;
    end

    % Assign color only for valid entries
    color = colors(f, :);

    bx = [T.body1_pos_x(idx), T.body2_pos_x(idx), T.body3_pos_x(idx), ...
          T.body4_pos_x(idx), T.body5_pos_x(idx)];
    by = [T.body1_pos_y(idx), T.body2_pos_y(idx), T.body3_pos_y(idx), ...
          T.body4_pos_y(idx), T.body5_pos_y(idx)];

    % Normalize positions
    bx = bx - bx(1); bx = bx - min(bx); if max(bx) > 0, bx = bx / max(bx) * 2; end
    by = by - by(1); by = by - min(by); if max(by) > 0, by = by / max(by); end

    % Plot springs and store handle for legend
    h = plot([bx(spring_pairs(1,1)), bx(spring_pairs(1,2))], ...
             [by(spring_pairs(1,1)), by(spring_pairs(1,2))], ...
             '-', 'Color', color, 'LineWidth', 1.5);
    plot_handles(end+1) = h; %#ok<*SAGROW>

    for j = 2:size(spring_pairs,1)
        plot([bx(spring_pairs(j,1)), bx(spring_pairs(j,2))], ...
             [by(spring_pairs(j,1)), by(spring_pairs(j,2))], ...
             '-', 'Color', color, 'LineWidth', 1.5);
    end

    for j = 1:5
        fill(bx(j) + radius_body*cos(theta), by(j) + radius_body*sin(theta), ...
             color, 'FaceAlpha', 0.6, 'EdgeColor', 'k', 'LineWidth', 0.3);
    end

    % Extract experiment type from filename (case-insensitive, more specific first!)
    file_lower = lower(files{f});
    if contains(file_lower, 'global_only')
        legend_entries{end+1} = 'GLOBAL_ONLY';
    elseif contains(file_lower, 'global')
        legend_entries{end+1} = 'GLOBAL';
    elseif contains(file_lower, 'local')
        legend_entries{end+1} = 'LOCAL';
    elseif contains(file_lower, 'neighbour')
        legend_entries{end+1} = 'NEIGHBOUR';
    elseif contains(file_lower, 'selfish')
        legend_entries{end+1} = 'SELFISH';
    elseif contains(file_lower, 'homeo')
        legend_entries{end+1} = 'HOMEO';
    else
        legend_entries{end+1} = 'UNKNOWN';
    end
end

% Only show legend for valid plotted entries
legend(plot_handles, legend_entries, 'Location', 'best');


%}

% ==================================================================
%% ===== Plot 4/5 Actuation Signals Over Time (By EXP, One Strategy)
% ==================================================================

% === Define strategy group to display === CHECH FILE NAME TOP OF CODE
strategyToPlot = 'LOCAL';  % <<== input desired strategy 

% === Experiment names ===
expList = {'EXP1', 'EXP2', 'EXP3', 'EXP4'};
nExps = numel(expList);

% === Setup plot ===
figure('Name', ['Actuation Signals — Strategy: ' strategyToPlot]);
tiledlayout(2,2);

% === Colors for M1, M2, M3 ===
colorM1 = [0.8 0 0];    % red
colorM2 = [0 0.6 0];    % green
colorM3 = [0 0.2 0.8];  % blue

for e = 1:nExps
    expName = expList{e};
    
    nexttile;
    hold on;
    title([expName ' — ' strategyToPlot]);
    xlabel('Time (s)');
    ylabel('Actuation Signal');
    grid on;

    foundData = false;

    % === Loop through files ===
    for i = 1:length(files)
        fname = upper(files(i).name);

        % Check if this file belongs to the current EXP and strategy group
        pattern = [expName '_' upper(strategyToPlot) '.XLSX'];
        if strcmp(fname, pattern)

            % Read data
            T = readtable(fullfile(files(i).folder, files(i).name), 'Sheet', 'Tabelle1');
            t = T.current_time;

            m1 = T.M1_actuation_final;
            m2 = T.M2_actuation_final;
            m3 = T.M3_actuation_final;

            % Plot actuation signals
            plot(t, m1, 'Color', colorM1, 'LineWidth', 1.2);
            plot(t, m2, 'Color', colorM2, 'LineWidth', 1.2);
            plot(t, m3, 'Color', colorM3, 'LineWidth', 1.2);

            foundData = true;
        end
    end

    if ~foundData
        text(0.5, 0.5, 'No Data Found', ...
            'Units', 'normalized', 'HorizontalAlignment', 'center', ...
            'FontSize', 12, 'FontWeight', 'bold', 'Color', [0.5 0 0]);
    else
        legend({'M1', 'M2', 'M3'}, 'Location', 'northeast');
    end
end


% ===============================================================================
%% ======= PLOTS 5/5 INTERNAL DYNAMICS : Selfishness vs Control strategies ======
% ===============================================================================

% File paths (replace these with your actual paths)
filePaths = {
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons\Wa1Wb3Wc1\EXP1_LOCAL.xlsx',
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons\Wa1Wb3Wc1\EXP1_NEIGHBOUR.xlsx'
};

% Strategy keywords
strategies = {'NEIGHBOUR', 'SELFISH', 'GLOBAL', 'HOMEO'};

% Define consistent colors
color_localM1 = [0.4, 0.76, 0.99];  % Light blue
color_finalM1 = [0.0, 0.2, 0.8];    % Dark blue
color_localM2 = [1.0, 0.6, 0.0];    % Orange
color_finalM2 = [0.8, 0.1, 0.1];    % Red

for i = 1:length(filePaths)
    % Load data
    data = readtable(filePaths{i});
    
    % Extract relevant columns
    time      = data.current_time;
    localM1   = data.local_Wb_M1;
    localM2   = data.local_Wb_M2;
    finalM1   = data.final_Wb_M1;
    finalM2   = data.final_Wb_M2;
    
    % Unique time steps
    unique_times = unique(time);
    
    % Preallocate vectors
    mean_localM1 = zeros(size(unique_times));
    mean_localM2 = zeros(size(unique_times));
    mean_finalM1 = zeros(size(unique_times));
    mean_finalM2 = zeros(size(unique_times));
    
    % Compute mean values per time point
    for j = 1:length(unique_times)
        t = unique_times(j);
        idx = time == t;
        
        mean_localM1(j) = mean(localM1(idx));
        mean_localM2(j) = mean(localM2(idx));
        mean_finalM1(j) = mean(finalM1(idx));
        mean_finalM2(j) = mean(finalM2(idx));
    end
    
    % Create figure for this strategy
    fig = figure('Name', ['Wb Curves - File ' num2str(i)], 'NumberTitle', 'off');
    t = tiledlayout(2,1);
    
    % Extract strategy name from filename
    [~, fname, ~] = fileparts(filePaths{i});
    strategy_name = 'Unknown Strategy';
    for s = 1:length(strategies)
        if ~isempty(regexp(upper(fname), ['(^|[_-])' strategies{s} '($|[_-])'], 'once'))
            strategy_name = strategies{s};
            break;
        end
    end
    t.Title.String = ['Strategy: ' strategy_name];
    
    % ---- SUBPLOT 1: M1 ----
    nexttile;
    plot(unique_times, mean_localM1, '-', 'Color', color_localM1, 'LineWidth', 1.8, 'DisplayName', 'Local Wb M1');
    hold on;
    plot(unique_times, mean_finalM1, '-', 'Color', color_finalM1, 'LineWidth', 1.8, 'DisplayName', 'Final Wb M1');
    hold off;
    ylabel('Wb M1');
    legend('Location', 'best');
    grid on;

    % ---- SUBPLOT 2: M2 ----
    nexttile;
    plot(unique_times, mean_localM2, '-', 'Color', color_localM2, 'LineWidth', 1.8, 'DisplayName', 'Local Wb M2');
    hold on;
    plot(unique_times, mean_finalM2, '-', 'Color', color_finalM2, 'LineWidth', 1.8, 'DisplayName', 'Final Wb M2');
    hold off;
    xlabel('Time');
    ylabel('Wb M2');
    legend('Location', 'best');
    grid on;
end


% ========================================================================
%% PLOT 6/7: Intra-Strategy 
% ========================================================================

file = 'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons\Wa1Wb3Wc1\EXP1_LOCAL.xlsx';
T = readtable(file);

spring_pairs = [1 3; 3 2; 1 2; 2 5; 3 5; 5 4; 2 4];
theta = linspace(0, 2*pi, 20);
radius_body = 0.05;
shapes = {};

m1 = T.M1_actuation_final;
m2 = T.M2_actuation_final;
m3 = T.M3_actuation_final;
t  = T.current_time;

isZeroAll = (m1 == 0) & (m2 == 0) & (m3 == 0);
d = diff([0; isZeroAll; 0]);
starts = find(d == 1);
ends = find(d == -1) - 1;

for z = 1:length(starts)
    duration = t(ends(z)) - t(starts(z));
    if duration < 5, continue; end

    idx = ends(z);
    bx = [T.body1_pos_x(idx), T.body2_pos_x(idx), T.body3_pos_x(idx), T.body4_pos_x(idx), T.body5_pos_x(idx)];
    by = [T.body1_pos_y(idx), T.body2_pos_y(idx), T.body3_pos_y(idx), T.body4_pos_y(idx), T.body5_pos_y(idx)];

    % Normalize
    bx = bx - bx(1); bx = bx - min(bx); if max(bx) > 0, bx = bx / max(bx) * 2; end
    by = by - by(1); by = by - min(by); if max(by) > 0, by = by / max(by); end

    shapes{end+1} = [bx(:), by(:)];
end

% Compute pairwise Procrustes distances: for each pair of shapes get d
% the smaller d, the more similar the shapes (0=identical shapes)
n = numel(shapes);
D = zeros(n);
for i = 1:n
    for j = i+1:n
        [d, ~] = procrustes(shapes{i}, shapes{j});
        D(i,j) = d; % store distance in Matrix D
        D(j,i) = d;
        fprintf('INTRA - Procrustes distance between shape %d and %d: %.4f\n', i, j, d);
    end
end

% Cluster shapes: start with one shape, group with next best similar 
% and so on, cutoff at a certain 'is different threshold'
Z = linkage(squareform(D), 'average'); % group similar shapes, linkage looks at distance matrix and groups closest first, makes a 'tree'
T_clust = cluster(Z, 'cutoff', 0.01, 'criterion', 'distance'); % T_clust tells what cluster belongs to: look at tree and use cutoff to sort into clusters
% shape below it get together, above belong to diff clusters 

% Mean shape per cluster: take shape as reference, make shape fit as
% possible (rotate/resize using Procruste) and get the mean of all of them. 
unique_clusters = unique(T_clust);
mean_shapes = cell(length(unique_clusters), 1);
for k = 1:length(unique_clusters)
    idxs = find(T_clust == unique_clusters(k));
    ref = shapes{idxs(1)};
    aligned = zeros(length(idxs), size(ref,1), size(ref,2)); % save shapes in aligned array
    for i = 1:length(idxs)
        [~, Z] = procrustes(ref, shapes{idxs(i)});
        aligned(i,:,:) = Z;
    end
    mean_shapes{k} = squeeze(mean(aligned, 1));
end

% Plot mean shapes
figure;
hold on;
axis equal;
grid on;

title('Morphological Clusters — Single File');
ylabel('Relative vertical position (normalized)');
set(gca, 'XTick', []);  % remove X-axis tick labels

colors = lines(length(mean_shapes));
for i = 1:length(mean_shapes)
    shape = mean_shapes{i};
    bx = shape(:,1) + i * 3;  % space shapes out
    by = shape(:,2);
    for j = 1:size(spring_pairs,1)
        plot([bx(spring_pairs(j,1)), bx(spring_pairs(j,2))], ...
             [by(spring_pairs(j,1)), by(spring_pairs(j,2))], '-', 'Color', colors(i,:), 'LineWidth', 1.5);
    end
    for j = 1:5
        fill(bx(j) + radius_body*cos(theta), by(j) + radius_body*sin(theta), ...
             colors(i,:), 'FaceAlpha', 0.6, 'EdgeColor', 'k', 'LineWidth', 0.3);
    end
end



% ========================================================================
%% PLOT 7/7: Inter-Strategy Shape Diversity Summary (by Mean Area)
% ========================================================================

files = {
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons_NOGRAV\Wa1Wb1Wc0.33\EXP1_GLOBAL.xlsx',
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons_NOGRAV\Wa1Wb1Wc3\EXP2_GLOBAL.xlsx',
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons_NOGRAV\Wa3Wb3Wc1\EXP3_GLOBAL.xlsx',
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons_NOGRAV\Wa3Wb3Wc3\EXP2_GLOBAL.xlsx',
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons_NOGRAV\Wa3Wb3Wc3\target_air_500\EXP3_GLOBAL.xlsx',
};

spring_pairs = [1 3; 3 2; 1 2; 2 5; 3 5; 5 4; 2 4];
theta = linspace(0, 2*pi, 20);
radius_body = 0.05;
shapes = {};

for f = 1:length(files)
    T = readtable(files{f});
    m1 = T.M1_actuation_final;
    m2 = T.M2_actuation_final;
    m3 = T.M3_actuation_final;
    t  = T.current_time;

    isZeroAll = (m1 == 0) & (m2 == 0) & (m3 == 0);
    d = diff([0; isZeroAll; 0]);
    starts = find(d == 1);
    ends = find(d == -1) - 1;

    for z = 1:length(starts)
        duration = t(ends(z)) - t(starts(z));
        if duration < 5, continue; end

        idx = ends(z);
        bx = [T.body1_pos_x(idx), T.body2_pos_x(idx), T.body3_pos_x(idx), T.body4_pos_x(idx), T.body5_pos_x(idx)];
        by = [T.body1_pos_y(idx), T.body2_pos_y(idx), T.body3_pos_y(idx), T.body4_pos_y(idx), T.body5_pos_y(idx)];

        % Normalize
        bx = bx - bx(1); bx = bx - min(bx); if max(bx) > 0, bx = bx / max(bx) * 2; end
        by = by - by(1); by = by - min(by); if max(by) > 0, by = by / max(by); end

        shapes{end+1} = [bx(:), by(:)];
    end
end

% Compute pairwise Procrustes distances
n = numel(shapes);
D = zeros(n);
for i = 1:n
    for j = i+1:n
        [d, ~] = procrustes(shapes{i}, shapes{j});
        D(i,j) = d;
        D(j,i) = d;
        fprintf('INTER - Procrustes distance between shape %d and %d: %.4f\n', i, j, d);
    end
end

% Cluster shapes
Z = linkage(squareform(D), 'average');
T_clust = cluster(Z, 'cutoff', 0.01, 'criterion', 'distance');

% Mean shape per cluster
unique_clusters = unique(T_clust);
mean_shapes = cell(length(unique_clusters), 1);
for k = 1:length(unique_clusters)
    idxs = find(T_clust == unique_clusters(k));
    ref = shapes{idxs(1)};
    aligned = zeros(length(idxs), size(ref,1), size(ref,2));
    for i = 1:length(idxs)
        [~, aligned_shape] = procrustes(ref, shapes{idxs(i)});
        aligned(i,:,:) = aligned_shape;
    end
    mean_shapes{k} = squeeze(mean(aligned, 1));
end

% Plot mean shapes
figure;
hold on;
axis equal;
grid on;

title('Morphological Clusters — Across File');
ylabel('Relative vertical position (normalized)');
set(gca, 'XTick', []);  % remove X-axis tick labels

colors = lines(length(mean_shapes));
for i = 1:length(mean_shapes)
    shape = mean_shapes{i};
    bx = shape(:,1) + i * 3;
    by = shape(:,2);
    for j = 1:size(spring_pairs,1)
        plot([bx(spring_pairs(j,1)), bx(spring_pairs(j,2))], ...
             [by(spring_pairs(j,1)), by(spring_pairs(j,2))], '-', 'Color', colors(i,:), 'LineWidth', 1.5);
    end
    for j = 1:5
        fill(bx(j) + radius_body*cos(theta), by(j) + radius_body*sin(theta), ...
             colors(i,:), 'FaceAlpha', 0.6, 'EdgeColor', 'k', 'LineWidth', 0.3);
    end
end


% === Dendrogram with shape thumbnails ===
figure;

% Top subplot: Dendrogram
subplot(2, 1, 1);
[H, T, outperm] = dendrogram(Z, 0);
set(H, 'LineWidth', 1.2);
cutoff_val = 0.01; % your cutoff value
yline(cutoff_val, '--r', sprintf('Cut-off = %.3f', cutoff_val), 'Color', 'r');


% After dendrogram plotting
hold on;

numMerges = size(Z, 1);
for i = 1:numMerges
    % Each merge corresponds to a branch between two clusters/leaves
    % Find the x coordinates of the branch ends from H(i)
    % H is an array of Line objects for the dendrogram branches
    % The two points connected by the merge are in the XData(2) and XData(3)
    xCoords = [H(i).XData(2), H(i).XData(3)];
    xCenter = mean(xCoords);   % middle point between the branch ends
    yHeight = Z(i, 3);         % height of this merge = Procrustes distance

    % Offset the text a bit above the branch height to avoid overlap
    text(xCenter, yHeight + 0.05 * range(ylim), sprintf('%.3f', yHeight), ...
         'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', [0 0 0]);
end
hold off;


title('Dendrogram of Shape Similarities');
ylabel('Procrustes Distance');
xlabel('');

% Bottom subplot: Shape thumbnails in dendrogram order
subplot(2, 1, 2);
hold on;
axis off;
axis equal;

nShapes = numel(outperm);
theta = linspace(0, 2*pi, 20);
radius = 0.05;

for i = 1:nShapes
    idx = outperm(i);
    shape = shapes{idx};

    % Normalize and center the shape
    shape = shape - mean(shape, 1);
    shape = shape / max(vecnorm(shape, 2, 2)) * 0.5;

    % Offset along x-axis to spread shapes out
    xOffset = (i - 1) * 1.2;
    bx = shape(:,1) + xOffset;
    by = shape(:,2);

    % Plot springs
    for j = 1:size(spring_pairs,1)
        plot([bx(spring_pairs(j,1)), bx(spring_pairs(j,2))], ...
             [by(spring_pairs(j,1)), by(spring_pairs(j,2))], ...
             '-', 'Color', [0.4 0.4 0.4]);
    end

    % Plot bodies
    for j = 1:5
        fill(bx(j) + radius*cos(theta), by(j) + radius*sin(theta), ...
             [0.5 0.5 0.5], 'EdgeColor', 'k');
    end

    % Add shape index label
    text(mean(bx), min(by) - 0.15, sprintf('%d', idx), ...
         'HorizontalAlignment', 'center', 'FontSize', 8);
end

xlim([-1, nShapes * 1.2]);
ylim([-1, 1]);
title('Shapes in Dendrogram Leaf Order');


%{

% --- Plot 2/2: Dendrogram with Distances and Shape Miniatures
figure;
[H, T, outperm] = dendrogram(Z, 0);  % '0' = full tree
set(H, 'LineWidth', 1.2);

% Add branch height distance labels
numMerges = size(Z, 1);
for i = 1:numMerges
    % Approximate coordinates for label
    x = mean(H(i).XData(2:3));
    y = H(i).YData(2);

    % Distance label
    dist = Z(i, 3);
    text(x, y + 0.005, sprintf('%.3f', dist), ...
         'HorizontalAlignment', 'center', ...
         'FontSize', 8, 'Color', [0.2 0.2 0.2]);
end

title('Dendrogram of Shape Similarities');
xlabel('Shape Index');
ylabel('Procrustes Distance');

% Red dashed cut-off line (re-added correctly)
cutoff = 0.01;
yline(cutoff, '--r', 'Cut-off');

hold on;

% === Plot shape miniatures below each leaf ===
leaf_positions = get(gca, 'XTick');
theta = linspace(0, 2*pi, 20);
radius_body = 0.05;
scale = 0.1;
y_offset = -0.07;  % move below x-axis

for i = 1:length(outperm)
    idx = outperm(i);  % original shape index
    shape = shapes{idx};

    % Center and scale
    shape = shape - mean(shape, 1);
    shape = shape * scale;

    % Offset
    bx = shape(:,1) + leaf_positions(i);
    by = shape(:,2) + y_offset;

    % Draw lines
    for j = 1:size(spring_pairs,1)
        plot([bx(spring_pairs(j,1)), bx(spring_pairs(j,2))], ...
             [by(spring_pairs(j,1)), by(spring_pairs(j,2))], ...
             '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.8);
    end

    % Draw bodies
    for j = 1:5
        fill(bx(j) + radius_body*scale*cos(theta), ...
             by(j) + radius_body*scale*sin(theta), ...
             [0.7 0.7 0.7], 'EdgeColor', 'k', 'LineWidth', 0.2);
    end
end

%}

%{
% ===================================================================
%% ======= PLOTS 5/5: 3D Scatter of Selfishness vs Weight Ratios ======
% ===================================================================

% === Setup ===
baseFolder = 'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons';
weightFolders = dir(baseFolder);
weightFolders = weightFolders([weightFolders.isdir]); 
weightFolders = weightFolders(~ismember({weightFolders.name}, {'.', '..'}));

dataStruct = struct('exp', {}, 'strategy', {}, 'WaWb', {}, 'WcWb', {}, ...
                    'selfishPercent', {}, 'folder', {}, 'Wb', {});  % Include Wb

for w = 1:length(weightFolders)
    folderName = weightFolders(w).name;
    weightPattern = regexp(folderName, 'Wa(\d+)Wb(\d+)Wc(\d+)', 'tokens');
    if isempty(weightPattern), continue; end

    weights = str2double(weightPattern{1});
    Wa = weights(1); Wb = weights(2); Wc = weights(3);
    if Wb == 0, continue; end  % avoid division by zero

    WaWb = Wa / Wb;
    WcWb = Wc / Wb;

    fullPath = fullfile(baseFolder, folderName);
    excelFiles = dir(fullfile(fullPath, '*.xlsx'));

    for f = 1:length(excelFiles)
        fileName = upper(excelFiles(f).name);
        expMatch = regexp(fileName, '(EXP\d+)_([A-Z]+)', 'tokens');
        if isempty(expMatch), continue; end

        expID = expMatch{1}{1};
        strategy = expMatch{1}{2};
        if ~ismember(strategy, {'SELFISH', 'GLOBAL', 'HOMEO'}), continue; end

        Tdata = readtable(fullfile(fullPath, excelFiles(f).name));
        if ~all(ismember({'final_Wb_M1', 'final_Wb_M2'}, Tdata.Properties.VariableNames)), continue; end

        avgSelfish = mean([Tdata.final_Wb_M1 == 0, Tdata.final_Wb_M2 == 0], 'all') * 100;

        % Store result
        dataStruct(end+1).exp = expID;
        dataStruct(end).strategy = strategy;
        dataStruct(end).WaWb = WaWb;
        dataStruct(end).WcWb = WcWb;
        dataStruct(end).selfishPercent = avgSelfish;
        dataStruct(end).folder = folderName;
        dataStruct(end).Wb = Wb;  % Now Wb is defined
    end
end

% === Convert and Add Wa/Wc ===
if isempty(dataStruct)
    error('No valid data found.');
end

T = struct2table(dataStruct);
T.WaWc = T.WaWb ./ T.WcWb;  % Compute Wa/Wc ratio

strategies = {'SELFISH', 'GLOBAL', 'HOMEO'};
shapes = {'o', 's', '^'};
colors = lines(numel(strategies));
uniqueExps = unique(T.exp);

% === Plot 5.a: 3D Selfishness VS Wa/Wb VS Wc/Wb ===
figure('Name', '3D Scatter: Selfishness vs Wa/Wb & Wc/Wb');
for i = 1:min(4, numel(uniqueExps))
    subplot(2, 2, i);
    hold on; grid on; view(45, 25);
    title(sprintf('Experiment %s', uniqueExps{i}));
    xlabel('Wa / Wb'); ylabel('Wc / Wb'); zlabel('% Time Being Selfish');

    for s = 1:length(strategies)
        stratData = T(strcmp(T.strategy, strategies{s}) & strcmp(T.exp, uniqueExps{i}), :);
        scatter3(stratData.WaWb, stratData.WcWb, stratData.selfishPercent, ...
            80, colors(s, :), shapes{s}, 'filled');
    end
    legend(strategies, 'Location', 'best');
end

% === Plot 5.a: 3D Selfishness VS Wa/Wc VS Wb ===
figure('Name', '3D Scatter: Selfishness vs Wa/Wc & Wb');
for i = 1:min(4, numel(uniqueExps))
    subplot(2, 2, i);
    hold on; grid on; view(45, 25);
    title(sprintf('Experiment %s', uniqueExps{i}));
    xlabel('Wa / Wc'); ylabel('Wb'); zlabel('% Time Being Selfish');

    for s = 1:length(strategies)
        stratData = T(strcmp(T.strategy, strategies{s}) & strcmp(T.exp, uniqueExps{i}), :);
        scatter3(stratData.WaWc, stratData.Wb, stratData.selfishPercent, ...
            80, colors(s, :), shapes{s}, 'filled');
    end
    legend(strategies, 'Location', 'best');
end

% === Plot 5.c: 2D Selfishness VS Wa/Wb ===
figure('Name', '2D Scatter: Selfishness vs Wa/Wb');
for i = 1:min(4, numel(uniqueExps))
    subplot(2, 2, i);
    hold on; grid on;
    title(sprintf('Experiment %s', uniqueExps{i}));
    xlabel('Wa / Wb'); ylabel('% Time Being Selfish');

    for s = 1:length(strategies)
        stratData = T(strcmp(T.strategy, strategies{s}) & strcmp(T.exp, uniqueExps{i}), :);
        scatter(stratData.WaWb, stratData.selfishPercent, ...
            80, shapes{s}, 'MarkerFaceColor', colors(s, :), 'MarkerEdgeColor', colors(s, :));
    end
    legend(strategies, 'Location', 'best');
end

% === Plot 5.d: 2D Selfishness VS Wc/Wb ===
figure('Name', '2D Scatter: Selfishness vs Wc/Wb');
for i = 1:min(4, numel(uniqueExps))
    subplot(2, 2, i);
    hold on; grid on;
    title(sprintf('Experiment %s', uniqueExps{i}));
    xlabel('Wc / Wb'); ylabel('% Time Being Selfish');

    for s = 1:length(strategies)
        stratData = T(strcmp(T.strategy, strategies{s}) & strcmp(T.exp, uniqueExps{i}), :);
        scatter(stratData.WcWb, stratData.selfishPercent, ...
            80, shapes{s}, 'MarkerFaceColor', colors(s, :), 'MarkerEdgeColor', colors(s, :));
    end
    legend(strategies, 'Location', 'best');
end
%}

%{
% ===================================================================
%% ======= PLOT 3/3: beingSelfish dynamics under weight conditions 
% ===================================================================

% === Setup ===
baseFolder = 'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons';
weightFolders = dir(baseFolder);
weightFolders = weightFolders([weightFolders.isdir]); 
weightFolders = weightFolders(~ismember({weightFolders.name}, {'.', '..'}));

dataStruct = struct('exp', {}, 'strategy', {}, 'weightRatio', {}, 'selfishPercent', {}, 'folder', {});

for w = 1:length(weightFolders)
    folderName = weightFolders(w).name;
    weightPattern = regexp(folderName, 'Wa(\d+)Wb(\d+)Wc(\d+)', 'tokens');
    if isempty(weightPattern), continue; end

    weights = str2double(weightPattern{1});
    Wa = weights(1); Wb = weights(2); Wc = weights(3);
    if Wc == 0, continue; end  % avoid division by zero
    weightRatio = (Wa + Wb) / Wc;

    fullPath = fullfile(baseFolder, folderName);
    excelFiles = dir(fullfile(fullPath, '*.xlsx'));

    for f = 1:length(excelFiles)
        fileName = upper(excelFiles(f).name);
        expMatch = regexp(fileName, '(EXP\d+)_([A-Z]+)', 'tokens');
        if isempty(expMatch), continue; end
        
        expID = expMatch{1}{1}; 
        strategy = expMatch{1}{2}; 
        if ~ismember(strategy, {'SELFISH', 'GLOBAL', 'HOMEO'}), continue; end

        T = readtable(fullfile(fullPath, excelFiles(f).name));
        if ~all(ismember({'final_Wb_M1', 'final_Wb_M2'}, T.Properties.VariableNames)), continue; end

        avgSelfish = mean([T.final_Wb_M1 == 0, T.final_Wb_M2 == 0], 'all') * 100;

        % Store result
        dataStruct(end+1).exp = expID;
        dataStruct(end).strategy = strategy;
        dataStruct(end).weightRatio = weightRatio;
        dataStruct(end).selfishPercent = avgSelfish;
        dataStruct(end).folder = folderName;
    end
end

% === Plotting ===
if isempty(dataStruct)
    error('No valid data found.');
end

T = struct2table(dataStruct);
expList = unique(T.exp);
strategies = {'SELFISH', 'GLOBAL', 'HOMEO'};
shapes = {'o', 's', '^'}; 
colors = lines(numel(expList)); 

figure('Name', 'Effect of Weight Ratio on Being Selfish');
tiledlayout(2,2, 'Padding', 'compact');

for e = 1:numel(expList)
    nexttile;
    ylim([0 50]);
    expName = expList{e};
    expData = T(strcmp(T.exp, expName), :);
    hold on;

    for s = 1:length(strategies)
        strat = strategies{s};
        shape = shapes{s};
        stratData = expData(strcmp(expData.strategy, strat), :);
        scatter(stratData.weightRatio, stratData.selfishPercent, ...
            80, shape, 'filled', ...
            'DisplayName', strat);
    end

    title(['Experiment ', expName(end)]);
    xlabel('(Wa + Wb) / Wc');
    ylabel('% Time Being Selfish');
    legend('Location','best');
    grid on;
end

sgtitle('Being Selfish vs Weight Ratio Across Strategies and Experiments');
%}

%{
% ===================================================================
%% ======= PLOT 4/4: Refined - 2D Ratio + Color Gradient by Strategy
% ===================================================================

% === Setup ===
baseFolder = 'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons';
weightFolders = dir(baseFolder);
weightFolders = weightFolders([weightFolders.isdir]); 
weightFolders = weightFolders(~ismember({weightFolders.name}, {'.', '..'}));

dataStruct = struct('exp', {}, 'strategy', {}, 'WaWb', {}, 'WcWb', {}, 'selfishPercent', {}, 'folder', {});

for w = 1:length(weightFolders)
    folderName = weightFolders(w).name;
    weightPattern = regexp(folderName, 'Wa(\d+)Wb(\d+)Wc(\d+)', 'tokens');
    if isempty(weightPattern), continue; end

    weights = str2double(weightPattern{1});
    Wa = weights(1); Wb = weights(2); Wc = weights(3);
    if Wb == 0, continue; end  % avoid division by zero
    WaWb = Wa / Wb;
    WcWb = Wc / Wb;

    fullPath = fullfile(baseFolder, folderName);
    excelFiles = dir(fullfile(fullPath, '*.xlsx'));

    for f = 1:length(excelFiles)
        fileName = upper(excelFiles(f).name);
        expMatch = regexp(fileName, '(EXP\d+)_([A-Z]+)', 'tokens');
        if isempty(expMatch), continue; end
        
        expID = expMatch{1}{1}; 
        strategy = expMatch{1}{2}; 
        if ~ismember(strategy, {'SELFISH', 'GLOBAL', 'HOMEO'}), continue; end

        T = readtable(fullfile(fullPath, excelFiles(f).name));
        if ~all(ismember({'final_Wb_M1', 'final_Wb_M2'}, T.Properties.VariableNames)), continue; end

        avgSelfish = mean([T.final_Wb_M1 == 0, T.final_Wb_M2 == 0], 'all') * 100;

        % Store result
        dataStruct(end+1).exp = expID;
        dataStruct(end).strategy = strategy;
        dataStruct(end).WaWb = WaWb;
        dataStruct(end).WcWb = WcWb;
        dataStruct(end).selfishPercent = avgSelfish;
        dataStruct(end).folder = folderName;
    end
end

% === Plotting ===
if isempty(dataStruct)
    error('No valid data found.');
end

T = struct2table(dataStruct);
expList = unique(T.exp);
strategies = {'SELFISH', 'GLOBAL', 'HOMEO'};
shapes = {'o', 's', '^'};  % marker shapes
cmap = parula(256);  % or try hot, jet, viridis (if installed)
colorLimits = [0 100];  % selfish % range

figure('Name', 'Selfishness Landscape: Wa/Wb vs Wc/Wb');

tiledlayout(2,2, 'Padding', 'compact');

for e = 1:numel(expList)
    nexttile;
    expName = expList{e};
    expData = T(strcmp(T.exp, expName), :);
    hold on;

    for s = 1:length(strategies)
        strat = strategies{s};
        shape = shapes{s};
        stratData = expData(strcmp(expData.strategy, strat), :);

        % Normalize color by selfish % (0–100 → 1–256 index)
        colorIdx = round(rescale(stratData.selfishPercent, 1, 256));
        scatterColors = cmap(colorIdx, :);

        scatter(stratData.WaWb, stratData.WcWb, ...
            100, scatterColors, shape, ...
            'filled', 'DisplayName', strat, 'MarkerEdgeColor','k');
    end

    title(['Experiment ', expName(end)]);
    xlabel('Wa / Wb');
    ylabel('Wc / Wb');
    axis equal;
    grid on;
    legend('Location','best');
end

% Colorbar for selfish %
cb = colorbar('Position',[0.93 0.11 0.015 0.815]);
cb.Label.String = '% Time Being Selfish';
cb.Limits = colorLimits;
colormap(cmap);

sgtitle('Refined View: Selfishness by Strategy and Weight Ratios');

%}
