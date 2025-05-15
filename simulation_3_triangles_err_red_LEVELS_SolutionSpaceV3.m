clear; clc;

% Folder with your Excel files
dataFolder = 'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons\Wa3Wb3Wc3';
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
%% ======= PLOT 1.a/2: Final Shapes Grouped by Strategy (Updated, Label Underneath)
% ===============================================
figure('Name','1.a/2: Final Shapes Grouped by Strategy');
hold on;
axis equal;
xlabel('X');
ylabel('Y');
title('1.a/2: Diversity of Solutions by Strategy');
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

    bx = [T.body1_pos_x(end), T.body2_pos_x(end), T.body3_pos_x(end), ...
          T.body4_pos_x(end), T.body5_pos_x(end)];
    by = [T.body1_pos_y(end), T.body2_pos_y(end), T.body3_pos_y(end), ...
          T.body4_pos_y(end), T.body5_pos_y(end)];

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
%% ======= PLOT 1.b/2: All Final Shapes Grouped by Strategy (Grey if failed)
% ===============================================
figure('Name','1.b/2: All Final Shapes by Strategy (Success = Colour, Failed = Grey)');
hold on;
axis equal;
xlabel('X');
ylabel('Y');
title('1.b/2: All Final Solutions by Strategy (Success = Colour, Failed = Grey)');
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

    % Final body positions
    bx = [T.body1_pos_x(end), T.body2_pos_x(end), T.body3_pos_x(end), ...
          T.body4_pos_x(end), T.body5_pos_x(end)];
    by = [T.body1_pos_y(end), T.body2_pos_y(end), T.body3_pos_y(end), ...
          T.body4_pos_y(end), T.body5_pos_y(end)];

    % Normalize
    bx = bx - bx(1);
    bx = bx - min(bx);
    if max(bx) > 0
        bx = bx / max(bx) * 2;
    end
    by = by - by(1);
    by = by - min(by);
    if max(by) > 0
        by = by / max(by) * 1;
    end

    % Apply offset
    bx = bx + dx;
    by = by + dy;

    % Color: grey if failed
    if T.success_log(end) == 0
        color = [0.7, 0.7, 0.7];
    end

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



% ===============================================
%% ======= PLOT 2/2: Stability 
% ===============================================

% === Setup ===
figure('Name','Actuation Silence Periods (All Zeros)');
hold on;
title('Periods Where M1, M2, M3 = 0 (Grouped by Strategy)');
ylabel('Start Time (s)');
xlabel('Strategy Group');
grid on;

groupSilentStart = cell(1, numGroups);
groupSilentDuration = cell(1, numGroups);
groupSilentStart_expNames = cell(1, numGroups);
groupSilentStart_expColors = cell(1, numGroups);

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

    % Experiment name detection & color assignment (match plot 2/5 logic)
    expMatch = regexp(fname, '(EXP[_\s]?0*\d+)', 'tokens', 'ignorecase');
    if ~isempty(expMatch)
        expName = upper(strrep(expMatch{1}{1}, '_', ''));
    else
        expName = sprintf('File%d', i);
    end

    if ~ismember(expName, uniqueExpNames)
        uniqueExpNames{end+1} = expName;
    end
    expIdx = find(strcmp(uniqueExpNames, expName));
    expColor = expColors(expIdx, :);

    % Check for simultaneous zero actuation
    m1 = T.M1_actuation_final;
    m2 = T.M2_actuation_final;
    m3 = T.M3_actuation_final;
    t = T.current_time;

    isZeroAll = (m1 == 0) & (m2 == 0) & (m3 == 0);

    % Find longest contiguous period of all-zero actuation
    d = diff([0; isZeroAll; 0]);
    starts = find(d == 1);
    ends = find(d == -1) - 1;

    if isempty(starts)
        continue;
    end

    durations = t(ends) - t(starts);
    [maxDur, idxMax] = max(durations);
    onsetTime = t(starts(idxMax));

    fprintf('File: %s, onsetTime: %.3f seconds\n', files(i).name, onsetTime);

    % Store for plotting
    groupSilentStart{groupIdx}(end+1) = onsetTime;
    groupSilentDuration{groupIdx}(end+1) = maxDur;
    groupSilentStart_expNames{groupIdx}{end+1} = expName;
    groupSilentStart_expColors{groupIdx}(end+1, :) = expColor;
end

% === Plotting ===
groupSpacing = 5;
x_ticks = [];
x_labels = {};

anyFound = false;
plottedExps = {};
legendEntries = {};
legendColors = [];

for k = 1:numGroups
    starts = groupSilentStart{k};
    durations = groupSilentDuration{k};

    if exist('groupSilentStart_expNames', 'var') && length(groupSilentStart_expNames) >= k
        expNamesGroup = groupSilentStart_expNames{k};
    else
        expNamesGroup = {};
    end
    if exist('groupSilentStart_expColors', 'var') && length(groupSilentStart_expColors) >= k
        expColorsGroup = groupSilentStart_expColors{k};
    else
        expColorsGroup = [];
    end

    n = numel(starts);
    if n == 0
        % No data, but still add x tick for group label
        x_ticks(end+1) = (k - 1) * groupSpacing + 1.5;
        x_labels{end+1} = groupNames{k};
        continue;
    end

    anyFound = true;

    for j = 1:n
        x = (k - 1) * groupSpacing + j;

        % Use default color if expColorsGroup empty or too short
        if isempty(expColorsGroup) || size(expColorsGroup,1) < j
            colorToUse = strategyColors(k,:);
        else
            colorToUse = expColorsGroup(j,:);
        end

        startTime = starts(j);
        duration = durations(j);
        endTime = startTime + duration;
        midTime = (startTime + endTime) / 2;
        halfDuration = duration / 2;
        
        % Plot the vertical error bar: line from start to end
        line([x, x], [startTime, endTime], ...
             'Color', colorToUse, ...
             'LineWidth', 1.5);
        
        % Plot the center marker at mid-point (like classic error bar center)
        plot(x, midTime, 'o', ...
             'MarkerFaceColor', colorToUse, ...
             'MarkerEdgeColor', 'k', ...
             'MarkerSize', 8);
        
        % Optional: add caps at top and bottom like classic error bars
        capWidth = 0.2;
        line([x - capWidth, x + capWidth], [startTime, startTime], ...
             'Color', colorToUse, 'LineWidth', 1.2);
        line([x - capWidth, x + capWidth], [endTime, endTime], ...
             'Color', colorToUse, 'LineWidth', 1.2);



        % Add legend entry if not already added
        if ~isempty(expNamesGroup) && j <= length(expNamesGroup)
            expLabel = expNamesGroup{j};
            if ~ismember(expLabel, plottedExps)
                legendEntries{end+1} = expLabel;
                legendColors(end+1, :) = colorToUse;
                plottedExps{end+1} = expLabel;
            end
        end
    end

    % Group label in center under bars
    mid_x = (k - 1) * groupSpacing + n / 2 + 0.5;
    x_ticks(end+1) = mid_x;
    x_labels{end+1} = groupNames{k};
end

if anyFound
    xticks(x_ticks);
    xticklabels(x_labels);
    ylim([0 65]);

    % Add legend bars for experiment colors
    hold on;
    hLegend = gobjects(1, numel(legendEntries));
    for i = 1:numel(legendEntries)
        hLegend(i) = plot(nan, nan, 'o', ...
            'MarkerFaceColor', legendColors(i,:), ...
            'MarkerEdgeColor', 'k', ...
            'MarkerSize', 8);
    end
    legend(hLegend, legendEntries, 'Location', 'eastoutside');
else
    text(0.5, 0.5, 'No simultaneous zero actuation found.', ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'FontSize', 12, 'FontWeight', 'bold', 'Color', [0.5 0 0]);
end

xlim([0, max(x_ticks)+groupSpacing]);
