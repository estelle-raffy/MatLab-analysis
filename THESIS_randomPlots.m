close all;
clear all;

% Data
methods = {'GLOBAL ONLY', 'LOCAL', 'NEIGHBOUR', 'GLOBAL', 'HOMEO', 'SELFISH'};
successes = [4, 3, 5, 5, 4, 2];
total_runs = [12, 36, 108, 108, 108, 108];
percentages = (successes ./ total_runs) * 100;

% Create bar plot
figure;
bar(percentages, 'FaceColor', [0.2 0.6 0.8]);
set(gca, 'XTickLabel', methods, 'XTick', 1:numel(methods), 'XTickLabelRotation', 45);
ylabel('Success Rate (%)');
title('Success Rate per Control Strategy');
grid on;

% Annotate bars with exact percentages (1 digit after the decimal)
for i = 1:length(percentages)
    text(i, percentages(i) + 0.5, sprintf('%.1f%%', percentages(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end

%% Check actuation of M3 compared to Global error vs time 
% Load data from Excel file
filename = 'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons_NOGRAV\Wa0.33Wb0.33Wc0.33\EXP1_GLOBAL_ONLY.xlsx';  % Replace with your actual file name
data = readtable(filename);

% Extract relevant columns
time = data.current_time;
error = data.global_error;
actuation = 10 * data.M3_actuation_final;  % Scale actuation by 10

% Plot
figure;
plot(time, error, 'm-', 'LineWidth', 1.5); hold on;                    % Dashed line for error
plot(time, actuation, 'b-', 'LineWidth', 1.5);                         % Solid line for scaled actuation
yline(42, 'r:', 'LineWidth', 1.5);  
xlabel('Time (s)');
ylabel('Value');
title('Global Error vs M3 Actuation (Scaled x10)');
legend('Global Error', 'M3 Actuation x10');
grid on;


%% DATA ANALYSIS FOR GLOBAL ONLY %%
% Define file paths
filePaths = {
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons_NOGRAV\Wa0.33Wb0.33Wc0.33\EXP1_GLOBAL_ONLY.xlsx',
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons_NOGRAV\Wa0.33Wb0.33Wc0.33\EXP2_GLOBAL_ONLY.xlsx',
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons_NOGRAV\Wa0.33Wb0.33Wc0.33\EXP3_GLOBAL_ONLY.xlsx',
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons_NOGRAV\Wa0.33Wb0.33Wc0.33\target_air_500\EXP3_GLOBAL_ONLY.xlsx',
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons_NOGRAV\Wa0.33Wb0.33Wc1\EXP1_GLOBAL_ONLY.xlsx',
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons_NOGRAV\Additional Evidence\EXP2_GLOBAL_ONLY_Wc1_150s.xlsx',
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons_NOGRAV\Additional Evidence\EXP3_GLOBAL_ONLY_Wc1_90s.xlsx',
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons_NOGRAV\Wa0.33Wb0.33Wc1\target_air_500\EXP3_GLOBAL_ONLY.xlsx',
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons_NOGRAV\Wa0.33Wb0.33Wc3\EXP1_GLOBAL_ONLY.xlsx',
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons_NOGRAV\Wa0.33Wb0.33Wc3\EXP2_GLOBAL_ONLY.xlsx',
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons_NOGRAV\Wa0.33Wb0.33Wc3\EXP3_GLOBAL_ONLY.xlsx',
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons_NOGRAV\Wa0.33Wb0.33Wc3\target_air_500\EXP3_GLOBAL_ONLY.xlsx',
};

% Define experiment colors
expColors = containers.Map({'1','2','3','4'}, ...
    {[0 0.45 0.74], [0.85 0.33 0.1], [0.93 0.69 0.13], [0.49 0.18 0.56]});

% Define Wc marker styles
wcMarkers = containers.Map({'0.33','1','3'}, {'o','^','s'});

% Create figure
figure;
set(gcf, 'Position', [100, 100, 1000, 800]);

% Containers to track plotted legend entries per subplot
plottedLegendEntries = cell(4,1);
for idx = 1:4
    plottedLegendEntries{idx} = containers.Map;
end

for i = 1:length(filePaths)
    data = readtable(filePaths{i});
    current_time = data.current_time;
    global_error = data.global_error;

    % Extract Wc_final (assumed numeric)
    if ismember('Wc_final', data.Properties.VariableNames)
        wc_val = data.Wc_final(1);
    else
        wc_val = NaN;
    end
    wc_key = num2str(wc_val);

    % Assign marker shape
    if isKey(wcMarkers, wc_key)
        marker = wcMarkers(wc_key);
    else
        marker = '.';  % fallback marker
    end

    % Get file name (without extension)
    [~, name, ~] = fileparts(filePaths{i});

    % Determine subplot index and experiment number:
    if contains(filePaths{i}, 'target_air_500')
        subplotIdx = 4;  % EXP4 group for target_air_500
        expNum = '4';
    else
        expMatch = regexp(name, 'EXP(\d)', 'tokens');
        if ~isempty(expMatch)
            expNum = expMatch{1}{1};
        else
            expNum = '1'; % default
        end
        subplotIdx = str2double(expNum);
    end

    % Plot in the appropriate subplot
    subplot(4,1,subplotIdx);
    hold on;

    % Plot line with color and markers spaced over data
    h = plot(current_time, global_error, ...
        'Color', expColors(expNum), ...
        'LineWidth', 1.5, ...
        'Marker', marker, ...
        'MarkerIndices', 1:25:length(current_time)); 

    % Add legend entry for marker type only once per subplot
    legendKey = ['Wc = ' wc_key];
    if ~isKey(plottedLegendEntries{subplotIdx}, legendKey)
        plottedLegendEntries{subplotIdx}(legendKey) = h;
    end

    % Fix axis limits
    xlim([0 150]);
    ylim([0 250]);

    % Add red dotted horizontal line at y=42
    yline(42, 'r--', 'LineWidth', 1.5);


    % Format subplot
    title(['EXP' num2str(subplotIdx)], 'FontWeight', 'bold');
    xlabel('Time');
    ylabel('Global Error');
    grid on;
end

% Add one legend per subplot on the right
for subplotIdx = 1:4
    subplot(4,1,subplotIdx);
    keysList = keys(plottedLegendEntries{subplotIdx});
    legendHandles = values(plottedLegendEntries{subplotIdx});
    legendLabels = keysList;

    % Place legend outside to the right
    legend([legendHandles{:}], legendLabels, ...
        'Location', 'eastoutside', ...
        'Box', 'off');
end

%% DATA ANALYSIS FOR LOCAL Wa/Wc Combinations
filePaths = {
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons_NOGRAV\Wa1Wb1Wc1\EXP1_LOCAL.xlsx',
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons_NOGRAV\Wa0.33Wb1Wc3\EXP1_LOCAL.xlsx',
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons_NOGRAV\Wa1Wb1Wc3\EXP1_LOCAL.xlsx',
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons_NOGRAV\Wa3Wb1Wc0.33\EXP1_LOCAL.xlsx',
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons_NOGRAV\Wa1Wb1Wc1\EXP3_LOCAL.xlsx',
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons_NOGRAV\Wa0.33Wb1Wc3\EXP3_LOCAL.xlsx',
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons_NOGRAV\Wa1Wb1Wc3\EXP3_LOCAL.xlsx',
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons_NOGRAV\Wa3Wb1Wc0.33\EXP3_LOCAL.xlsx',
};

% === DEFINE COMBINATIONS TO TRACK ===
targetCombos = {
    struct('Wa', 1,    'Wc', 1,    'label', 'Success (Wa=1, Wc=1)',    'marker', 'o', 'color', [0.2 0.7 0.2]);
    struct('Wa', 0.33, 'Wc', 3,    'label', 'Low Wa (0.33, Wc=3)',     'marker', 's', 'color', [0 0.45 0.74]);
    struct('Wa', 1,    'Wc', 3,    'label', 'High Wc (Wa=1, Wc=3)',    'marker', '^', 'color', [0.85 0.33 0.1]);
    struct('Wa', 3,    'Wc', 0.33, 'label', 'High Wa (3, Wc=0.33)',    'marker', 'd', 'color', [0.49 0.18 0.56]);
};

% === FIGURE SETUP ===
figure('Name', 'LOCAL Strategy - Global Error Analysis');
set(gcf, 'Position', [100, 100, 1000, 900]);

legendTracker = cell(4,1); % Track legend per subplot
for idx = 1:4
    legendTracker{idx} = containers.Map;
end

% === PARSE AND PLOT EACH FILE ===
for i = 1:length(filePaths)
    [~, name, ~] = fileparts(filePaths{i});
    data = readtable(filePaths{i});
    
    if ~all(ismember({'Wa', 'Wc_final', 'current_time', 'global_error'}, data.Properties.VariableNames))
        warning('Skipping %s: Missing expected columns.', name);
        continue;
    end

    Wa = round(data.Wa(1), 2);
    Wc = round(data.Wc_final(1), 2);
    t = data.current_time;
    e = data.global_error;

    % Determine which target combo this belongs to
    comboIdx = -1;
    for k = 1:length(targetCombos)
        if abs(Wa - targetCombos{k}.Wa) < 0.01 && abs(Wc - targetCombos{k}.Wc) < 0.01
            comboIdx = k;
            break;
        end
    end
    if comboIdx == -1
        continue; % not one of the 4 selected cases
    end

    % === Determine EXP# (subplot index) ===
    if contains(name, 'target_air_500', 'IgnoreCase', true)
        subplotIdx = 4;
    else
        expMatch = regexp(name, 'EXP(\d)', 'tokens');
        if ~isempty(expMatch)
            subplotIdx = str2double(expMatch{1}{1});
        else
            subplotIdx = 1; % default fallback
        end
    end

    % === PLOT ===
    subplot(4,1,subplotIdx);
    hold on;

    plot(t, e, ...
        'Color', targetCombos{comboIdx}.color, ...
        'LineWidth', 1.5, ...
        'Marker', targetCombos{comboIdx}.marker, ...
        'MarkerIndices', 1:25:length(t));

    title(['EXP' num2str(subplotIdx)], 'FontWeight', 'bold');
    xlabel('Time (s)');
    ylabel('Global Error');
    grid on;
    xlim([0, 60]);
    ylim([0, 110]);

    % Add reference line
    yline(42, 'r--', 'LineWidth', 1.2);

    % Track legend (only add one per combo)
    label = targetCombos{comboIdx}.label;
    if ~isKey(legendTracker{subplotIdx}, label)
        legendTracker{subplotIdx}(label) = ...
            plot(NaN, NaN, ...
            'Color', targetCombos{comboIdx}.color, ...
            'Marker', targetCombos{comboIdx}.marker, ...
            'LineStyle', '-', ...
            'LineWidth', 1.5);
    end
end

% === ADD LEGENDS TO EACH SUBPLOT ===
for subplotIdx = 1:4
    subplot(4,1,subplotIdx);
    legKeys = keys(legendTracker{subplotIdx});
    legVals = values(legendTracker{subplotIdx});
    if ~isempty(legKeys)
        legend([legVals{:}], legKeys, 'Location', 'eastoutside', 'Box', 'off');
    end
end


% === ADDITIONAL FULL EXP1 PLOT ===
figure('Name','Full Plot - EXP1','Position', [200, 200, 1000, 400]);
hold on;

legendEntries = containers.Map;  % <-- ✅ Initialize here

for i = 1:length(filePaths)
    if contains(filePaths{i}, 'EXP1') && ~contains(filePaths{i}, 'target_air_500')
        data = readtable(filePaths{i});
        current_time = data.current_time;
        global_error = data.global_error;

        if ismember('Wc_final', data.Properties.VariableNames)
            wc_val = data.Wc_final(1);
        else
            wc_val = NaN;
        end
        wc_key = num2str(wc_val);

        if isKey(wcMarkers, wc_key)
            marker = wcMarkers(wc_key);
        else
            marker = '.';
        end

        [~, name, ~] = fileparts(filePaths{i});
        expMatch = regexp(name, 'EXP(\d)', 'tokens');
        if ~isempty(expMatch)
            expNum = expMatch{1}{1};
        else
            expNum = '1';
        end

        h = plot(current_time, global_error, ...
            'Color', expColors(expNum), ...
            'LineWidth', 1.5, ...
            'Marker', marker, ...
            'MarkerIndices', 1:25:length(current_time));

        % ✅ Now this works:
        legendEntries(wc_key) = h;
    end
end


%% Actuation singal Comparison: M1, M2 and M3 (Wc) for the 4 Wa/Wc combinations in the previous graph

filePaths = {
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons_NOGRAV\Wa1Wb1Wc1\EXP1_LOCAL.xlsx',
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons_NOGRAV\Wa0.33Wb1Wc3\EXP1_LOCAL.xlsx',
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons_NOGRAV\Wa1Wb1Wc3\EXP1_LOCAL.xlsx',
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons_NOGRAV\Wa3Wb1Wc0.33\EXP1_LOCAL.xlsx',
};

% === Define Wa/Wc combinations ===
targetCombos = {
    struct('Wa', 1,    'Wc', 1,    'label', 'Wa=1, Wc=1',    'color', [0.2 0.7 0.2]);
    struct('Wa', 0.33, 'Wc', 3,    'label', 'Wa=0.33, Wc=3', 'color', [0 0.45 0.74]);
    struct('Wa', 1,    'Wc', 3,    'label', 'Wa=1, Wc=3',    'color', [0.85 0.33 0.1]);
    struct('Wa', 3,    'Wc', 0.33, 'label', 'Wa=3, Wc=0.33', 'color', [0.49 0.18 0.56]);
};

% === Figure Setup ===
figure('Name', 'Actuation Comparison (M1, M2, M3)');
set(gcf, 'Position', [100, 100, 1000, 900]);

for i = 1:length(filePaths)
    data = readtable(filePaths{i});
    [~, name, ~] = fileparts(filePaths{i});
    
    if ~all(ismember({'Wa', 'Wc_final', 'current_time', ...
                     'M1_actuation_final', 'M2_actuation_final', 'M3_actuation_final'}, ...
                     data.Properties.VariableNames))
        warning('Skipping %s: Missing required columns.', name);
        continue;
    end
    
    Wa = round(data.Wa(1), 2);
    Wc = round(data.Wc_final(1), 2);
    t = data.current_time;
    m1 = data.M1_actuation_final;
    m2 = data.M2_actuation_final;
    m3 = data.M3_actuation_final;

    % Identify which subplot this Wa/Wc combo maps to
    comboIdx = -1;
    for k = 1:length(targetCombos)
        if abs(Wa - targetCombos{k}.Wa) < 0.01 && abs(Wc - targetCombos{k}.Wc) < 0.01
            comboIdx = k;
            break;
        end
    end
    if comboIdx == -1
        warning('Combination Wa=%.2f, Wc=%.2f not in target list.', Wa, Wc);
        continue;
    end

    % === PLOT in Subplot ===
    subplot(4,1,comboIdx);
    hold on;
    
    % Plot M1 (Wa)
    plot(t, m1, '-', 'LineWidth', 1.5, ...
        'DisplayName', 'M1 Actuation');
    
    % Plot M2
    plot(t, m2, '-', 'LineWidth', 1.5, ...
        'DisplayName', 'M2 Actuation');

    % Plot M3 (Wc)
    plot(t, m3, '-', 'LineWidth', 1.5, ...
        'DisplayName', 'M3 Actuation');
    
    title(['Actuation for ' targetCombos{comboIdx}.label], 'FontWeight', 'bold');
    xlabel('Time (s)');
    ylabel('Actuation Value');
    grid on;
    lgd = legend('Location', 'eastoutside', 'Box', 'off');
    xlim([0, max(t)]);
end


% ==============================================================
%% ======= SAME BUT WA / WB / WC FOR DECENTRALISED CONTROL
% ==============================================================

filePaths = {
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons_NOGRAV\Wa1Wb1Wc3\EXP2_SELFISH.xlsx',
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons_NOGRAV\Wa3Wb3Wc3\EXP2_SELFISH.xlsx',
};

% Define combinations with Wa, Wb, and Wc
targetCombos = {
    struct('Wa', 1, 'Wb', 1, 'Wc', 3, 'label', 'Wa=1, Wb=1, Wc=3', 'color', [0 0.45 0.74]);
    struct('Wa', 3, 'Wb', 3, 'Wc', 3, 'label', 'Wa=3, Wb=3, Wc=3', 'color', [0.85 0.33 0.1]);
};

% === Figure Setup ===
figure('Name', 'Actuation Comparison for Varying Wa/Wb (Wc=3)');
set(gcf, 'Position', [100, 100, 900, 700]);

for i = 1:length(filePaths)
    data = readtable(filePaths{i});
    [~, name, ~] = fileparts(filePaths{i});

    if ~all(ismember({'Wa', 'final_Wb_M1', 'Wc_final', 'current_time', ...
                      'M1_actuation_final', 'M2_actuation_final', 'M3_actuation_final'}, ...
                      data.Properties.VariableNames))
        warning('Skipping %s: Missing required columns.', name);
        continue;
    end

    Wa = round(data.Wa(1), 2);
    Wb = round(data.final_Wb_M1(1), 2);
    Wc = round(data.Wc_final(1), 2);
    t = data.current_time;
    m1 = data.M1_actuation_final;
    m2 = data.M2_actuation_final;
    m3 = data.M3_actuation_final;

    % Match combination
    comboIdx = -1;
    for k = 1:length(targetCombos)
        if abs(Wa - targetCombos{k}.Wa) < 0.01 && ...
           abs(Wb - targetCombos{k}.Wb) < 0.01 && ...
           abs(Wc - targetCombos{k}.Wc) < 0.01
            comboIdx = k;
            break;
        end
    end
    if comboIdx == -1
        warning('Combination Wa=%.2f, Wb=%.2f, Wc=%.2f not matched.', Wa, Wb, Wc);
        continue;
    end

    % === PLOT in Subplot ===
    subplot(2,1,comboIdx);
    hold on;

    plot(t, m1, '-', 'LineWidth', 1.5, 'DisplayName', 'M1 Actuation');
    plot(t, m2, '-', 'LineWidth', 1.5, 'DisplayName', 'M2 Actuation');
    plot(t, m3, '-', 'LineWidth', 1.5, 'DisplayName', 'M3 Actuation');

    title(['Actuation for ' targetCombos{comboIdx}.label], 'FontWeight', 'bold');
    xlabel('Time (s)');
    ylabel('Actuation Value');
    grid on;
    legend('Location', 'eastoutside', 'Box', 'off');
    xlim([0, max(t)]);
end

% ==============================================================
%% ======= COMPARING LOCAL/FINAL WB FOR TASK 1 BETWEEN N/S/G 
% ==============================================================

% === FILE PATHS (EDIT THESE) ===
filePaths = {
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons_NOGRAV\Wa1Wb1Wc0.33\EXP1_NEIGHBOUR.xlsx',
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons_NOGRAV\Wa1Wb1Wc0.33\EXP1_SELFISH.xlsx',
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons_NOGRAV\Wa1Wb1Wc0.33\EXP1_GLOBAL.xlsx',
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons_NOGRAV\Wa1Wb1Wc0.33\EXP1_HOMEO.xlsx',
};

% === STRATEGY IDENTIFIERS ===
strategies = {'NEIGHBOUR', 'SELFISH', 'GLOBAL', 'HOMEO'};
strategyColors = containers.Map(...
    {'M1_local', 'M1_final', 'M2_local', 'M2_final'}, ...
    {[0 0.2 0.8], [0.6 0.8 1], [1 0 0], [1 0.7 0.4]});  % now: M1_local = dark blue, M1_final = light blue, M2_local = red, M2_final = orange

% === SETUP FIGURE ===
figure('Name', 'Wb Evolution per Strategy');
set(gcf, 'Position', [200, 100, 1200, 900]);

% === LOOP OVER STRATEGIES ===
for s = 1:length(strategies)
    strategy = strategies{s};

    % === FIND FILES FOR STRATEGY ===
    strategyFiles = filePaths(contains(filePaths, strategy, 'IgnoreCase', true));

    if isempty(strategyFiles)
        warning('No files found for strategy: %s', strategy);
        continue;
    end

    % === SETUP SUBPLOT ===
    subplot(length(strategies),1,s);  % Changed from 2,2,s to 4,1,s to stack vertically
    hold on;
    title(strategy, 'FontWeight', 'bold', 'FontSize', 12);
    xlabel('Time (s)');
    ylabel('W_b value');
    grid on;

    % === PLOT EACH FILE IN STRATEGY GROUP ===
    for f = 1:length(strategyFiles)
        data = readtable(strategyFiles{f});

        if ~all(ismember({'current_time', 'local_Wb_M1', 'local_Wb_M2', 'final_Wb_M1', 'final_Wb_M2'}, ...
                         data.Properties.VariableNames))
            warning('Missing columns in: %s', strategyFiles{f});
            continue;
        end

        t = data.current_time;
        lM1 = data.local_Wb_M1;
        lM2 = data.local_Wb_M2;
        fM1 = data.final_Wb_M1;
        fM2 = data.final_Wb_M2;

        % Plot lines
        plot(t, lM1, '-', 'Color', strategyColors('M1_local'), 'LineWidth', 2, 'DisplayName', 'local W_b M1');
        plot(t, fM1, '-', 'Color', strategyColors('M1_final'), 'LineWidth', 2, 'DisplayName', 'final W_b M1');
        plot(t, lM2, '-', 'Color', strategyColors('M2_local'), 'LineWidth', 2, 'DisplayName', 'local W_b M2');
        plot(t, fM2, '-', 'Color', strategyColors('M2_final'), 'LineWidth', 2, 'DisplayName', 'final W_b M2');
    end

    legend({'local W_b M1', 'final W_b M1', 'local W_b M2', 'final W_b M2'}, ...
           'Location', 'eastoutside', 'Box', 'off');
    ylim([0 1.05]); % adjust if needed
    hold off;
end




% ==============================================================
% ======= Energy Tracking per File with Zero-Actuation Highlight
% ==============================================================

%% === Define file paths manually ===
filePaths = {
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons_NOGRAV\Wa1Wb1Wc1\EXP1_GLOBAL_ONLY.xlsx',
    'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_comparisons_NOGRAV\Wa1Wb1Wc1\EXP1_LOCAL.xlsx',
};

% === PARAMETERS ===
stiffness = 100;  % Spring stiffness constant
strategyGroups = {'LOCAL', 'NEIGHBOUR', 'SELFISH', 'GLOBAL', 'HOMEO', 'GLOBAL_ONLY'};
numFiles = length(filePaths);

% === FIGURE SETUP ===
figure('Name','Energy per File');

for i = 1:numFiles
    filePath = filePaths{i};
    [~, fileName, ext] = fileparts(filePath);  % Get just the name for title
    fileName = [fileName, ext];  % Add extension back for display
    T = readtable(filePath);

    % Extract area and time
    A1 = T.area_M1;
    A2 = T.area_M2;
    A3 = T.area_M3;
    t = T.current_time;

    % Reference areas at t=0
    A1_ref = A1(1);
    A2_ref = A2(1);
    A3_ref = A3(1);

    % Compute energy over time
    E = 0.5 * stiffness * ((A1 - A1_ref).^2 + (A2 - A2_ref).^2 + (A3 - A3_ref).^2);

    % Identify actuation columns
    m1 = T.M1_actuation_final;
    m2 = T.M2_actuation_final;
    m3 = T.M3_actuation_final;
    isZeroAll = (m1 == 0) & (m2 == 0) & (m3 == 0);

    % Find all ≥5s periods of zero actuation
    d = diff([0; isZeroAll; 0]);
    starts = find(d == 1);
    ends = find(d == -1) - 1;
    validStarts = [];
    validEnds = [];

    for z = 1:length(starts)
        duration = t(ends(z)) - t(starts(z));
        if duration >= 5
            validStarts(end+1) = starts(z);
            validEnds(end+1) = ends(z);
        end
    end

    % Determine strategy from filename
    fname = upper(fileName);
    groupLabel = 'UNKNOWN';
    for k = 1:length(strategyGroups)
        if contains(fname, strategyGroups{k}, 'IgnoreCase', true)
            groupLabel = strategyGroups{k};
            break;
        end
    end

    % === PLOT ===
    subplot(numFiles, 1, i);
    plot(t, E, 'b', 'LineWidth', 1.5);
    xlim([0, 60]);
    ylim([0, 6e7]);

    hold on;

    % Highlight all valid zero-actuation periods
    for j = 1:length(validStarts)
        tStart = t(validStarts(j));
        tEnd = t(validEnds(j));
        idxRange = (t >= tStart & t <= tEnd);
        area(t(idxRange), E(idxRange), 'FaceColor', [1, 0.7, 0.7], 'EdgeColor', 'none');
    end

    % Extract strategy (groupLabel) - already done
    % Extract EXP number from filename using regexp
    expMatch = regexp(fname, 'EXP[_\s]?0*(\d+)', 'tokens', 'ignorecase');
    if ~isempty(expMatch)
        expNum = expMatch{1}{1};  % Get the number part as a string
    else
        expNum = 'Unknown';
    end
    
    % Update title to "STRATEGY - EXP#"
    title(sprintf('%s — EXP%s', groupLabel, expNum), 'Interpreter', 'none');

    ylabel('Energy');
    if i == numFiles
        xlabel('Time (s)');
    end
    grid on;
    xlim([0, 60]);
end
