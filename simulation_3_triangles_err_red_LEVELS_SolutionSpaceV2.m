clear; clc;

% Folder with your individual Excel files
dataFolder = 'C:\\Users\\om21104\\OneDrive - University of Bristol\\Desktop\\Project SC\\Results\\3_Modules_NAIVE\\LEVELS\\THESIS_Tests\\Solution_space_excel_files_homeo';
files = dir(fullfile(dataFolder, '*.xlsx'));

% Control strategy identifiers (used in file names)
strategyGroups = {'LOCAL', 'SELFISH', 'GLOBAL_ONLY', 'GLOBAL'};
groupNames = {'Local', 'Neigh/Selfish', 'Global Only', 'Global'};
numGroups = numel(strategyGroups);
strategyColors = lines(numGroups); % bright and distinct colors

% ===============================
%% ======= PLOT 1/3: Final Shapes (Springs + Bodies)
% ===============================

figure;
hold on;
title('1/3: Final Shapes of Successful Solutions');
axis equal;
xlabel('X'); ylabel('Y');
grid on;

% Define spring connections
spring_pairs = [1 3; 3 2; 1 2; 2 5; 3 5; 5 4; 2 4];

theta = linspace(0, 2*pi, 20);
radius_body = 0.05;

for i = 1:length(files)
    T = readtable(fullfile(files(i).folder, files(i).name), 'Sheet', 'Tabelle1');
    if T.success_log(end) == 0, continue; end

    % Final body positions
    bx = [T.body1_pos_x(end), T.body2_pos_x(end), T.body3_pos_x(end), ...
          T.body4_pos_x(end), T.body5_pos_x(end)];
    by = [T.body1_pos_y(end), T.body2_pos_y(end), T.body3_pos_y(end), ...
          T.body4_pos_y(end), T.body5_pos_y(end)];

    % Anchor to first body
    bx = bx - bx(1);
    by = by - by(1);

    % Color by strategy
    fname = upper(files(i).name);
    groupIdx = [];

    if contains(fname, 'GLOBAL_ONLY')
        groupIdx = find(strcmp(strategyGroups, 'GLOBAL_ONLY'));
    elseif contains(fname, 'SELFISH')
        groupIdx = find(strcmp(strategyGroups, 'SELFISH'));
    elseif contains(fname, 'GLOBAL')
        groupIdx = find(strcmp(strategyGroups, 'GLOBAL'));
    elseif contains(fname, 'LOCAL')
        groupIdx = find(strcmp(strategyGroups, 'LOCAL'));
    end
    if isempty(groupIdx), continue; end
    color = strategyColors(groupIdx, :);

    % Plot springs
    for j = 1:size(spring_pairs, 1)
        i1 = spring_pairs(j, 1);
        i2 = spring_pairs(j, 2);
        plot([bx(i1), bx(i2)], [by(i1), by(i2)], '-', 'Color', color, 'LineWidth', 2);
    end

    % Plot bodies
    for j = 1:5
        fill(bx(j) + radius_body*cos(theta), by(j) + radius_body*sin(theta), ...
             color, 'FaceAlpha', 0.4, 'EdgeColor', 'k');
    end
end

legend(groupNames);

% ===============================
%% ======= PLOT 2/3: Area Dynamics vs Time (4 subplots)
% ===============================
%{
figure('Name','2/3: Module Dynamics Over Time');
sgtitle('2/3: Evolution of Module Areas for Each Strategy');

for g = 1:numGroups
    subplot(2,2,g);
    hold on; grid on;
    title(groupNames{g});
    xlabel('Time'); ylabel('Area Size');
    
    cmap = lines(3); % colors for 3 modules
    
    for i = 1:length(files)
        T = readtable(fullfile(files(i).folder, files(i).name), 'Sheet', 'Tabelle1');
        if T.success_log(end) == 0, continue; end

        fname = upper(files(i).name);
        % Fix group detection (SELFISH before GLOBAL)
        if contains(fname, 'GLOBAL_ONLY')
            groupIdx = find(strcmp(strategyGroups, 'GLOBAL_ONLY'));
        elseif contains(fname, 'SELFISH')
            groupIdx = find(strcmp(strategyGroups, 'SELFISH'));
        elseif contains(fname, 'GLOBAL')
            groupIdx = find(strcmp(strategyGroups, 'GLOBAL'));
        elseif contains(fname, 'LOCAL')
            groupIdx = find(strcmp(strategyGroups, 'LOCAL'));
        else
            continue;
        end

        if groupIdx ~= g, continue; end

        % Area over time
        A1 = T.area_M1;
        A2 = T.area_M2;
        A3 = T.area_M3;

        % Actuation final (absolute)
        a1 = abs(T.M1_actuation_final(end));
        a2 = abs(T.M2_actuation_final(end));
        a3 = abs(T.M3_actuation_final(end));
        [~, dominantIdx] = max([a1, a2, a3]);
        color = cmap(dominantIdx, :);

        % Normalize time
        t = 1:length(A1);
        plot(t, A1, '-', 'Color', color, 'LineWidth', 0.8, 'HandleVisibility','off');
        plot(t, A2, '-', 'Color', color, 'LineWidth', 0.8, 'HandleVisibility','off');
        plot(t, A3, '-', 'Color', color, 'LineWidth', 0.8, 'HandleVisibility','off');
    end
end
%}

% ===============================
%% ======= PLOT 3/3: Error (bars) + Local Frustration (line)
% ===============================

avgError = zeros(1, numGroups);
stdError = zeros(1, numGroups);
avgFrustration = zeros(1, numGroups);
count = zeros(1, numGroups);

for i = 1:length(files)
    T = readtable(fullfile(files(i).folder, files(i).name), 'Sheet', 'Tabelle1');
    if T.success_log(end) == 0, continue; end

    fname = upper(files(i).name);

    if contains(fname, 'GLOBAL_ONLY')
        groupIdx = find(strcmp(strategyGroups, 'GLOBAL_ONLY'));
    elseif contains(fname, 'SELFISH')
        groupIdx = find(strcmp(strategyGroups, 'SELFISH'));
    elseif contains(fname, 'GLOBAL')
        groupIdx = find(strcmp(strategyGroups, 'GLOBAL'));
    elseif contains(fname, 'LOCAL')
        groupIdx = find(strcmp(strategyGroups, 'LOCAL'));
    else
        continue;
    end

    globalErr = mean(T.global_error);
    frustration = mean((T.M1_local_frustration + T.M2_local_frustration)/2);

    avgError(groupIdx) = avgError(groupIdx) + globalErr;
    avgFrustration(groupIdx) = avgFrustration(groupIdx) + frustration;
    stdError(groupIdx) = stdError(groupIdx) + std(T.global_error);
    count(groupIdx) = count(groupIdx) + 1;
end

avgError = avgError ./ count;
avgFrustration = avgFrustration ./ count;
stdError = stdError ./ count;

figure;
yyaxis left;
b = bar(1:numGroups, avgError, 0.5);
b.FaceColor = [0.6 0.85 1];
hold on;
errorbar(1:numGroups, avgError, stdError, 'k.', 'LineWidth', 1.2);
ylabel('Average Global Error');

yyaxis right;
plot(1:numGroups, avgFrustration, '-o', 'Color', [1 0.3 0.3], 'LineWidth', 2);
ylabel('Average Local Frustration');

set(gca, 'XTick', 1:numGroups, 'XTickLabel', groupNames);
xlabel('Control Strategy');
title('3/3: Error vs Frustration per Control Strategy');
grid on;

% ===============================
%% ======= PLOT 4/4: Integrated Information (Φ)
% ===============================

avgPhi = zeros(1, numGroups);
stdPhi = zeros(1, numGroups);
count = zeros(1, numGroups);

for i = 1:length(files)
    T = readtable(fullfile(files(i).folder, files(i).name), 'Sheet', 'Tabelle1');
    if T.success_log(end) == 0, continue; end

    fname = upper(files(i).name);

    if contains(fname, 'GLOBAL_ONLY')
        groupIdx = find(strcmp(strategyGroups, 'GLOBAL_ONLY'));
    elseif contains(fname, 'SELFISH')
        groupIdx = find(strcmp(strategyGroups, 'SELFISH'));
    elseif contains(fname, 'GLOBAL')
        groupIdx = find(strcmp(strategyGroups, 'GLOBAL'));
    elseif contains(fname, 'LOCAL')
        groupIdx = find(strcmp(strategyGroups, 'LOCAL'));
    else
        continue;
    end

    M1Error = T.M1_local_err;
    M2Error = T.M2_local_err;

    phi1 = temporal_mutual_information(M1Error);
    phi2 = temporal_mutual_information(M2Error);

    phi = (phi1 + phi2) / 2;

    avgPhi(groupIdx) = avgPhi(groupIdx) + phi;
    stdPhi(groupIdx) = stdPhi(groupIdx) + std([phi1, phi2]);
    count(groupIdx) = count(groupIdx) + 1;
end

avgPhi = avgPhi ./ count;
stdPhi = stdPhi ./ count;

figure;
barHandle = bar(1:numGroups, avgPhi, 0.5);
hold on;
errorbar(1:numGroups, avgPhi, stdPhi, 'k.', 'LineWidth', 1.2);
set(barHandle, 'FaceColor', [0.7 0.9 0.4]);
set(gca, 'XTick', 1:numGroups, 'XTickLabel', groupNames);
ylabel('Average Integrated Information (Φ)');
xlabel('Control Strategy');
title('4/4: Integrated Information vs Control Strategy');
grid on;

% ===============================
%% Helper Function: Temporal Mutual Information
% ===============================
function phi = temporal_mutual_information(signal)
    signal = signal(~isnan(signal));
    if numel(signal) < 2
        phi = NaN;
        return;
    end
    signal = (signal - mean(signal)) / std(signal);
    x_t = signal(1:end-1);
    x_t1 = signal(2:end);
    r = corr(x_t, x_t1);

    if abs(r) < 1
        phi = -0.5 * log(1 - r^2);
    else
        phi = 0;
    end
end
