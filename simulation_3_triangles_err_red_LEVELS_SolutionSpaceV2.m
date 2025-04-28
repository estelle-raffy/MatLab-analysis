clear; clc;

% Folder with your Excel files
dataFolder = 'C:\\Users\\om21104\\OneDrive - University of Bristol\\Desktop\\Project SC\\Results\\3_Modules_NAIVE\\LEVELS\\THESIS_Tests\\Solution_space_excel_files_homeo';
files = dir(fullfile(dataFolder, '*.xlsx'));

%% DEBUG SELFISH 
disp('Checking filenames for SELFISH...');
for i = 1:length(files)
    if contains(upper(files(i).name), 'SELFISH')
        disp(['Found SELFISH file: ', files(i).name]);
    end
end

% Define strategy groups
strategyGroups = {'LOCAL', 'SELFISH', 'GLOBAL_ONLY', 'GLOBAL'};
groupNames = {'Local', 'Neigh/Selfish', 'Global Only', 'Global'};
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
%% ======= PLOT 1/4: Final Shapes Grouped by Strategy
% ===============================================

figure;
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

for i = 1:length(files)
    T = readtable(fullfile(files(i).folder, files(i).name), 'Sheet', 'Tabelle1');
    if T.success_log(end) == 0, continue; end

    % Final body positions
    bx = [T.body1_pos_x(end), T.body2_pos_x(end), T.body3_pos_x(end), ...
          T.body4_pos_x(end), T.body5_pos_x(end)];
    by = [T.body1_pos_y(end), T.body2_pos_y(end), T.body3_pos_y(end), ...
          T.body4_pos_y(end), T.body5_pos_y(end)];

    % Normalization
    bx = bx - bx(1);
    by = by - by(1);

    bx = bx - min(bx);
    bx = bx / max(bx) * 2; % normalized width
    by = by - min(by);
    by = by / max(by) * 1; % normalized height

    % Determine group
    fname = upper(files(i).name);
    groupIdx = [];
    for k = 1:numGroups
        if contains(fname, strategyGroups{k})
            groupIdx = k;
            break;
        end
    end
    if isempty(groupIdx), warning('Unknown group in file: %s', files(i).name); continue; end
    color = strategyColors(groupIdx, :);

    % Compute offset
    idx = groupSolutionCount(groupIdx);
    dx = (groupIdx - 1) * spacing_x;
    dy = -mod(idx, max_per_column) * spacing_y;
    dx = dx + floor(idx / max_per_column) * (spacing_x/2); % if too many, shift right

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

    % Plot balls
    for j = 1:5
        fill(bx(j) + radius_body*cos(theta), by(j) + radius_body*sin(theta), ...
            color, 'FaceAlpha', 0.6, 'EdgeColor', 'k', 'LineWidth', 0.5);
    end
end

% Add group names as labels
for k = 1:numGroups
    xpos = (k-1) * spacing_x + 1;
    ypos = 2; % Adjust vertically if needed
    text(xpos, ypos, groupNames{k}, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 12);
end

% ===============================================
%% ======= PLOT 2/4: Area Dynamics Over Time
% ===============================================

figure('Name','2/4: Module Dynamics Over Time');
sgtitle('2/4: Evolution of Module Areas for Each Strategy');

for g = 1:numGroups
    subplot(2,2,g);
    hold on; grid on;
    title(groupNames{g});
    xlabel('Time');
    ylabel('Area Size');

    cmap = lines(2); % colors: blue for M1, orange for M2

    hM1 = []; % handle for M1 dominant curves
    hM2 = []; % handle for M2 dominant curves

    for i = 1:length(files)
        T = readtable(fullfile(files(i).folder, files(i).name), 'Sheet', 'Tabelle1');
        
        % DEBUG: Check SELFISH success
        fname = upper(files(i).name);
        if contains(fname, 'SELFISH')
            disp(['SELFISH file: ', files(i).name, ' | success_log(end) = ', num2str(T.success_log(end))]);
        end

        if T.success_log(end) == 0, continue; end

        fname = upper(files(i).name);
        match = contains(fname, strategyGroups{g});
        if ~match, continue; end

        A1 = T.area_M1;
        A2 = T.area_M2;

        a1 = abs(T.M1_actuation_final(end));
        a2 = abs(T.M2_actuation_final(end));
        [~, dominantIdx] = max([a1, a2]);
        color = cmap(dominantIdx, :);

        t = 1:length(A1);
        
        % Plot curves and store the handle for legend
        if dominantIdx == 1
            h = plot(t, A1, '-', 'Color', color, 'LineWidth', 1);
            plot(t, A2, '-', 'Color', color, 'LineWidth', 1);
            if isempty(hM1)
                hM1 = h; % Save first M1 dominant curve for legend
            end
        else
            h = plot(t, A1, '-', 'Color', color, 'LineWidth', 1);
            plot(t, A2, '-', 'Color', color, 'LineWidth', 1);
            if isempty(hM2)
                hM2 = h; % Save first M2 dominant curve for legend
            end
        end
    end

    % === NEW: Proper Legend with Correct Color Handles ===
    if ~isempty(hM1) && ~isempty(hM2)
        legend([hM1, hM2], {'M1 Dominant', 'M2 Dominant'}, 'Location', 'best');
    elseif ~isempty(hM1)
        legend(hM1, {'M1 Dominant'}, 'Location', 'best');
    elseif ~isempty(hM2)
        legend(hM2, {'M2 Dominant'}, 'Location', 'best');
    end
end

% ===================================================
% === FUTURE: If you want to include M3 as dominant ===
% Uncomment this block and comment the old one above
%
% cmap = lines(3); % colors for M1, M2, M3
% a1 = abs(T.M1_actuation_final(end));
% a2 = abs(T.M2_actuation_final(end));
% a3 = abs(T.M3_actuation_final(end)); % NEW
% [~, dominantIdx] = max([a1, a2, a3]);
% color = cmap(dominantIdx, :);
%
% Handle separately for hM1, hM2, hM3 and adjust legend
%
% ===================================================


% ===============================================
%% ======= PLOT 3/4: Global Error vs Local Frustration
% ===============================================

avgError = zeros(1, numGroups);
stdError = zeros(1, numGroups);
avgFrustration = zeros(1, numGroups);
count = zeros(1, numGroups);

for i = 1:length(files)
    T = readtable(fullfile(files(i).folder, files(i).name), 'Sheet', 'Tabelle1');
    if T.success_log(end) == 0, continue; end

    fname = upper(files(i).name);

    groupIdx = [];
    for k = 1:numGroups
        if contains(fname, strategyGroups{k})
            groupIdx = k;
            break;
        end
    end
    if isempty(groupIdx), continue; end

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
title('3/4: Error vs Frustration');
grid on;

%% ======= PLOT 4/4: Integrated Information (Φ)
% ===============================================

avgPhi = zeros(1, numGroups);  % Initialize arrays for average Φ
stdPhi = zeros(1, numGroups);  % Standard deviation for error bars
count = zeros(1, numGroups);   % Count to average Φ over the number of files

% Loop through each file to compute Φ for each strategy
for i = 1:length(files)
    T = readtable(fullfile(files(i).folder, files(i).name), 'Sheet', 'Tabelle1');
    if T.success_log(end) == 0, continue; end  % Skip files with unsuccessful runs

    fname = upper(files(i).name);  % Extract file name to determine group

    % Determine the control strategy group based on the file name
    groupIdx = [];
    for k = 1:numGroups
        if contains(fname, strategyGroups{k})
            groupIdx = k;
            break;
        end
    end
    if isempty(groupIdx), continue; end  % Skip if no strategy group is found

    % Extract relevant columns from the table for the current file
    M1Error = T.M1_new_local_err;  % Local error for M1
    M2Error = T.M2_new_local_err;  % Local error for M2
    M3Error = T.global_error;    % Global error for M3 (for GLOBAL and GLOBAL_ONLY)

    % Determine if M1 and M2 are behaving selfishly
    beingSelfishM1 = T.beingSelfishM1;
    beingSelfishM2 = T.beingSelfishM2;

    % Compute mutual information for LOCAL strategy (M1 and M2 errors)
    % In LOCAL, M1 and M2 are mechanically coupled, so we calculate MI between them
    phi1 = temporal_mutual_information(M1Error);
    phi2 = temporal_mutual_information(M2Error);
    phi_local = (phi1 + phi2) / 2;

    % For SELFISH strategy, consider errors and neighbor differences
    M1_neigh_diff = T.M1_neigh_diff;  % Neighboring difference for M1
    M2_neigh_diff = T.M2_neigh_diff;  % Neighboring difference for M2

    % Calculate mutual information for SELFISH strategy
    if any(beingSelfishM1 == 1 | beingSelfishM2 == 1)  % If either is selfish, use LOCAL behavior
        phi_selfish = (phi1 + phi2) / 2;  % Revert to LOCAL
    else  % Collective behavior
        phi_selfish = temporal_mutual_information([M1Error, M2Error, M1_neigh_diff, M2_neigh_diff]);
    end

    % For GLOBAL strategy, include M3 error
    phi_global = temporal_mutual_information([M1Error, M2Error, M3Error]);

    % For GLOBAL_ONLY strategy, only use M3 error
    phi_global_only = temporal_mutual_information(M3Error);

    % Store Φ values based on the control strategy
    if groupIdx == 1  % LOCAL strategy
        avgPhi(groupIdx) = avgPhi(groupIdx) + phi_local;
        stdPhi(groupIdx) = stdPhi(groupIdx) + std([phi1, phi2]);
    elseif groupIdx == 2  % SELFISH strategy
        avgPhi(groupIdx) = avgPhi(groupIdx) + phi_selfish;
        stdPhi(groupIdx) = stdPhi(groupIdx) + std([phi1, phi2]);
    elseif groupIdx == 3  % GLOBAL strategy
        avgPhi(groupIdx) = avgPhi(groupIdx) + phi_global;
        stdPhi(groupIdx) = stdPhi(groupIdx) + std([phi1, phi2, phi_global]);
    elseif groupIdx == 4  % GLOBAL_ONLY strategy
        avgPhi(groupIdx) = avgPhi(groupIdx) + phi_global_only;
        stdPhi(groupIdx) = stdPhi(groupIdx) + std(phi_global_only);
    end
    count(groupIdx) = count(groupIdx) + 1;  % Increment the count for averaging
end

% Average the results
avgPhi = avgPhi ./ count;
stdPhi = stdPhi ./ count;

% Plot the integrated information (Φ) for each strategy
figure;
barHandle = bar(1:numGroups, avgPhi, 0.5);
hold on;
errorbar(1:numGroups, avgPhi, stdPhi, 'k.', 'LineWidth', 1.2);
set(barHandle, 'FaceColor', [0.7 0.9 0.4]);  % Set bar color
set(gca, 'XTick', 1:numGroups, 'XTickLabel', groupNames);  % Set strategy names
ylabel('Average Integrated Information (Φ)');
xlabel('Control Strategy');
title('4/4: Integrated Information (Φ)');
grid on;


% ===============================================
%% Helper Function
% ===============================================
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
