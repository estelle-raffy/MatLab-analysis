clear; clc;
% Folder with your individual Excel files
dataFolder = 'C:\\Users\\om21104\\OneDrive - University of Bristol\\Desktop\\Project SC\\Results\\3_Modules_NAIVE\\LEVELS\\THESIS_Tests\\Solution_space_excel_files_homeo';
files = dir(fullfile(dataFolder, '*.xlsx'));
% Shape-only feature extraction (5 bodies = 10 coords)
shapeFeatureNames = {};
for b = 1:5
    shapeFeatureNames{end+1} = sprintf('body%d_pos_x', b);
    shapeFeatureNames{end+1} = sprintf('body%d_pos_y', b);
end
% === LOAD AND EXTRACT DATA ===
shapeFeatures = [];
frustrationLevels = [];
successFlags = [];
dominantWeights = [];
pointLabels = {};
colors = [];
shapes = [];
strategyMap = containers.Map();
expLabels = {};
for i = 1:length(files)
    % Load data
    filePath = fullfile(files(i).folder, files(i).name);
    T = readtable(filePath, 'Sheet', 'Tabelle1');
    last = height(T); % Final time step
    % Extract SHAPE features only
    f = [];
    for j = 1:length(shapeFeatureNames)
        f(end+1) = T.(shapeFeatureNames{j})(last);
    end
    shapeFeatures(end+1,:) = f;
    % Local frustration (avg of M1 + M2)
    frust = (T.M1_local_frustration(last) + T.M2_local_frustration(last)) / 2;
    frustrationLevels(end+1) = frust;
    % Success status
    successFlags(end+1) = T.success_log(last);
    % Dominant weight by majority vote
    weights_all = [T.Wa, T.final_Wb_M1, T.final_Wb_M2, T.Wc_final];
    [~, idxs] = max(weights_all, [], 2); % find dominant at each timestep
    majorityVote = mode(idxs);
    dominantWeights(end+1) = majorityVote;
    % ID and strategy from filename
    nameParts = split(files(i).name, {'_', '.'});
    taskID = upper(nameParts{1});
    strategy = upper(nameParts{2});
    pointLabels{end+1} = sprintf('%s_%s', taskID, strategy);
    expLabels{end+1} = sprintf('%s (%s)', taskID, ['W' num2str(majorityVote)]); % Label with dominant weight
    % Color (by strategy)
    if contains(strategy, 'LOCAL') && ~contains(strategy, 'GLOBAL')
        col = [1 1 0]; % yellow
    elseif contains(strategy, 'SELFISH')
        col = [1.0 0.6 0.1]; % orange
    elseif contains(strategy, 'GLOBAL_ONLY')
        col = [0.8 0.3 0.8]; % purple
    elseif contains(strategy, 'GLOBAL')
        col = [0.3 0.6 1.0]; % blue
    else
        col = [0.5 0.5 0.5]; % fallback
    end
    colors(end+1,:) = col;
    % Shape (o or x)
    if T.success_log(last) == 1
        shapes{end+1} = 'o';
    else
        shapes{end+1} = 'x';
    end
    % Store for arrow mapping
    strategyMap(pointLabels{end}) = i;
end
% === PCA: Shape Only ===
[coeff, score, ~, ~, explained] = pca(shapeFeatures);
X = score(:,1:2);
% Normalize for visual weights
wMin = min(dominantWeights);
wMax = max(dominantWeights);
normWeights = (dominantWeights - wMin) / (wMax - wMin);
bubbleSizes = 50 + 100 * normWeights;
alphas = 0.2 + 0.8 * normWeights;
% === PLOT: Solution Space ===
figure; hold on;
title('Solution Space (Shape-Based PCA)');
xlabel(sprintf('PC1 (%.1f%%)', explained(1)));
ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
% Frustration gradient
fMin = min(frustrationLevels);
fMax = max(frustrationLevels);
normFrustration = (frustrationLevels - fMin) / (fMax - fMin);
cmap = [0 1 0; 1 0 0]; % green to red
interpColors = interp1([0 1], cmap, normFrustration');
% Plot each point
for i = 1:length(X)
    scatter(X(i,1), X(i,2), bubbleSizes(i), interpColors(i,:), shapes{i}, ...
        'MarkerFaceColor', interpColors(i,:), ...
        'MarkerEdgeColor', colors(i,:), ...
        'MarkerFaceAlpha', alphas(i), ...
        'MarkerEdgeAlpha', alphas(i), ...
        'LineWidth', 1.5);
    text(X(i,1) + 0.02, X(i,2), expLabels{i}, 'FontSize', 8, 'Color', [0.3 0.3 0.3]);
end
% === PCA Convex Hulls per strategy ===
strategies = {'LOCAL', 'SELFISH', 'GLOBAL', 'GLOBAL_ONLY'};
strategyColors = {
    [1 1 0]; [1.0 0.6 0.1]; [0.3 0.6 1.0]; [0.8 0.3 0.8]
};
strategyVolumes = containers.Map;
for s = 1:length(strategies)
    strat = strategies{s};
    stratIdx = find(contains(pointLabels, strat));
    if length(stratIdx) >= 3
        pts = X(stratIdx,:);
        k = convhull(pts(:,1), pts(:,2));
        fill(pts(k,1), pts(k,2), strategyColors{s}, ...
            'FaceAlpha', 0.0, 'EdgeColor', strategyColors{s}, 'LineWidth', 1.2); % FULLY TRANSPARENT FILL
        strategyVolumes(strat) = polyarea(pts(k,1), pts(k,2));
    else
        strategyVolumes(strat) = 0;
    end
end
disp('=== PCA Region Areas per Strategy ===');
disp(strategyVolumes);
% === LEGENDS ===
h_success = scatter(nan, nan, 70, 'o', 'filled', 'MarkerEdgeColor', 'k');
h_fail = scatter(nan, nan, 70, 'x', 'MarkerEdgeColor', 'k');
h_local = plot(nan, nan, 'o', 'MarkerEdgeColor', [1 1 0], 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
h_selfish = plot(nan, nan, 'o', 'MarkerEdgeColor', [1.0 0.6 0.1], 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
h_global = plot(nan, nan, 'o', 'MarkerEdgeColor', [0.3 0.6 1.0], 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
h_global_only = plot(nan, nan, 'o', 'MarkerEdgeColor', [0.8 0.3 0.8], 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
legend([h_success h_fail h_local h_selfish h_global h_global_only], ...
    {'Success', 'Failure', 'Local', 'Selfish', 'Global', 'Global Only'}, ...
    'Location', 'bestoutside');
colormap(cmap); colorbar;
ylabel(colorbar, 'Avg. Local Frustration');
grid on; axis equal;
% Enable data cursor
datacursormode on;
