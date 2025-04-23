% === SETUP ===
clear; clc;

% Folder with your individual Excel files
dataFolder = 'C:/Users/om21104/OneDrive - University of Bristol/Desktop/Project SC/Results/3_Modules_NAIVE/LEVELS/THESIS_Tests/Solution_space_excel_files_homeo';
files = dir(fullfile(dataFolder, '*.xlsx'));

% Shape-only feature extraction (5 bodies = 10 coords)
shapeFeatureNames = {};
for b = 1:5
    shapeFeatureNames{end+1} = sprintf('body%d_pos_x', b);
    shapeFeatureNames{end+1} = sprintf('body%d_pos_y', b);
end

% === LOAD AND EXTRACT DATA ===
shapeFeatures = [];
finalShapes = {};
dominantWeights = [];
frustrationLevels = [];
pointLabels = {};
strategyColors = {};
expLabels = {};

for i = 1:length(files)
    filePath = fullfile(files(i).folder, files(i).name);
    T = readtable(filePath, 'Sheet', 'Tabelle1');
    last = height(T);

    % Only consider successful experiments
    if T.success_log(last) ~= 1
        continue;
    end

    % Extract shape features
    f = [];
    for j = 1:length(shapeFeatureNames)
        f(end+1) = T.(shapeFeatureNames{j})(last);
    end
    shapeFeatures(end+1,:) = f;
    finalShapes{end+1} = f;

    % Frustration
    frust = (T.M1_local_frustration(last) + T.M2_local_frustration(last)) / 2;
    frustrationLevels(end+1) = frust;

    % Dominant weight
    if ismember('dominant_weight', T.Properties.VariableNames)
        dominantWeight = T.dominant_weight(last);
    else
        dominantWeight = 1;
    end
    dominantWeights(end+1) = dominantWeight;

    % Strategy and ID
    nameParts = split(files(i).name, {'_', '.'});
    taskID = upper(nameParts{1});
    strategy = upper(nameParts{2});
    pointLabels{end+1} = sprintf('%s\\n%.2f', taskID, dominantWeight);
    expLabels{end+1} = taskID;

    % Strategy color
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
    strategyColors{end+1} = col;
end

% PCA on shape features
[coeff, score, ~, ~, explained] = pca(shapeFeatures);
X = score(:, 1:2);

% Normalize frustration for coloring
fMin = min(frustrationLevels);
fMax = max(frustrationLevels);
normFrust = (frustrationLevels - fMin) / (fMax - fMin);
cmap = [0 1 0; 1 0 0]; % green to red
interpColors = interp1([0 1], cmap, normFrust');

% === PLOT ===
figure; hold on;
title('Shape-based PCA of Final Successful System Configurations');
xlabel(sprintf('PC1 (%.1f%%)', explained(1)));
ylabel(sprintf('PC2 (%.1f%%)', explained(2)));

theta = linspace(0, 2*pi, 30);
radius_body = 0.03;

% Plot final shapes
for i = 1:size(X,1)
    shape = finalShapes{i};
    col = strategyColors{i};
    fillColor = interpColors(i,:);

    % Body positions relative to PCA location
    offsetX = X(i,1);
    offsetY = X(i,2);

    for b = 1:5
        bx = shape(2*b-1);
        by = shape(2*b);
        px = offsetX + 0.03 * bx;
        py = offsetY + 0.03 * by;

        fill(px + radius_body*cos(theta), py + radius_body*sin(theta), ...
             fillColor, 'EdgeColor', col, 'FaceAlpha', 0.8, 'LineWidth', 1.2);
    end

    text(offsetX + 0.02, offsetY, pointLabels{i}, 'FontSize', 8, 'Color', [0.2 0.2 0.2]);
end

colormap(cmap); colorbar;
ylabel(colorbar, 'Avg. Local Frustration');
grid on; axis equal;
