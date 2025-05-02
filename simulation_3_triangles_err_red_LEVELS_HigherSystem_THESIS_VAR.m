%% Load the Excel File 
filePath = 'C:\Users\om21104\OneDrive - University of Bristol\Documents\Python Scripts\3_triangles_LEVELS_THESIS-TESTS_homeo.xlsx';
rawTable = readtable(filePath, 'Sheet', 'Sheet');

%% Extract Positions for Bodies
time = rawTable.current_time;
y0 = rawTable.body4_pos_x; y1 = rawTable.body4_pos_y;
y4 = rawTable.body1_pos_x; y5 = rawTable.body1_pos_y;
y6 = rawTable.body2_pos_x; y7 = rawTable.body2_pos_y;
y8 = rawTable.body3_pos_x; y9 = rawTable.body3_pos_y;
y10 = rawTable.body5_pos_x; y11 = rawTable.body5_pos_y;
y2 = rawTable.target_x; y3 = rawTable.target_y;

%% Define Circle and Spring Properties
body_diameter = 21;
radius_body = body_diameter / 2;
radius_target = body_diameter / 2;
radius_threshold = 45;
sampling_rate = 10;

%% Create Figure
figure; hold on;
history_bodies = cell(1, 5);
gif_filename = 'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\3_Modules_NAIVE\LEVELS\THESIS_Tests\FUNCTIONALITIES_COMPARISONS\Target_experiments_LOCAL\my_animation.gif';

indices = 1:sampling_rate:length(time);
last_index = indices(end);

%% Loop through time steps
for idx = 1:length(indices)
    i = indices(idx);
    clf; hold on;

    % Update history
    history_bodies{1} = [history_bodies{1}; y4(i), y5(i)];
    history_bodies{2} = [history_bodies{2}; y6(i), y7(i)];
    history_bodies{3} = [history_bodies{3}; y8(i), y9(i)];
    history_bodies{4} = [history_bodies{4}; y0(i), y1(i)];
    history_bodies{5} = [history_bodies{5}; y10(i), y11(i)];

    % Target
    theta = linspace(0, 2*pi, 100);
    x_target_circle = y2(i) + radius_target * cos(theta);
    y_target_circle = y3(i) + radius_target * sin(theta);
    fill(x_target_circle, y_target_circle, 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    
    % Threshold
    rectangle('Position', [y2(i) - radius_threshold, y3(i) - radius_threshold, ...
              2 * radius_threshold, 2 * radius_threshold], ...
              'Curvature', [1, 1], 'EdgeColor', 'y', 'LineWidth', 2);

    % Bodies
    bodies_x = [y4(i), y6(i), y8(i), y0(i), y10(i)];
    bodies_y = [y5(i), y7(i), y9(i), y1(i), y11(i)];
    
    for j = 1:length(bodies_x)
        color = 'b'; if j == 4, color = 'g'; end
        fill(bodies_x(j) + radius_body * cos(theta), ...
             bodies_y(j) + radius_body * sin(theta), ...
             color, 'FaceAlpha', 0.6, 'EdgeColor', 'k');
    end

    % Springs
    spring_pairs = [1 3; 3 2; 1 2; 2 5; 3 5; 5 4; 2 4];
    for j = 1:size(spring_pairs, 1)
        plot([bodies_x(spring_pairs(j,1)), bodies_x(spring_pairs(j,2))], ...
             [bodies_y(spring_pairs(j,1)), bodies_y(spring_pairs(j,2))], ...
             'k-', 'LineWidth', 2);
    end

    % Past trajectories
    for j = 1:length(history_bodies)
        if size(history_bodies{j}, 1) > 1
            past_x = history_bodies{j}(:, 1);
            past_y = history_bodies{j}(:, 2);
            color = [0, 0, 1]; if j == 4, color = [0, 1, 0]; end
            plot(past_x, past_y, '-', 'Color', color, 'LineWidth', 1.5, 'HandleVisibility', 'off');
        end
    end

    % Past positions fading
    for j = 1:length(history_bodies)
        for k = 1:size(history_bodies{j}, 1)
            alpha = k / size(history_bodies{j}, 1);
            x_circle = history_bodies{j}(k, 1) + radius_body * cos(theta);
            y_circle = history_bodies{j}(k, 2) + radius_body * sin(theta);
            color = 'b'; if j == 4, color = 'g'; end
            fill(x_circle, y_circle, color, 'FaceAlpha', 1 - alpha, 'EdgeColor', 'k');
        end
    end

    % Add Obstacle (replacing the rectangle function)
    % Define the size of the obstacle (same as your Python code)
    obstacle_size = 40;
    half_size = obstacle_size / 2;

    % Create a square-like obstacle with 'fill' instead of 'rectangle'
    x_obstacle = [obstacle_x(i) - half_size, obstacle_x(i) + half_size, obstacle_x(i) + half_size, obstacle_x(i) - half_size];
    y_obstacle = [obstacle_y(i) - half_size, obstacle_y(i) - half_size, obstacle_y(i) + half_size, obstacle_y(i) + half_size];

    % Fill the obstacle with semi-transparency
    fill(x_obstacle, y_obstacle, 'm', 'FaceAlpha', 0.4, 'EdgeColor', 'm');


    % Plot Formatting
    xlim([-300 300]); ylim([0 300]); axis equal;
    title('EXP4 LOCAL Wa1 Wb0 Wc0 WnM1 1 WnM2 -1 M1M2 M2M3 Demand400 Target500-550 Thresh45 T60', ...
          sprintf('System Configuration at t = %.1f sec', time(i)));

    % Legend
    h1 = plot(NaN, NaN, 'bo', 'DisplayName', 'Bodies', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
    h2 = plot(NaN, NaN, 'ko-', 'DisplayName', 'Springs', 'LineWidth', 2);
    h3 = plot(NaN, NaN, 'ro', 'DisplayName', 'Target', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    h4 = plot(NaN, NaN, 'yo', 'DisplayName', 'Global Threshold', 'MarkerSize', 15, 'MarkerEdgeColor', 'y', 'LineWidth', 1.5);
    legend([h1, h2, h3, h4], 'Location', 'northeast');

    % === Show success/fail message on final frame ===
    if i == last_index
        success_log_at_end = rawTable.success_log(end);
        if success_log_at_end == 1
            success_status = 'SUCCEEDS';
            msg_color = 'green';
        else
            success_status = 'FAILS';
            msg_color = 'red';
        end
        text(-200, 200, success_status, 'Color', msg_color, ...
             'FontSize', 20, 'FontWeight', 'bold', 'BackgroundColor', 'white');
    end

    % Save frame to GIF
    frame = getframe(gcf);
    img = frame2im(frame);
    [imind, cm] = rgb2ind(img, 256);
    
    if idx == 1
        imwrite(imind, cm, gif_filename, 'gif', 'LoopCount', Inf, 'DelayTime', 0.4);
    else
        imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.4);
    end

    pause(0.1); % Optional visualization pause
end

%% Variance Analysis (unchanged from your script)
area_M1 = rawTable.area_M1;
area_M2 = rawTable.area_M2;
area_M3 = rawTable.area_M3;
success_log = rawTable.success_log;
time_vector = rawTable.current_time;

target_time = 60;
[~, idx_60] = min(abs(time_vector - target_time));

area_M1_segment = area_M1(1:idx_60);
area_M2_segment = area_M2(1:idx_60);
area_M3_segment = area_M3(1:idx_60);

var_M1 = var(area_M1_segment);
var_M2 = var(area_M2_segment);
var_M3 = var(area_M3_segment);

avg_variance = mean([var_M1, var_M2, var_M3]);

var_diff_M1_M2 = abs(var_M1 - var_M2);
var_diff_M1_M3 = abs(var_M1 - var_M3);
var_diff_M2_M3 = abs(var_M2 - var_M3);

std_variance = std([var_M1, var_M2, var_M3]);

if success_log(idx_60) == 1
    success_status = '✅ SUCCESS at t = 60';
else
    success_status = '❌ FAILURE at t = 60';
end

fprintf('--- Variance Summary (t = 0 to 60) ---\n');
fprintf('%s\n', success_status);
fprintf('Module 1 Variance: %.4f\n', var_M1);
fprintf('Module 2 Variance: %.4f\n', var_M2);
fprintf('Module 3 Variance: %.4f\n', var_M3);
fprintf('Average Variance: %.4f\n', avg_variance);
fprintf('Variance Differences: |M1-M2| = %.4f, |M1-M3| = %.4f, |M2-M3| = %.4f\n', var_diff_M1_M2, var_diff_M1_M3, var_diff_M2_M3);
fprintf('Standard Deviation of Variances: %.4f\n', std_variance);

figure;
subplot(4, 1, 1);
bar([var_M1, var_M2, var_M3]);
xlabel('Module'); ylabel('Variance of Area'); title('Variance of Each Module');

subplot(4, 1, 2);
bar(avg_variance, 'FaceColor', [0.5, 0, 0.5]);
xlabel('Metric'); ylabel('Average Variance of Areas'); title('Average Variance of System');

subplot(4, 1, 3);
bar([var_diff_M1_M2, var_diff_M1_M3, var_diff_M2_M3]);
xlabel('Variance differences M1-M2, M1-M3, M2-M3');
ylabel('Variance Differences');
title('Variance Differences Between Modules');

subplot(4, 1, 4);
bar(std_variance, 'FaceColor', [0, 0.5, 0]);
xlabel('Metric'); ylabel('Standard Deviation of Variances');
title('Standard Deviation of Variances');
