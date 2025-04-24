%% === Setup ===
clear; clc; close all;
% File path
filePath = 'C:\Users\om21104\OneDrive - University of Bristol\Desktop\Project SC\Results\Hardware\Baseline vs 3 starting positions\Tests.xlsx';
% Sheet names
sheets = {'Baseline', 'Stimulus 3900', 'Stimulus 3500'};
labels = {'Baseline', 'Stimulus A (3900)', 'Stimulus B (3500)'};
colors = [0.2 0.6 1; 1 0.4 0.4; 0.6 0.8 0.3]; % Blue, Red, Green
% Initialize results
mean_actuation = zeros(3,3);  % Rows = conditions, Columns = motors
std_actuation = zeros(3,3);
% Loop through each sheet
for i = 1:3
    data = readtable(filePath, 'Sheet', sheets{i});
    data.Motor = categorical(data.Motor);
    motor_ids = categories(data.Motor);
    for j = 1:length(motor_ids)
        motor_data = data(data.Motor == motor_ids{j}, :);
        mean_actuation(i,j) = mean(motor_data.Final_Actuation);
        std_actuation(i,j) = std(motor_data.Final_Actuation);
    end
end
%% === Plot: Error Bar Chart ===
figure('Name','Final Actuation per Motor with Error Bars','Color','w');
hold on;
x = 1:3; % Motors
bar_width = 0.25;
offsets = [-bar_width, 0, bar_width];
% Pre-store bar positions to adjust limits
all_bar_positions = [];
% Draw bars first (with edge color to make baseline visible)
for cond = 1:3
    bar_pos = x + offsets(cond);
    b = bar(bar_pos, mean_actuation(cond,:), bar_width, ...
        'FaceColor', colors(cond,:), ...
        'EdgeColor', 'k', 'LineWidth', 0.8, ...
        'DisplayName', labels{cond});
    all_bar_positions = [all_bar_positions, bar_pos];
end
% Now draw error bars *on top* of the bars
for cond = 1:3
    bar_pos = x + offsets(cond);
    errorbar(bar_pos, mean_actuation(cond,:), std_actuation(cond,:), ...
        'k.', 'LineWidth', 1.5, 'CapSize', 10);
end
% Improve axis
xticks(x);
xticklabels({'Motor 0', 'Motor 1', 'Motor 2'});
ylabel('Mean Final Actuation');
title('Final Actuation per Motor with Error Bars');
legend('Location','northeastoutside');
ylim padded;  % Automatically expand Y-axis slightly
grid on;
box on;