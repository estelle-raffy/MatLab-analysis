%% Load the Excel File
filePath = 'C:\Users\om21104\OneDrive - University of Bristol\Documents\Python Scripts\ExperimentsResults_Arduino-Python_tests';
rawTable = readtable(filePath, 'Sheet', 'Sheet');

data = rawTable;
data.Motor = categorical(data.Motor);
motor_ids = unique(data.Motor);

% Color settings
colors = [1 0 0; 0 1 0; 0.5 0 0.5]; % Red, Green, Purple for motors
weight_colors = [0 0.4470 0.7410; 0.9290 0.6940 0.1250]; % Blue (selfish), Yellow (cooperative)

%% === HYPOTHESIS 1 ===
figure('Name','Hypothesis 1: Global Control Transitions');

% --- Plot 1.1: Global Frustration vs System Error (categorical color)
subplot(2,1,1)

% Global weight and local weight
g_weight = data.Global_Neighbour_Weight;
l_weight = data.Local_Neighbour_Weight;

% Identify when global weight matches local weight (same) or differs (different)
same_weights = g_weight == l_weight;  % When global weight matches local weight
different_weights = g_weight ~= l_weight;  % When global weight differs from local weight

% Scatter plot for global weight = local weight (blue)
scatter(data.Global_Frustration(same_weights), data.Error_python(same_weights), 40, 'b', 'o');
hold on;

% Scatter plot for global weight ≠ local weight (red)
scatter(data.Global_Frustration(different_weights), data.Error_python(different_weights), 40, 'r', 'filled');

% Labels and title
xlabel('Global Frustration');
ylabel('System Error (Error\_python)');
title('1.1: Global Frustration vs System Error');

% Legends for matching and non-matching weights
legend('Global = Local (same)', 'Global ≠ Local (different)');

% Grid and hold off to finish the plot
grid on; 
hold off;


% --- UPDATED Plot 1.2: Local vs Global Weights Per Motor
subplot(2,1,2)
hold on;
motor_labels = {'Motor 0', 'Motor 1', 'Motor 2'};
for i = 1:length(motor_ids)
    m_label = motor_labels{i};
    motor_data = data(data.Motor == motor_ids(i), :);

    % Local weight = what selfishness function chose
    plot(motor_data.Sample_Number, motor_data.Local_Neighbour_Weight, '-', ...
        'Color', colors(i,:), 'LineWidth', 2, ...
        'DisplayName', [m_label ' Local']);

    % Global weight = what the system allowed
    plot(motor_data.Sample_Number, motor_data.Global_Neighbour_Weight, '--', ...
        'Color', colors(i,:), 'LineWidth', 2, ...
        'DisplayName', [m_label ' Global']);
end

xlabel('Sample Number'); ylabel('Neighbour Weight');
yticks([0 1]); ylim([-0.1 1.1]);
yticklabels({'Selfish (0)', 'Mimic (1)'});
title('1.2: Local vs Global Neighbour Weights per Motor Over Time');
legend('Location','eastoutside');
grid on; hold off;

%% === HYPOTHESIS 2 ===
figure('Name','Hypothesis 2: Local Conflict vs System Performance');

subplot(2,1,1)
samples = unique(data.Sample_Number);
mean_frustration = arrayfun(@(s) mean(data.Local_Frustration(data.Sample_Number == s)), samples);
error_vals = arrayfun(@(s) data.Error_python(find(data.Sample_Number == s, 1)), samples);
scatter(mean_frustration, error_vals, 40, 'b', 'filled');
xlabel('Mean Local Frustration per Timestep');
ylabel('System Error (Error\_python)');
title('2.1: Mean Local Frustration vs System Error');
legend('Samples'); % NEW: Added legend to clarify blue dots
grid on;

% --- Plot 2.2: Final Actuation vs Local Error (Color = Local Weight)
subplot(2,1,2)
lw = data.Local_Neighbour_Weight;
scatter(data.Local_Error(lw==0), data.Final_Actuation(lw==0), 40, weight_colors(1,:), 'filled');
hold on;
scatter(data.Local_Error(lw==1), data.Final_Actuation(lw==1), 40, weight_colors(2,:), 'filled');
xlabel('Local Error'); ylabel('Final Actuation');
title('2.2: Final Actuation vs Local Error (Color = Local Control)');
legend('Local W = 0 (selfish)', 'Local W = 1 (cooperative)');
grid on; hold off;

%% === HYPOTHESIS 3 ===
figure('Name','Hypothesis 3: Behavioral Diversity from Local Control');

% --- Plot 3.1: Histogram of Final Actuation by Local Weight
subplot(2,1,1)
hold on;
edges = -60:10:60;
for i = 1:length(unique(lw))
    val = unique(lw);
    idx = data.Local_Neighbour_Weight == val(i);
    histogram(data.Final_Actuation(idx), edges, ...
        'DisplayName', ['Local W = ' num2str(val(i))], ...
        'FaceAlpha', 0.6, 'FaceColor', weight_colors(i,:));
end
xlabel('Final Actuation Value'); ylabel('Frequency');
title('3.1: Final Actuation Histogram by Local Control');
legend; hold off;

% --- Plot 3.2: Cumulative Actuation per Motor
subplot(2,1,2)
hold on;
for i = 1:length(motor_ids)
    idx = data.Motor == motor_ids(i);
    cumAct = cumsum(data.Final_Actuation(idx));
    plot(data.Sample_Number(idx), cumAct, '-', ...
        'DisplayName', ['Motor ' char(motor_ids(i))], ...
        'Color', colors(i,:));
end
xlabel('Time (Sample Number)'); ylabel('Cumulative Actuation');
title('3.2: Cumulative Actuation per Motor Over Time');
legend; grid on; hold off;
