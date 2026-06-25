clear; clc; close all;

%% Fixed parameters

% Fixed number of vertices
n = 4000;

% Values of c covering the theoretical transition region
% More sampling points are used to obtain a smoother curve.
c_values = linspace(-1, 2, 30);
num_c = length(c_values);

% Number of Monte Carlo samples
num_samples = 800;

% Preallocate memory for the simulation results
simu_result = zeros(1, num_c);

%% Theoretical upper and lower bounds

% Theoretical lower bound
theoretical_value_low = ...
    exp(-2 .* exp(-c_values)) .* ...
    (1 ...
    + 2 .* exp(-c_values) ...
    + exp(-2 .* c_values));

% Theoretical upper bound
theoretical_value_up = ...
    1 ...
    - exp(-2 .* exp(-c_values)) .* exp(-2 .* c_values);


%% Main simulation loop: vary c
for c_idx = 1:num_c
    c = c_values(c_idx);

    fprintf('Processing c = %.6f\n', c);

    % Edge probability
    % In MATLAB, log(n) denotes the natural logarithm.
    p = (log(n) + c) / n;

    % Ensure that p is a valid probability
    p = max(0, min(p, 1));

    %% Monte Carlo simulation
    count = 0;

    for sample = 1:num_samples

        % Generate a random adjacency matrix.
        % Every possible entry, including diagonal entries,
        % is independently generated with probability p.
        A = double(rand(n) < p);

        % Convert the adjacency matrix to the system matrix convention
        A_sys = A';

        % Compute the generic rank using a random numerical realization
        generic_rank = rank(A_sys .* rand(n));

        % Construct the cost matrix for minimum-cost matching
        A_cost = 5001 * ones(n);

        % Existing edges are assigned zero cost
        A_cost(A_sys == 1) = 0;

        % A missing diagonal edge is assigned unit cost
        diagonal_indices = logical(eye(n));
        A_cost(diagonal_indices) = 1 - diag(A_sys);

        % Compute a minimum-cost perfect matching
        [~, min_totalcost] = munkres(A_cost);

        % Check the structural diagonalizability condition
        if generic_rank == n - min_totalcost
            count = count + 1;
        end
    end

    % Empirical probability of structural diagonalizability
    simu_result(c_idx) = count / num_samples;

    fprintf(['Completed c = %.2f: p = %.6f, ', ...
             'simulation probability = %.3f\n'], ...
             c, p, simu_result(c_idx));
end

%% Plot the results
h_figure = figure( ...
    'Position', [100, 100, 1400, 800], ...
    'Color', 'w');

hold on;
grid on;

% Plot the theoretical lower bound
plot( ...
    c_values, ...
    theoretical_value_low, ...
    'b-', ...
    'LineWidth', 4, ...
    'DisplayName', 'Theoretical Lower Bound');

% Plot the theoretical upper bound
plot( ...
    c_values, ...
    theoretical_value_up, ...
    'r-', ...
    'LineWidth', 4, ...
    'DisplayName', 'Theoretical Upper Bound');

% Plot the Monte Carlo simulation results
plot( ...
    c_values, ...
    simu_result, ...
    'o-', ...
    'Color', [227, 207, 87] / 255, ...
    'LineWidth', 4, ...
    'MarkerSize', 12, ...
    'MarkerFaceColor', 'y', ...
    'DisplayName', 'Simulation Results');

%% Axis labels
xlabel( ...
    'Parameter c', ...
    'FontSize', 22, ...
    'FontWeight', 'bold');

ylabel( ...
    'Structural Diagonalizability Probability', ...
    'FontSize', 22, ...
    'FontWeight', 'bold');

%% Axis settings
ax = gca;
ax.FontSize = 22;
ax.LineWidth = 1.8;
ax.GridAlpha = 0.4;
ax.LabelFontSizeMultiplier = 1.1;

xlim([min(c_values), max(c_values)]);
ylim([0, 1.05]);

box on;

%% Legend settings
legend( ...
    'Location', 'southeast', ...
    'FontSize', 22);

hold off;

%% Print the complete simulation results
fprintf('\nSimulation summary:\n');
fprintf('Number of vertices: n = %d\n', n);
fprintf('Number of Monte Carlo samples: %d\n\n', num_samples);

fprintf('c-value\t\tEdge probability p\tSimulation probability\n');
fprintf('------------------------------------------------------------\n');

for c_idx = 1:num_c
    c = c_values(c_idx);
    p = (log(n) + c) / n;
    p = max(0, min(p, 1));

    fprintf('%.6f\t%.8f\t\t%.4f\n', ...
        c, p, simu_result(c_idx));
end