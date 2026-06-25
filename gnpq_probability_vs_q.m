clear; clc; close all;

%% Global font settings
set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultTextFontName', 'Times New Roman');
set(groot, 'defaultLegendFontName', 'Times New Roman');

%% Fixed parameters

% Fixed number of vertices
n = 4000;

% Fixed parameter c
c_fixed = 0;

% Probability of each non-self-loop directed edge
% In MATLAB, log(n) denotes the natural logarithm.
p = (log(n) + c_fixed) / n;

% Ensure that p is a valid probability
p = max(0, min(p, 1));

% Self-loop probabilities
q_values = 0:0.1:1;
num_q = length(q_values);

% Number of Monte Carlo samples
num_samples = 800;

% Fix the random seed for reproducibility
rng(1);

% Preallocate memory for the simulation results
simu_result = zeros(1, num_q);

%% Theoretical upper and lower bounds for G(n,p,q)

% Theoretical lower bound:
theoretical_value_low = ...
    exp(-2 .* (1 - q_values) .* exp(-c_fixed)) .* ...
    (1 ...
    + 2 .* (1 - q_values) .* exp(-c_fixed) ...
    + (1 - q_values).^2 .* exp(-2 * c_fixed));

% Theoretical upper bound:
theoretical_value_up = ...
    1 ...
    - (1 - q_values).^2 .* exp(-2 * c_fixed) .* ...
      exp(-2 .* (1 - q_values) .* exp(-c_fixed));

%% Main simulation loop: vary q
for q_idx = 1:num_q
    q = q_values(q_idx);

    fprintf('Processing q = %.1f\n', q);
    fprintf('Fixed n = %d, p = %.8f\n', n, p);

    %% Monte Carlo simulation
    count = 0;

    for sample = 1:num_samples

        % Generate all matrix entries initially with probability p.
        % The diagonal entries will subsequently be replaced according
        % to the prescribed self-loop probability q.
        A = double(rand(n) < p);

        % Generate the diagonal entries independently with probability q
        diagonal_support = double(rand(n, 1) < q);
        A(1:n+1:end) = diagonal_support;

        % Convert the adjacency matrix to the system matrix convention
        A_sys = A';

        % Compute the generic rank of the zero-nonzero pattern
        generic_rank = sprank(sparse(A_sys));

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
    simu_result(q_idx) = count / num_samples;

    fprintf(['Completed q = %.1f: ', ...
             'simulation probability = %.4f\n\n'], ...
             q, simu_result(q_idx));
end

%% Plot the results
h_figure = figure( ...
    'Position', [100, 100, 1400, 800], ...
    'Color', 'w');

hold on;
grid on;

% Plot the theoretical lower bound
plot( ...
    q_values, ...
    theoretical_value_low, ...
    'b-', ...
    'LineWidth', 4, ...
    'DisplayName', 'Theoretical Lower Bound');

% Plot the theoretical upper bound
plot( ...
    q_values, ...
    theoretical_value_up, ...
    'r-', ...
    'LineWidth', 4, ...
    'DisplayName', 'Theoretical Upper Bound');

% Plot the Monte Carlo simulation results
plot( ...
    q_values, ...
    simu_result, ...
    'o-', ...
    'Color', [227, 207, 87] / 255, ...
    'LineWidth', 4, ...
    'MarkerSize', 12, ...
    'MarkerFaceColor', 'y', ...
    'DisplayName', 'Simulation Results');

%% Axis labels
xlabel( ...
    'Probability q', ...
    'FontSize', 22, ...
    'FontWeight', 'bold');

ylabel( ...
    'Structural Diagonalizability Probability', ...
    'FontSize', 22, ...
    'FontWeight', 'bold');

%% Axis settings
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 22;
ax.LineWidth = 1.8;
ax.GridAlpha = 0.4;
ax.LabelFontSizeMultiplier = 1.1;

xlim([0, 1]);
ylim([0.4, 1.02]);

xticks(0:0.1:1);
yticks(0.4:0.1:1);

box on;

%% Legend settings
legend( ...
    'Location', 'southeast', ...
    'FontSize', 22, ...
    'FontName', 'Times New Roman');

hold off;

%% Print the complete simulation results
fprintf('\nSimulation summary for G(n,p,q):\n');
fprintf('Number of vertices: n = %d\n', n);
fprintf('Fixed parameter: c = %.1f\n', c_fixed);
fprintf('Non-self-loop probability: p = %.8f\n', p);
fprintf('Number of Monte Carlo samples: %d\n\n', num_samples);

fprintf(['q-value\t\tLower bound\t\tUpper bound\t\t', ...
         'Simulation probability\n']);
fprintf(['------------------------------------------------', ...
         '------------------------\n']);

for q_idx = 1:num_q
    fprintf('%.1f\t\t%.4f\t\t\t%.4f\t\t\t%.4f\n', ...
        q_values(q_idx), ...
        theoretical_value_low(q_idx), ...
        theoretical_value_up(q_idx), ...
        simu_result(q_idx));
end

%% Save the high-resolution figure

% Windows Desktop path
desktop_path = fullfile(getenv('USERPROFILE'), 'Desktop');

% Save the figure to the current directory if the Desktop folder
% cannot be found.
if ~exist(desktop_path, 'dir')
    desktop_path = pwd;
end

output_file = fullfile( ...
    desktop_path, ...
    'gnpq_probability_vs_q_c0.png');

exportgraphics( ...
    h_figure, ...
    output_file, ...
    'Resolution', 600);

fprintf('\nThe figure has been saved to:\n%s\n', output_file);