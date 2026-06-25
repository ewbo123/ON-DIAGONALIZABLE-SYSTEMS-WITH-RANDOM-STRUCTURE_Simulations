clear; clc; close all;

%% Global font settings
set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultTextFontName', 'Times New Roman');
set(groot, 'defaultLegendFontName', 'Times New Roman');
set(groot, 'defaultAxesFontSize', 14);
set(groot, 'defaultTextFontSize', 14);

%% Parameter settings

% Fixed values of c
selected_c = [-0.5, 0, 0.5, 1.0, 1.5, 2.0];
num_selected_c = length(selected_c);

% Actual graph sizes used in the simulations
%
% These values do not need to be equally spaced.
% Their positions in the figure are determined by the logarithmic axis.
n_values = [100, 500, 1000, 3000, 6000, 10000];
num_n = length(n_values);

% Number of Monte Carlo samples for each pair (c,n)
num_samples = 800;

% Fix the random seed for reproducibility
rng(1);

% Check whether munkres.m is available
if exist('munkres', 'file') ~= 2
    error(['The function munkres.m was not found. ', ...
           'Please add munkres.m to the current folder or MATLAB path.']);
end

% Preallocate memory for the simulation results
sim_probabilities = zeros(num_selected_c, num_n);

%% Main simulation loop

for c_idx = 1:num_selected_c
    c = selected_c(c_idx);

    fprintf('\n====================================================\n');
    fprintf('Processing c = %.2f\n', c);
    fprintf('====================================================\n');

    for n_idx = 1:num_n
        n = n_values(n_idx);

        % Edge probability
        % In MATLAB, log(n) denotes the natural logarithm.
        p = (log(n) + c) / n;

        % Ensure that p is a valid probability
        p = max(0, min(p, 1));

        fprintf('\nGraph size n = %d, edge probability p = %.8f\n', ...
            n, p);

        count = 0;

        %% Monte Carlo simulation
        for sample = 1:num_samples

            % Generate a random adjacency matrix from G(n,p).
            % Every possible directed edge, including self-loops,
            % is independently present with probability p.
            A = rand(n) < p;

            % Convert the adjacency matrix to the system matrix convention
            A_sys = A.';
            clear A;

            % Compute the generic rank of the zero-nonzero pattern
            generic_rank = sprank(sparse(A_sys));

            % Construct the cost matrix for minimum-cost matching
            %
            % Existing edges: cost 0
            % Missing diagonal edges: cost 1
            % Other missing edges: sufficiently large cost
            large_cost = n + 1;
            A_cost = large_cost * ones(n);
            A_cost(A_sys) = 0;

            diagonal_indices = logical(eye(n));
            A_cost(diagonal_indices) = 1 - double(diag(A_sys));

            % Compute a minimum-cost perfect matching
            [~, min_totalcost] = munkres(A_cost);

            % Check the structural diagonalizability condition
            if generic_rank == n - min_totalcost
                count = count + 1;
            end

            % Display intermediate progress
            if mod(sample, 50) == 0 || sample == num_samples
                fprintf(['  Completed sample %d/%d, ', ...
                         'current probability = %.4f\n'], ...
                         sample, num_samples, count / sample);
            end
        end

        % Empirical probability of structural diagonalizability
        sim_probabilities(c_idx, n_idx) = count / num_samples;

        fprintf(['Completed c = %.2f, n = %d: ', ...
                 'simulation probability = %.4f\n'], ...
                 c, n, sim_probabilities(c_idx, n_idx));

        % Save a checkpoint after each graph size
        save( ...
            'gnp_logscale_checkpoint.mat', ...
            'selected_c', ...
            'n_values', ...
            'num_samples', ...
            'sim_probabilities');
    end
end

%% Theoretical upper and lower bounds

theoretical_low = ...
    exp(-2 .* exp(-selected_c)) .* ...
    (1 ...
    + 2 .* exp(-selected_c) ...
    + exp(-2 .* selected_c));

theoretical_up = ...
    1 ...
    - exp(-2 .* exp(-selected_c)) .* ...
      exp(-2 .* selected_c);

%% Colors and markers

colors = [
    0.0000, 0.4470, 0.7410;
    0.8500, 0.3250, 0.0980;
    0.9290, 0.6940, 0.1250;
    0.4660, 0.6740, 0.1880;
    0.4940, 0.1840, 0.5560;
    1.0000, 0.5000, 0.0000
];

markers = {'o', 's', '^', 'd', 'v', '*'};

%% Create the six-panel figure

h_subplot = figure( ...
    'Position', [100, 100, 1800, 1120], ...
    'Color', 'w');

t = tiledlayout(h_subplot, 2, 3);
t.TileSpacing = 'compact';
t.Padding = 'loose';

for c_idx = 1:num_selected_c
    ax = nexttile(t);

    hold(ax, 'on');
    grid(ax, 'on');

    c = selected_c(c_idx);
    current_color = colors(c_idx, :);

    % Use a true logarithmic horizontal axis
    ax.XScale = 'log';

    %% Theoretical shaded range

    x_fill = [
        min(n_values), ...
        max(n_values), ...
        max(n_values), ...
        min(n_values)
    ];

    y_fill = [
        theoretical_low(c_idx), ...
        theoretical_low(c_idx), ...
        theoretical_up(c_idx), ...
        theoretical_up(c_idx)
    ];

    fill( ...
        ax, ...
        x_fill, ...
        y_fill, ...
        current_color, ...
        'FaceAlpha', 0.20, ...
        'EdgeColor', current_color, ...
        'LineStyle', '--', ...
        'LineWidth', 1.2, ...
        'DisplayName', 'Theoretical Range');

    %% Simulation curve

    plot( ...
        ax, ...
        n_values, ...
        sim_probabilities(c_idx, :), ...
        [markers{c_idx}, '-'], ...
        'Color', current_color, ...
        'LineWidth', 2.6, ...
        'MarkerSize', 8, ...
        'MarkerFaceColor', current_color, ...
        'DisplayName', 'Simulation Results');

    %% Panel title

    title( ...
        ax, ...
        sprintf( ...
            'c = %.1f (Theory: %.3f - %.3f)', ...
            c, ...
            theoretical_low(c_idx), ...
            theoretical_up(c_idx)), ...
        'FontSize', 15, ...
        'FontWeight', 'bold', ...
        'FontName', 'Times New Roman');

    %% Axis settings

    ax.FontName = 'Times New Roman';
    ax.FontSize = 12.5;
    ax.FontWeight = 'bold';
    ax.LineWidth = 1.5;
    ax.GridAlpha = 0.30;

    ax.XMinorTick = 'on';
    ax.XMinorGrid = 'on';
    ax.YMinorGrid = 'on';

    xlim(ax, [min(n_values), max(n_values)]);
    ylim(ax, [0, 1.05]);

    % Display the actual simulated graph sizes as tick labels.
    % Their positions are still determined by the logarithmic scale.
    xticks(ax, n_values);
    xticklabels(ax, string(n_values));
    xtickangle(ax, 0);

    % Remove repeated labels from individual panels
    xlabel(ax, '');
    ylabel(ax, '');

    legend( ...
        ax, ...
        'Location', 'southeast', ...
        'FontSize', 12, ...
        'FontName', 'Times New Roman', ...
        'FontWeight', 'bold');

    box(ax, 'on');
    hold(ax, 'off');
end

%% Shared axis labels

xlabel( ...
    t, ...
    'Number of Vertices n (log scale)', ...
    'FontSize', 18, ...
    'FontWeight', 'bold', ...
    'FontName', 'Times New Roman');

ylabel( ...
    t, ...
    'Probability of Structural Diagonalizability', ...
    'FontSize', 18, ...
    'FontWeight', 'bold', ...
    'FontName', 'Times New Roman');

%% Print simulation results

fprintf('\nSimulation summary for G(n,p):\n');
fprintf('Number of Monte Carlo samples: %d\n', num_samples);
fprintf('Graph sizes: ');
fprintf('%d ', n_values);
fprintf('\n\n');

fprintf('c-value');

for n_idx = 1:num_n
    fprintf('\tn=%d', n_values(n_idx));
end

fprintf('\n');

for c_idx = 1:num_selected_c
    fprintf('%.1f', selected_c(c_idx));

    for n_idx = 1:num_n
        fprintf('\t%.4f', sim_probabilities(c_idx, n_idx));
    end

    fprintf('\n');
end

%% Print theoretical bounds

fprintf('\nTheoretical upper and lower bounds:\n');
fprintf('c-value\tLower bound\tUpper bound\n');
fprintf('----------------------------------------\n');

for c_idx = 1:num_selected_c
    fprintf( ...
        '%.1f\t\t%.4f\t\t%.4f\n', ...
        selected_c(c_idx), ...
        theoretical_low(c_idx), ...
        theoretical_up(c_idx));
end

%% Save the six-panel figure

desktop_path = fullfile(getenv('USERPROFILE'), 'Desktop');

if ~exist(desktop_path, 'dir')
    desktop_path = pwd;
end

output_file = fullfile( ...
    desktop_path, ...
    'subplot_gnp_true_logscale.png');

exportgraphics( ...
    h_subplot, ...
    output_file, ...
    'Resolution', 600);

fprintf('\nThe six-panel figure has been saved to:\n%s\n', output_file);