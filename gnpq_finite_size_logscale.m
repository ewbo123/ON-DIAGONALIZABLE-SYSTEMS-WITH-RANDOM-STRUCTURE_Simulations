clear; clc; close all;

%% Global font settings

set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultTextFontName', 'Times New Roman');
set(groot, 'defaultLegendFontName', 'Times New Roman');
set(groot, 'defaultAxesFontSize', 14);
set(groot, 'defaultTextFontSize', 14);

%% Parameter settings

% Fixed parameter c
c_fixed = 0;

% Selected self-loop probabilities
selected_q = [0, 0.2, 0.4, 0.6, 0.8, 1.0];
num_selected_q = length(selected_q);

% Graph sizes used in the simulations
n_values = 500:500:4000;
num_n = length(n_values);

% Number of Monte Carlo samples for each pair (q,n)
num_samples = 800;

% Fix the random seed for reproducibility
rng(1);

% Check whether munkres.m is available
if exist('munkres', 'file') ~= 2
    error(['The function munkres.m was not found. ', ...
           'Please add munkres.m to the current folder or MATLAB path.']);
end

% Preallocate memory for simulation results
sim_probabilities = zeros(num_selected_q, num_n);

%% Main simulation loop

for q_idx = 1:num_selected_q

    q = selected_q(q_idx);

    fprintf('\n====================================================\n');
    fprintf('Processing q = %.1f\n', q);
    fprintf('====================================================\n');

    for n_idx = 1:num_n

        n = n_values(n_idx);

        % Probability of each non-self-loop directed edge
        %
        % In MATLAB, log(n) denotes the natural logarithm.
        p = (log(n) + c_fixed) / n;

        % Ensure that p is a valid probability
        p = max(0, min(p, 1));

        fprintf('\nGraph size n = %d\n', n);
        fprintf('Non-self-loop probability p = %.8f\n', p);
        fprintf('Self-loop probability q = %.1f\n', q);

        count = 0;

        %% Monte Carlo simulation

        for sample = 1:num_samples

            % Generate all entries initially with probability p.
            % The diagonal entries are subsequently replaced using
            % the prescribed self-loop probability q.
            A = rand(n) < p;

            % Generate the diagonal entries independently with
            % self-loop probability q.
            diagonal_support = rand(n, 1) < q;
            A(1:n+1:end) = diagonal_support;

            % Convert the adjacency matrix to the system-matrix convention
            A_sys = A.';
            clear A;

            % Compute the generic rank of the zero-nonzero pattern
            generic_rank = sprank(sparse(A_sys));

            % Construct the cost matrix for minimum-cost matching:
            %
            % Existing edges:              cost 0
            % Missing diagonal edges:      cost 1
            % Missing non-diagonal edges:  sufficiently large cost
            large_cost = n + 1;
            A_cost = large_cost * ones(n);

            % Existing edges have zero cost
            A_cost(A_sys) = 0;

            % Missing diagonal entries have unit cost
            diagonal_indices = logical(eye(n));
            A_cost(diagonal_indices) = ...
                1 - double(diag(A_sys));

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
                         sample, ...
                         num_samples, ...
                         count / sample);
            end
        end

        % Empirical probability of structural diagonalizability
        sim_probabilities(q_idx, n_idx) = ...
            count / num_samples;

        fprintf(['Completed q = %.1f, n = %d: ', ...
                 'simulation probability = %.4f\n'], ...
                 q, ...
                 n, ...
                 sim_probabilities(q_idx, n_idx));

        % Save a checkpoint after each completed (q,n) pair
        save( ...
            'gnpq_logscale_checkpoint.mat', ...
            'c_fixed', ...
            'selected_q', ...
            'n_values', ...
            'num_samples', ...
            'sim_probabilities');
    end
end

%% Theoretical upper and lower bounds for G(n,p,q)

% Theoretical lower bound
theoretical_low = ...
    exp(-2 .* (1 - selected_q) .* exp(-c_fixed)) .* ...
    (1 ...
    + 2 .* (1 - selected_q) .* exp(-c_fixed) ...
    + (1 - selected_q).^2 .* exp(-2 * c_fixed));

% Theoretical upper bound
theoretical_up = ...
    1 ...
    - (1 - selected_q).^2 .* exp(-2 * c_fixed) .* ...
      exp(-2 .* (1 - selected_q) .* exp(-c_fixed));

%% Colors, markers, and logarithmic tick settings

% Blue, red, yellow, green, purple, and orange
colors = [
    0.0000, 0.4470, 0.7410;
    0.8500, 0.3250, 0.0980;
    0.9290, 0.6940, 0.1250;
    0.4660, 0.6740, 0.1880;
    0.4940, 0.1840, 0.5560;
    1.0000, 0.5000, 0.0000
];

markers = {'o', 's', '^', 'd', 'v', '*'};

% Representative ticks for the logarithmic horizontal axis.
%
% The ratio between adjacent displayed values is 2, so these
% tick positions are equally spaced on a logarithmic axis.
tick_values = [500, 1000, 2000, 4000];

tick_labels = arrayfun( ...
    @(x) sprintf('%d', x), ...
    tick_values, ...
    'UniformOutput', false);

%% Create the six-panel figure

h_subplot = figure( ...
    'Position', [100, 100, 1800, 1120], ...
    'Color', 'w');

t = tiledlayout(h_subplot, 2, 3);

t.TileSpacing = 'compact';
t.Padding = 'loose';

for q_idx = 1:num_selected_q

    ax = nexttile(t);

    hold(ax, 'on');
    grid(ax, 'on');

    q = selected_q(q_idx);
    current_color = colors(q_idx, :);

    %% Theoretical shaded range

    x_fill = [
        min(n_values), ...
        max(n_values), ...
        max(n_values), ...
        min(n_values)
    ];

    y_fill = [
        theoretical_low(q_idx), ...
        theoretical_low(q_idx), ...
        theoretical_up(q_idx), ...
        theoretical_up(q_idx)
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

    %% Monte Carlo simulation curve

    plot( ...
        ax, ...
        n_values, ...
        sim_probabilities(q_idx, :), ...
        [markers{q_idx}, '-'], ...
        'Color', current_color, ...
        'LineWidth', 2.6, ...
        'MarkerSize', 8, ...
        'MarkerFaceColor', current_color, ...
        'DisplayName', 'Simulation Results');

    %% Use a true logarithmic horizontal axis

    set(ax, 'XScale', 'log');

    %% Panel title

    title( ...
        ax, ...
        sprintf( ...
            'q = %.1f (Theory: %.3f--%.3f)', ...
            q, ...
            theoretical_low(q_idx), ...
            theoretical_up(q_idx)), ...
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

    % All eight data points remain in the figure.
    % Only four representative logarithmic tick labels are shown.
    xticks(ax, tick_values);
    xticklabels(ax, tick_labels);
    xtickangle(ax, 0);

    % Remove repeated axis labels from individual panels
    xlabel(ax, '');
    ylabel(ax, '');

    legend( ...
        ax, ...
        'Location', 'southeast', ...
        'FontSize', 14, ...
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
    'Structural Diagonalizability Probability', ...
    'FontSize', 18, ...
    'FontWeight', 'bold', ...
    'FontName', 'Times New Roman');

%% Print the complete simulation results

fprintf('\nSimulation summary for G(n,p,q):\n');
fprintf('Fixed parameter c: %.1f\n', c_fixed);
fprintf('Number of Monte Carlo samples: %d\n', num_samples);

fprintf('Graph sizes: ');
fprintf('%d ', n_values);
fprintf('\n');

fprintf('Self-loop probabilities: ');
fprintf('%.1f ', selected_q);
fprintf('\n\n');

fprintf('q-value');

for n_idx = 1:num_n
    fprintf('\tn=%d', n_values(n_idx));
end

fprintf('\n');

for q_idx = 1:num_selected_q

    fprintf('%.1f', selected_q(q_idx));

    for n_idx = 1:num_n
        fprintf('\t%.4f', ...
            sim_probabilities(q_idx, n_idx));
    end

    fprintf('\n');
end

%% Print theoretical bounds

fprintf('\nTheoretical upper and lower bounds for G(n,p,q):\n');
fprintf('Fixed parameter c = %.1f\n', c_fixed);
fprintf('q-value\tLower bound\tUpper bound\n');
fprintf('----------------------------------------\n');

for q_idx = 1:num_selected_q
    fprintf( ...
        '%.1f\t\t%.4f\t\t%.4f\n', ...
        selected_q(q_idx), ...
        theoretical_low(q_idx), ...
        theoretical_up(q_idx));
end

%% Save the six-panel figure only

desktop_path = fullfile(getenv('USERPROFILE'), 'Desktop');

% Save to the current folder if the Desktop folder cannot be found
if ~exist(desktop_path, 'dir')
    desktop_path = pwd;
end

output_file = fullfile( ...
    desktop_path, ...
    'subplot_gnpq_true_logscale.png');

exportgraphics( ...
    h_subplot, ...
    output_file, ...
    'Resolution', 600);

fprintf('\nThe six-panel figure has been saved to:\n%s\n', ...
    output_file);