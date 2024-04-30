clc
clear all;
close all;

% Define the network and initial conditions
N = 5; % Number of nodes
x0 = [1; 2; 3; 4; 5]; % Initial states of nodes
G = 0.75*[0 1 0 0 1; 1 0 1 0 0; 0 1 0 0 1; 0 0 0 0 1; 1 0 1 1 0]; % Adjacency matrix
A = double(G); % Coupling matrix
k_max = 100; % Number of iterations
epsilon = 1/3; % Step size

% Initialization
x_alpha = x0; % Initialize alpha states
x_beta = -randn*ones(N,1); % Initialize beta states (could be random)
x = zeros(N, k_max); % Record of state evolution
z = zeros(1, k_max);

% Dynamics simulation with observer update
for k = 1:(k_max - 1)
    x_alpha_temp = x_alpha;
    for i = 1:N
        neighbors = find(A(i,:) > 0); % Find neighbors
        if isempty(neighbors)
            continue;  % Skip if no neighbors
        end
        % Calculate the sum of the differences
        sum_diff = sum(x_alpha_temp(neighbors) - x_alpha_temp(i));
        
        % Update x_alpha for node i
        x_alpha(i) = x_alpha_temp(i) + epsilon * sum_diff;
    end
    
    x(:,k + 1) = x_alpha + x_beta; % Total state is the sum of alpha and beta
    
    z = zeros(1, k_max);
    z(1) = 0;

    for k = 1:(k_max - 1)
        sum_weighted_diffs = sum(A(1, :) .* (x(:, k) - x(1, k)).');
        z(k + 1) = z(k) + x(1, k + 1) - (x(1, k) + epsilon * sum_weighted_diffs);
    end
end

% Plot results
figure;
plot(1:k_max, x');
title('State Evolution Over Time');
xlabel('Time Step');
ylabel('State Value');
legend('Node 1', 'Node 2', 'Node 3', 'Node 4', 'Node 5');
hold on;
plot(1:k_max, z, 'r*-', 'DisplayName', 'z[k]');
hold on;
x1_initial = x0(1) * ones(1, k_max);
plot(1:k_max, x1_initial, 'k--', 'DisplayName', 'x1[0]');
hold off;
title('State Evolution Over Time with Reference Lines');
xlabel('Time Step');
ylabel('State Value');
legend('show');


%% Fig 6
noise_levels = 0.1:0.05:0.9;

simulations = 100;
Est_Err = zeros(length(noise_levels), simulations);
Avg_Err = zeros(length(noise_levels), simulations);

for noise_idx = 1:length(noise_levels)
    noise_level = noise_levels(noise_idx);
    
    for sim = 1:simulations
        x_alpha = x0;
        x = zeros(N, k_max);
        for k = 1:k_max
            x_alpha_temp = x_alpha;
            for i = 1:N
                neighbors = find(A(i,:) > 0);
                sum_diff = sum(x_alpha_temp(neighbors) - x_alpha_temp(i));
                noise = noise_level * randn;
                x_alpha(i) = x_alpha_temp(i) + epsilon * sum_diff + noise;
            end
            x(:,k) = x_alpha + x_beta;
        end
        final_states = x(:, end);
        theoretical_average = mean(x0);
        Est_Err(noise_idx, sim) = norm((x0 - final_states),2);
        Avg_Err(noise_idx, sim) = norm((final_states - theoretical_average * ones(N, 1)),2);
    end
end

mean_Est_Err = mean(Est_Err, 2);
std_Est_Err = std(Est_Err, 0, 2);
mean_Avg_Err = mean(Avg_Err, 2);
std_Avg_Err = std(Avg_Err, 0, 2);

% Plot results
figure;
errorbar(noise_levels, mean_Avg_Err, std_Avg_Err, 'DisplayName', 'Avg Err');
hold on;
errorbar(noise_levels, mean_Est_Err, std_Est_Err, 'DisplayName', 'Est Err');
hold off;
legend('show');
xlabel('Noise level c');
ylabel('Error');
title('Effect of Noise on Privacy and Accuracy with Error Bars');