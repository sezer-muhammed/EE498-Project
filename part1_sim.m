% Define parameters
mu = 1/82.45;
mu_star = 1 - mu;

% Define initial conditions
x1_0 = 1.2;
x2_0 = 0;
x3_0 = 0;
x4_0 = -1.0494;

% Define simulation time span
tspan = [0, 10]; % Adjust the end time as needed

% Define the state space equations
f = @(t, x) [
    x(2);
    2 * x(4) + x(1) - mu_star * (x(1) + mu) / sqrt((x(1) + mu)^2 + x(3)^2)^3 - mu * (x(1) - mu_star) / sqrt((x(1) - mu_star)^2 + x(3)^2)^3;
    x(4);
    - 2 * x(2) + x(3) - mu_star * x(3) / sqrt((x(1) + mu)^2 + x(3)^2)^3 - mu * x(3) / sqrt((x(1) - mu_star)^2 + x(3)^2)^3
];

% Define different configurations to compare
configurations = {
    struct('MaxStep', 'auto', 'RelTol', 1e-5, 'Solver', @ode45, 'DisplayName', 'Configuration 1'),
    struct('MaxStep', 'auto', 'RelTol', 1e-2, 'Solver', @ode45, 'DisplayName', 'Configuration 2'),
    struct('MaxStep', 'auto', 'RelTol', 1e-5, 'Solver', @ode23, 'DisplayName', 'Configuration 3'),
    struct('MaxStep', 0.1, 'RelTol', 1e-5, 'Solver', @ode45, 'DisplayName', 'Configuration 4 (Fixed Step Size)')
};

figure;
hold on;

% Simulate and plot different configurations
for i = 1:numel(configurations)
    config = configurations{i};

    % Set simulation parameters
    options = odeset('MaxStep', config.MaxStep, 'RelTol', config.RelTol);
    
    % Solve the state space equations
    [t, x] = config.Solver(f, tspan, [x1_0, x2_0, x3_0, x4_0], options);

    % Extract x and y coordinates
    x_coordinates = x(:, 1);
    y_coordinates = x(:, 3);

    % Plot the trajectory
    plot(x_coordinates, y_coordinates, 'DisplayName', config.DisplayName);
end

hold off;
xlabel('x');
ylabel('y');
title('Apollo Satellite Trajectory');
legend('Location', 'best');
