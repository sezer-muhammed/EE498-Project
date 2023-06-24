clc;
close all;

%% Variables
mu = 1/82.45;
mu_star = 1 - mu;

%% Runge-Kutta (p=4)
h = 0.001;   % Step size
time = 0:h:10;  % Time span

x = zeros(4, length(time)); % Rows are states, columns represent time

x(1,1) = 1.2;
x(2,1) = 0;
x(3,1) = 0;
x(4,1) = -1.0494;

F_x_one = @(t,x_two) x_two;
F_x_two = @(t,x_one,x_three,x_four) 2*x_four + x_one - mu_star*(x_one+mu)./sqrt(((x_one+mu)^2+x_three^2).^3) - mu*(x_one-mu_star)./sqrt(((x_one-mu_star)^2+x_three^2).^3);
F_x_three = @(t,x_four) x_four;
F_x_four = @(t,x_one,x_two,x_three) -2*x_two + x_three - mu_star*x_three./sqrt(((x_one+mu)^2+x_three^2).^3) - mu*x_three./sqrt(((x_one-mu_star)^2+x_three^2).^3);

k_matrix = zeros(4,4); % k matrix (rows are states, columns are k indices)

for k = 1:(length(time)-1)
    k_matrix(1,1) = h * F_x_one(time(k) , x(2,k));
    k_matrix(2,1) = h * F_x_two(time(k) , x(1,k) , x(3,k) , x(4,k));
    k_matrix(3,1) = h * F_x_three(time(k) , x(4,k));
    k_matrix(4,1) = h * F_x_four(time(k) , x(1,k), x(2,k), x(3,k));

    k_matrix(1,2) = h * F_x_one(time(k) + h/2   , x(2,k) + k_matrix(2,1)/2);
    k_matrix(2,2) = h * F_x_two(time(k) + h/2   , x(1,k) + k_matrix(1,1)/2, x(3,k) + k_matrix(3,1)/2 , x(4,k) + k_matrix(4,1)/2);
    k_matrix(3,2) = h * F_x_three(time(k) + h/2 , x(4,k) + k_matrix(4,1)/2);
    k_matrix(4,2) = h * F_x_four(time(k) + h/2  , x(1,k) + k_matrix(1,1)/2, x(2,k) + k_matrix(2,1)/2 , x(3,k) + k_matrix(3,1)/2);

    k_matrix(1,3) = h * F_x_one(time(k) + h/2   , x(2,k) + k_matrix(2,2)/2);
    k_matrix(2,3) = h * F_x_two(time(k) + h/2   , x(1,k) + k_matrix(1,2)/2, x(3,k) + k_matrix(3,2)/2 , x(4,k) + k_matrix(4,2)/2);
    k_matrix(3,3) = h * F_x_three(time(k) + h/2 , x(4,k) + k_matrix(4,2)/2);
    k_matrix(4,3) = h * F_x_four(time(k) + h/2  , x(1,k) + k_matrix(1,2)/2, x(2,k) + k_matrix(2,2)/2 , x(3,k) + k_matrix(3,2)/2);

    k_matrix(1,4) = h * F_x_one(time(k) + h/2   , x(2,k) + k_matrix(2,3));
    k_matrix(2,4) = h * F_x_two(time(k) + h/2   , x(1,k) + k_matrix(1,3), x(3,k) + k_matrix(3,3) , x(4,k) + k_matrix(4,3));
    k_matrix(3,4) = h * F_x_three(time(k) + h/2 , x(4,k) + k_matrix(4,3));
    k_matrix(4,4) = h * F_x_four(time(k) + h/2  , x(1,k) + k_matrix(1,3), x(2,k) + k_matrix(2,3) , x(3,k) + k_matrix(3,3));

    x(1,k+1) = x(1,k) + (k_matrix(1,1) + 2*k_matrix(1,2) + 2*k_matrix(1,3) + k_matrix(1,4))/6;
    x(2,k+1) = x(2,k) + (k_matrix(2,1) + 2*k_matrix(2,2) + 2*k_matrix(2,3) + k_matrix(2,4))/6;
    x(3,k+1) = x(3,k) + (k_matrix(3,1) + 2*k_matrix(3,2) + 2*k_matrix(3,3) + k_matrix(3,4))/6;
    x(4,k+1) = x(4,k) + (k_matrix(4,1) + 2*k_matrix(4,2) + 2*k_matrix(4,3) + k_matrix(4,4))/6;
end


%% Variable Step Size
MAX_ITER = 9999;
atol = 1e-2 * ones(4,1);
rtol = 1e-2 * ones(4,1);
h = 0.0001;
h_min = 0.000001;
time = 0;
tf = 10;
if time + h > tf
    h = tf - time;
end

beta = 0.95;
fac_zero = 0.2;
fac_one = 1.5;
step_sizes = zeros(1,MAX_ITER);
x = zeros(4,MAX_ITER); % Rows are states, columns represent time
x_hat = zeros(4,1);
eta = zeros(4,1);
x(1,1) = 1.2;
x(2,1) = 0;
x(3,1) = 0;
x(4,1) = -1.0494;

k_matrix = zeros(4,4); % k matrix (rows are states, columns are k indices)
k = 1;

for ITER = 1:MAX_ITER

    k_matrix(1,1) = h * F_x_one(time , x(2,k));
    k_matrix(2,1) = h * F_x_two(time , x(1,k) , x(3,k) , x(4,k));
    k_matrix(3,1) = h * F_x_three(time , x(4,k));
    k_matrix(4,1) = h * F_x_four(time , x(1,k), x(2,k), x(3,k));

    k_matrix(1,2) = h * F_x_one(time + h/2   , x(2,k) + k_matrix(2,1)/2);
    k_matrix(2,2) = h * F_x_two(time + h/2   , x(1,k) + k_matrix(1,1)/2, x(3,k) + k_matrix(3,1)/2 , x(4,k) + k_matrix(4,1)/2);
    k_matrix(3,2) = h * F_x_three(time + h/2 , x(4,k) + k_matrix(4,1)/2);
    k_matrix(4,2) = h * F_x_four(time + h/2  , x(1,k) + k_matrix(1,1)/2, x(2,k) + k_matrix(2,1)/2 , x(3,k) + k_matrix(3,1)/2);

    k_matrix(1,3) = h * F_x_one(time + h/2   , x(2,k) + k_matrix(2,2)/2);
    k_matrix(2,3) = h * F_x_two(time + h/2   , x(1,k) + k_matrix(1,2)/2, x(3,k) + k_matrix(3,2)/2 , x(4,k) + k_matrix(4,2)/2);
    k_matrix(3,3) = h * F_x_three(time + h/2 , x(4,k) + k_matrix(4,2)/2);
    k_matrix(4,3) = h * F_x_four(time + h/2  , x(1,k) + k_matrix(1,2)/2, x(2,k) + k_matrix(2,2)/2 , x(3,k) + k_matrix(3,2)/2);

    k_matrix(1,4) = h * F_x_one(time + h/2   , x(2,k) + k_matrix(2,3));
    k_matrix(2,4) = h * F_x_two(time + h/2   , x(1,k) + k_matrix(1,3), x(3,k) + k_matrix(3,3) , x(4,k) + k_matrix(4,3));
    k_matrix(3,4) = h * F_x_three(time + h/2 , x(4,k) + k_matrix(4,3));
    k_matrix(4,4) = h * F_x_four(time + h/2  , x(1,k) + k_matrix(1,3), x(2,k) + k_matrix(2,3) , x(3,k) + k_matrix(3,3));

    x_hat(1) = x(1,k) + (k_matrix(1,1) + 2*k_matrix(1,2) + 2*k_matrix(1,3) + k_matrix(1,4))/6;
    x_hat(2) = x(2,k) + (k_matrix(2,1) + 2*k_matrix(2,2) + 2*k_matrix(2,3) + k_matrix(2,4))/6;
    x_hat(3) = x(3,k) + (k_matrix(3,1) + 2*k_matrix(3,2) + 2*k_matrix(3,3) + k_matrix(3,4))/6;
    x_hat(4) = x(4,k) + (k_matrix(4,1) + 2*k_matrix(4,2) + 2*k_matrix(4,3) + k_matrix(4,4))/6;

    eta(1) = (k_matrix(1,2) - 2*k_matrix(1,3) + k_matrix(1,4))/6;
    eta(2) = (k_matrix(2,2) - 2*k_matrix(2,3) + k_matrix(2,4))/6;
    eta(3) = (k_matrix(3,2) - 2*k_matrix(3,3) + k_matrix(3,4))/6;
    eta(4) = (k_matrix(4,2) - 2*k_matrix(4,3) + k_matrix(4,4))/6;

    sigma =  sqrt(sum((abs(eta)./(atol + rtol.*abs(x_hat))).^2 )/4 );

    h_bar = h* min(fac_one,max(fac_zero,beta*sigma^(-0.25)));

    if sigma <= 1
        if time == tf
            x(:,k+1) = x_hat;
            break
        end
        time = time + h;
        x(:,k+1) = x_hat;
        if time + h_bar > tf
            h_bar = tf - time;
        end
        step_sizes(k) = h;
        k = k + 1;
    end

    h = h_bar;
    if h < h_min
        break
    end

end

min_step_size_idx = find(step_sizes(1:k-1) == min(step_sizes(1:k-1)));

% Plot trajectory with variable step size
figure;
plot(x(1,1:k), x(3,1:k), 'LineWidth', 1.5);
title('Variable Step Size Simulation');
xlabel('X Coordinate (meters)');
ylabel('Y Coordinate (meters)');
grid on;
axis equal;

% Mark specific points on the trajectory
hold on;

% Mark minimum step size locations
plot(x(1,min_step_size_idx), x(3,min_step_size_idx), 'b*', 'MarkerSize', 10);


% Plot velocity versus time
velocity = sqrt(x(2,1:k).^2 + x(4,1:k).^2);
figure;
plot(1:k, velocity, 'LineWidth', 1.5);
title('Variable Step Size Simulation: Velocity versus Time');
xlabel('Number of Iterations');
ylabel('Velocity (meters/second)');
grid on;

% Mark specific points on the velocity plot
hold on;

% Plot step size versus number of iterations
figure;
plot(1:k-1, step_sizes(1:k-1), 'LineWidth', 1.5);
title('Step Size versus Number of Iterations');
xlabel('Number of Iterations');
ylabel('Step Size (seconds)');
grid on;
hold on;
plot(min_step_size_idx, step_sizes(min_step_size_idx), 'b*', 'MarkerSize', 10);
