% Define Constants 
T_eng = 0.46;
K_eng = @(t) 0.732; % Assuming delta K(t) is zero
F_s = @(s) 1.5*s/(s^2 + 3*s + 4);
T_brk = 0.193;
K_brk = @(t) 0.979; 
m = 1; % Specify mass
r_accel = 1; % Specify r_accel
a_thr_off = 1; % Specify a_thr_off
T_hw = 1; % Specify T_hw
rho_accel = 1/(1 + r_accel/m);

% Discretization Parameters
Ts = 0.01; % Sampling time
Tfinal = 10; % Final simulation time
N = Tfinal/Ts; % Number of time steps

% Initialization
x = zeros(3,N); % State matrix
y = zeros(3,N); % Output matrix
u = zeros(1,N); % Control signal
t = 0:Ts:Tfinal-Ts; % Time vector

x(:,1) = [0; 2; 3]; % Replace with your own initial condition
% Constants
Q = 100 * eye(3);
Q(1,1) = 0
R = 1;

% Control loop (RHC)
for k = 1:N-1
    % System Equations
    Af = @(u) (u >= a_thr_off) * -1/T_eng + (u < a_thr_off) * -1/T_brk;
    Bf = @(u) (u >= a_thr_off) * K_eng(t(k))/T_eng + (u < a_thr_off) * K_brk(t(k))/T_brk;
    A = [0, 1, -T_hw; 0, 0, -1; 0, 0, Af(u(k))];
    B = [0; 0; rho_accel*Bf(u(k))];
    
    % Euler's Method (Discretization)
    x(:,k+1) = x(:,k) + Ts * (A * x(:,k) + B * u(k));
    
    % Output Equation
    y(:,k) = [1, 0, 0; 0, 1, 0; 0, 0, 1] * x(:,k);
    
    % RHC Controller
    u_min = -2; % Specify the minimum control input
    u_max = 2; % Specify the maximum control input
    nu = 100; % Number of control inputs to test
    U = linspace(u_min, u_max, nu); % Discretize the control input
    J = zeros(1, nu); % Cost function
    
    for i = 1:nu
        % Predicted state
        x_pred = x(:,k+1) + Ts * (A * x(:,k+1) + B * U(i));
        
        % Cost function
        J(i) = x_pred'*Q*x_pred + U(i)'*R*U(i);
    end
    
    % Optimal control action
    [~, idx] = min(J);
    u(k+1) = U(idx);
end

% Plot results
figure;
subplot(2,1,1);
plot(t, x);
title('States');
subplot(2,1,2);
plot(t, y);
title('Outputs');

