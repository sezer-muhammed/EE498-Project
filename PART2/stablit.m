clear all;
close all;

%% Define parameters
T_eng = 0.46;
T_brk = 0.193;
K_c = 0.732;
K_brk = 0.979;
a_thr_off = 0; % set your threshold value here
T_hw = 0; % set your T_hw value here
dt = 0.01; % define your time step
t = 0:dt:10; % define your time range

%% Initialize system states and output
x = zeros(3, length(t));
x(3, 1) = 0.1;
y = zeros(1, length(t));
u = zeros(1, length(t));

%% Define filter for Delta K(s)
Fs = tf([1.5, 0], [1, 3, 4]);
Fs_discrete = c2d(Fs, dt, 'tustin'); % convert continuous filter to discrete filter
[filter_num, filter_den] = tfdata(Fs_discrete, 'v'); % get filter coefficients

%% Define average system parameters
A_avg = [0, 1, 0; 0, 0, 1; 0, 0, -((1/T_eng + 1/T_brk)/2)]; 
B_avg = [0; 0; ((K_c/T_eng + K_brk/T_brk)/2)];
C_avg = [1, 0, 0]; 
D_avg = 0;

%% Create average system
sys_avg = ss(A_avg, B_avg, C_avg, D_avg, dt); % dt is your sampling time

%% Initialize MPC Controller using average system
MPCobj = mpc(sys_avg);
MPCobj.Model.Nominal.U = 0; % assuming that the nominal input is 0
MPCobj.Model.Nominal.Y = 0; % assuming that the nominal output is 0



%% Initialize the MPC state
xmpc = mpcstate(MPCobj); 

%% Time variant simulation
for i=2:length(t)
    % Compute Delta K(t)
    if i == 2
        delta_K = filter(filter_num, filter_den, [u(i-1), 0]);
    else
        delta_K = filter(filter_num, filter_den, [u(i-1), 0], delta_K);
    end
    
    % Compute A_f(t) and B_f(t)
    if u(i-1) >= a_thr_off
        A_f = -1 / T_eng;
        B_f = (K_c + delta_K(end)) / T_eng;
    else
        A_f = -1 / T_brk;
        B_f = K_brk / T_brk;
    end
    
    % Create state-space matrices
    A = [0, 1, 0; 0, 0, 1; 0, 0, A_f];
    B = [0; 0; B_f];
    C = [1, 0, 0]; 
    D = 0; % assuming D=0 as it's not specified
    
    % Update system states
    x(:, i) = x(:, i-1) + dt * (A * x(:, i-1) + B * u(i-1));
    
    % Update the MPC state
    xmpc.Plant = x(:, i);
    
    % Compute system output
    y(i) = C * x(:, i) + D * u(i-1);

    % Apply MPC control
    u(i) = mpcmove(MPCobj, xmpc, y(i), 0); % Assuming reference signal is 0
end


%% Plot output and states
figure;
subplot(4,1,1);
plot(t, y);
xlabel('Time (s)');
ylabel('Output y(t)');
title('System Output over Time');
grid on;

subplot(4,1,2);
plot(t, x(1,:));
xlabel('Time (s)');
ylabel('State x1');
title('State x1 over Time');
grid on;

subplot(4,1,3);
plot(t, x(2,:));
xlabel('Time (s)');
ylabel('State x2');
title('State x2 over Time');
grid on;

subplot(4,1,4);
plot(t, x(3,:));
xlabel('Time (s)');
ylabel('State x3');
title('State x3 over Time');
grid on;
