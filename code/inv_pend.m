clear variables

%%

set(groot, 'defaultAxesTickLabelInterpreter', 'latex'); 
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');

set(groot, 'defaultAxesFontSize', 20);

% Set global Font Size for Text (includes X/Y Labels and Titles)
set(groot, 'defaultTextFontSize', 20);

% Set global Font Size for Legends
set(groot, 'defaultLegendFontSize', 14);

set(groot, 'defaultLineLineWidth', 2);

%% This linese of code are used to reveale all labels in the Simulink diagrams for plotting

hb = find_system(gcs,'Type','Block');
handles = cell2mat(get_param(hb,'Handle'));
arrayfun(@(h) set_param(h,'ShowName','on'), handles);
arrayfun(@(h) set_param(h,'ShowName','on','HideAutomaticName','off'), handles);

%% Parameters

% Inverted Pendulum system
Rm = 2.6;           % Motor armature resistance [Ohm]
Km = 0.00767;       % Motor Torque constant [Nm/A]
Kb = 0.00767;       % Motor back EMF constant [V/(rad/s)]
Kg = 3.7;           % Motor gear ratio
M = 0.455;          % Cart mass [kg]
l = 0.305;          % Rod length (l = l_tot/2) [m]
m = 0.210;          % Rod mass [kg]
r = 0.635 * 10^-2;  % Radius of motor output gear [m]
g = 9.81;           % Acceleration due to gravity [m/s^2]


% Filter
Ts = 0.005;         % Sampling time [5ms]
f = 5;              % 2.5 5 10 20 [Hz]  !!! 5 is used in real setup !!!
Wc = f*2*pi;        % Cut-off frequency [rad/s]


% Quantizer
Fx = 0.456/4.41;    % Conversion factor from volts to meters [m/V]
Fa = (pi/2)/6.328;  % Conversion factor from volts to radians [rad/V]          

%% State Space

A = [0, 0, 1, 0;
     0, 0, 0, 1;
     0, 0, 0, 0];

A(3,2) = (-m*g)/M;        A(3,3) = (-Kg^2*Km*Kb)/(M*Rm*r^2);
A(4,2) = ((M+m)*g)/(M*l); A(4,3) = (Kg^2*Km*Kb)/(M*Rm*r^2*l);


B = [0; 0; (Km*Kg)/(M*Rm*r); (-Km*Kg)/(r*Rm*M*l)];

C = [1, 0, 0, 0;
     0, 1, 0, 0];

sys = ss(A, B, C, [0; 0]);

%% Open loop analysis

% 1) Stability 
eigen = eig(A); % there exists pole in the RHP

% 2) Controllability
Co = ctrb(A, B);
rank(Co); % The system is controllable

% 3) Observability
Ob = obsv(A, C);
rank(Ob); % The system is observable

% 4) Transmission zero
trans_zeros = tzero(sys); %no transmission zeros

%% Closed Loop model with LQR

%%%%%%%%%%%%%%%%%%%%%%% First Closed loop system experiments %%%%%%%%%%%%%%
% Q and R obtained by trial and error experiments
Q = diag([5 8 0 0]); R = 0.03;

K = lqr(A, B, Q, R);

sys_cl = ss(A - B*K, B*K(1), eye(4), zeros(4,1));  % BKxd -> BK(1) with a reference signal 

% See the lqr_tuning.m for the plots of the tuning experiments

% performance
 
t = 0:0.001:5;              % Simulation time
x_ref = 0.1*ones(size(t));  % setpoint

% Columns: SettlingTime, maxangle, time of maxangle, maxControlVoltage
performance = zeros(1,4);


[~,t,x] = lsim(sys_cl, x_ref, t); % simulation 

x_des = [x_ref; zeros(3,length(t))]; % desired state 
u = -K*(x' - x_des); % control law
u = u';

% Step response characteristics for cart position x
infox = stepinfo(x(:,1), t, x_ref(end));

performance(:) = [ ...
    infox.SettlingTime, ...
    max(x(:,2)), ...
    t(x(:,2) == max(x(:,2))), ...
    max(abs(u)), ...
]

figure(1); clf;
tl = tiledlayout(4,1);

ax1 = nexttile; hold on; grid on;
ylabel('$x$ [m]','Interpreter','latex')

ax2 = nexttile; hold on; grid on;
ylabel('$\alpha$ [rad]','Interpreter','latex')

ax3 = nexttile; hold on; grid on;
ylabel('$\dot{x}$ [m/s]','Interpreter','latex')

ax4 = nexttile; hold on; grid on;
ylabel('$\dot{\alpha}$ [rad/s]','Interpreter','latex')
xlabel('Time [s]')

plot(ax1,t,x(:,1))
plot(ax2,t,x(:,2))
plot(ax3,t,x(:,3))
plot(ax4,t,x(:,4))

figure(2); clf;
hold on; grid on;
ylabel('$u$ [V]','Interpreter','latex')
xlabel('Time [s]')
plot(t,u)

%saveas(figure(1), '../report/figures/step_states_3.png')
%saveas(figure(2), '../report/figures/step_volt_3.png')




%% Sensor error
% Data obtained from the sensors on the setup

a = open('angle_sensor_data.mat');
a = a.angle_sensor_data;
a_data = a.Data;
a_time = a.Time;

figure
plot(a_time, a_data);
xlabel '$t$';
ylabel '$V$';

angle_noise_var = var(a_data);


p = open('position_sensor_data.mat');
p = p.position_sensor_data;
p_data = p.Data;
p_time = p.Time;

figure
plot(p_time, p_data);
xlabel '$t$';
ylabel '$V$';

position_noise_var = var(p_data);

%%

%% Plots for report section "realistic closed loop model"

% Q and R from closed loop tuning
% real_states1 
% real_volt1


Q = diag([5 8 1 3]); R = 0.03;
K = lqr(A, B, Q, R);

f = 2; Wc = f*2*pi; 

t = 0:Ts:5;

%states = out.real_states;
%save('real_states2.mat', 'states')
%voltage = out.real_volt;
%save('real_volt2.mat', 'voltage')

% 1 -> LQR with Q1 and R1
% 2 -> LQR with Q2 and R2
states = open("real_states2.mat"); % choose 1 or 2 for correct dataset
voltage = open("real_volt2.mat"); % choose 1 or 2 for correct dataset
states = states.states;
voltage = voltage.voltage;

real_x = reshape(states.Data, 4, numel(t));

real_u = reshape(voltage.Data, 1, numel(t));


figure(1); clf;
tl = tiledlayout(4,1);

ax1 = nexttile; hold on; grid on;
ylabel('$x$ [m]','Interpreter','latex')

ax2 = nexttile; hold on; grid on;
ylabel('$\alpha$ [rad]','Interpreter','latex')

ax3 = nexttile; hold on; grid on;
ylabel('$\dot{x}$ [m/s]','Interpreter','latex')

ax4 = nexttile; hold on; grid on;
ylabel('$\dot{\alpha}$ [rad/s]','Interpreter','latex')
xlabel('Time [s]')

plot(ax1,t,real_x(1,:))
plot(ax2,t,real_x(2,:))
plot(ax3,t,real_x(3,:))
plot(ax4,t,real_x(4,:))

figure(2); clf;
hold on; grid on;
ylabel('$u$ [V]','Interpreter','latex')
xlabel('Time [s]')
plot(t,real_u)

% Step response characteristics for cart position x
infox = stepinfo(real_x(1,:), t, 0, 0.05);

performance(:) = [ ...
    infox.SettlingTime, ...
    max(real_x(2,:)), ...
    t(real_x(2,:) == max(real_x(2,:))), ...
    max(abs(real_u)), ...
]

saveas(figure(1), '../report/figures/real_states2.png')
saveas(figure(2), '../report/figures/real_volt2.png');

%%
%%%%%%%%%%%%%%%%%%% Realistic Closed loop system experiments %%%%%%%%%%%%%%
% All tests were done with f = 5;

f = 2; Wc = f*2*pi; 

% Q_base (best) f = 5;
Q = diag([10 15 2 6]); R = 0.03;       

% Q_slow
%Q = diag([3 8 1 3]); R = 0.03; 

% Q_fast_position
%Q = diag([30 15 2 6]); R = 0.03;  

% Q_strong_angle
%Q = diag([10 50 2 8]); R = 0.03;  

% Q_more_damping
%Q = diag([10 15 8 20]); R = 0.03; 

% R_high
%Q = diag([10 15 2 6]); R = 0.8; 

% R_low
%Q = diag([10 15 2 6]); R = 0.003; 


% LQR                  
K = lqr(A, B, Q, R);
    
% Closed loop system
sys_cl = ss(A - B*K, B, eye(4), zeros(4,1));




%% Performance

% Step response of the closed loop system
[x, t] = step(sys_cl, 5); 

subplot(2,1,1);
plot(t, x(:,1)); title('Cart Position (m)'); ylabel('meters');
grid on;
subplot(2,1,2);
plot(t, x(:,2)); title('Pendulum Angle (rad)'); ylabel('radians');
grid on;

% Performance: RiseTime, SettlingTime, PeakTime, Overshoot, ...
perf = stepinfo(sys_cl);

% Voltage
u = zeros(length(t), 1);
for i = 1:length(t)
    u(i) = -K * x(i,:)';
end

figure;
plot(t, u, 'b', 'LineWidth', 2);
hold on; grid on;
xlabel 'Time [sec]';
ylabel 'Voltage [V]';

% Closed loop poles 
poles_cl = eig(A-B*K)
pzmap(sys_cl)

% log how the dominant poles evolve as the Q and R matrices are tweaked 
% Same with performance stats of the step response

%% Realistic Simulation

