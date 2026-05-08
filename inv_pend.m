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
Wc = 2;            % Cut-off frequency [Hz]

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

% Assignment suggestions for Q and R. !!! TUNING !!!
Q = diag([4 4 3 3]); % Q(1,1) -> Cart Tracking
R = 0.003;              % Q(2,2) -> Pendulum 
                        % R -> Control Input
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

% Low pass filter
% Backward difference derivatives
% Actuator Saturation
% Quantization

pos = out.position_sys_cl;

t = pos(:, 1);
x = pos(:, 2);

info = stepinfo(x,t);
setTime = info.SettlingTime