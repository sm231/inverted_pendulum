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
f = 5;           % 2.5 5 10 20 
Wc = f*2*pi;            % Cut-off frequency [rad/s]


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

%% Closed Loop model with LQR

%%%%%%%%%%%%%%%%%%%%%%% First Closed loop system experiments %%%%%%%%%%%%%%
% Q and R given by KUL
%Q = diag([0.25 4 0 0]); R = 0.003;
%K = lqr(A, B, Q, R);
%sys_cl = ss(A - B*K, B, eye(4), zeros(4,1)); 

%% Effect of Q11
t = 0:0.001:5;
x_ref = 0.1*ones(size(t));

Q11s = [0.1 0.25 0.5 1 2 5 10 20];

% Columns: Q1, SettlingTime, max(abs(u))
performance_Q11 = zeros(length(Q11s),3);

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

figure(2); clf;
hold on; grid on;
ylabel('$u$ [V]','Interpreter','latex')
xlabel('Time [s]')

for i = 1:length(Q11s)

    Q11 = Q11s(i);

    Q = diag([Q11 4 0 0]);
    R = 0.003;

    K = lqr(A, B, Q, R);

    Acl = A - B*K;
    Bcl = B*K(1);

    sys_cl = ss(Acl, Bcl, eye(4), zeros(4,1));

    [~,t,x] = lsim(sys_cl, x_ref, t);

    x_des = [x_ref; zeros(3,length(t))];
    u = -K*(x' - x_des);
    u = u';

    % Step response characteristics for cart position x
    info = stepinfo(x(:,1), t, x_ref(end));

    performance_Q11(i,:) = [ ...
        Q11, ...
        info.SettlingTime, ...
        max(abs(u))
    ];

    figure(1)
    plot(ax1,t,x(:,1))
    plot(ax2,t,x(:,2))
    plot(ax3,t,x(:,3))
    plot(ax4,t,x(:,4))

    figure(2)
    plot(t,u)
end

figure(1)
legend(ax1,string(Q11s),'Location','bestoutside')
title(tl,'Effect of changing $Q_{11}$','Interpreter','latex')

figure(2)
legend(string(Q11s),'Location','best')
title('Control action for different $Q_{11}$','Interpreter','latex')

% Display results
performance_Q11


% poles
figure(3); clf;
hold on; grid on;

xlabel('Real axis')
ylabel('Imaginary axis')

title('Closed-loop poles for different $Q_{11}$','Interpreter','latex')

colors = lines(length(Q11s));

for i = 1:length(Q11s)

    Q11 = Q11s(i);

    Q = diag([Q11 4 0 0]);
    R = 0.003;

    K = lqr(A,B,Q,R);

    Acl = A - B*K;

    p = eig(Acl);

    plot(real(p),imag(p), ...
        'x', ...
        'Color',colors(i,:), ...
        'MarkerSize',10, ...
        'LineWidth',2);
end
legend("", "", "Q11 = " + string(Q11s),'Location','bestoutside')


saveas(figure(1), 'figures/q11_states.png')
saveas(figure(2), 'figures/q11_voltage.png')
saveas(figure(3), 'figures/q11_poles.png')

%% Effect of Q22
t = 0:0.001:5;
x_ref = 0.1*ones(size(t));

Q22s = [1 4 8 16 32];

% Columns: Q2, SettlingTime, max(abs(u))
performance_Q22 = zeros(length(Q22s),3);

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

figure(2); clf;
hold on; grid on;
ylabel('$u$ [V]','Interpreter','latex')
xlabel('Time [s]')

for i = 1:length(Q22s)

    Q22 = Q22s(i);

    Q = diag([0.25 Q22 0 0]);
    R = 0.003;

    K = lqr(A, B, Q, R);

    Acl = A - B*K;
    Bcl = B*K(1);

    sys_cl = ss(Acl, Bcl, eye(4), zeros(4,1));

    [~,t,x] = lsim(sys_cl, x_ref, t);

    x_des = [x_ref; zeros(3,length(t))];
    u = -K*(x' - x_des);
    u = u';

    % Step response characteristics for angle of the rod alpha
    info = stepinfo(x(:,1), t, x_ref(end));

    performance_Q22(i,:) = [ ...
        Q22, ...
        info.SettlingTime, ...
        max(abs(u))
    ];

    figure(1)
    plot(ax1,t,x(:,1))
    plot(ax2,t,x(:,2))
    plot(ax3,t,x(:,3))
    plot(ax4,t,x(:,4))

    figure(2)
    plot(t,u)
end

figure(1)
legend(ax1,string(Q22s),'Location','bestoutside')
title(tl,'Effect of changing $Q_{22}$','Interpreter','latex')

figure(2)
legend(string(Q22s),'Location','best')
title('Control action for different $Q_{22}$','Interpreter','latex')

% Display results
performance_Q22


% Poles
figure(3); clf;
hold on; grid on;

xlabel('Real axis')
ylabel('Imaginary axis')

title('Closed-loop poles for different $Q_{22}$','Interpreter','latex')

colors = lines(length(Q22s));

for i = 1:length(Q22s)

    Q22 = Q22s(i);

    Q = diag([0.25 Q22 0 0]);
    R = 0.003;

    K = lqr(A,B,Q,R);

    Acl = A - B*K;

    p = eig(Acl);

    plot(real(p),imag(p), ...
        'x', ...
        'Color',colors(i,:), ...
        'MarkerSize',10, ...
        'LineWidth',2);
end

legend("Q22 = " + string(Q22s),'Location','bestoutside')

saveas(figure(1), 'figures/q22_states.png')
saveas(figure(2), 'figures/q22_voltage.png')
saveas(figure(3), 'figures/q22_poles.png')

%% Effect of Q33

t = 0:0.001:5;
x_ref = 0.1*ones(size(t));

Q33s = [0 0.1 0.5 1 2 5 10];

% Columns: Q33, SettlingTime, MaxControlVoltage
performance_Q33 = zeros(length(Q33s),3);

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

figure(2); clf;
hold on; grid on;
ylabel('$u$ [V]','Interpreter','latex')
xlabel('Time [s]')

for i = 1:length(Q33s)

    Q33 = Q33s(i);

    Q = diag([0.25 4 Q33 0]);
    R = 0.003;

    K = lqr(A, B, Q, R);

    Acl = A - B*K;
    Bcl = B*K(1);

    sys_cl = ss(Acl, Bcl, eye(4), zeros(4,1));

    [~,t,x] = lsim(sys_cl, x_ref, t);

    x_des = [x_ref; zeros(3,length(t))];
    u = -K*(x' - x_des);
    u = u';

    info = stepinfo(x(:,1), t, x_ref(end));

    performance_Q33(i,:) = [ ...
        Q33, ...
        info.SettlingTime, ...
        max(abs(u)) ...
    ];

    figure(1)
    plot(ax1,t,x(:,1))
    plot(ax2,t,x(:,2))
    plot(ax3,t,x(:,3))
    plot(ax4,t,x(:,4))

    figure(2)
    plot(t,u)
end

figure(1)
legend(ax1,string(Q33s),'Location','bestoutside')
title(tl,'Effect of changing $Q_{33}$','Interpreter','latex')

figure(2)
legend(string(Q33s),'Location','best')
title('Control action for different $Q_{33}$','Interpreter','latex')

performance_Q33


% Poles
figure(3); clf;
hold on; grid on;

xlabel('Real axis')
ylabel('Imaginary axis')

title('Closed-loop poles for different $Q_{33}$','Interpreter','latex')

colors = lines(length(Q33s));

for i = 1:length(Q33s)

    Q33 = Q33s(i);

    Q = diag([0.25 4 Q33 0]);
    R = 0.003;

    K = lqr(A,B,Q,R);

    Acl = A - B*K;

    p = eig(Acl);

    plot(real(p),imag(p), ...
        'x', ...
        'Color',colors(i,:), ...
        'MarkerSize',10, ...
        'LineWidth',2);
end

legend("Q33 = " + string(Q33s),'Location','bestoutside')

saveas(figure(1), 'figures/q33_states.png')
saveas(figure(2), 'figures/q33_voltage.png')
saveas(figure(3), 'figures/q33_poles.png')


%% Effect of Q44

t = 0:0.001:5;
x_ref = 0.1*ones(size(t));

Q44s = [0 0.1 0.5 1 2 5 10];

% Columns: Q44, SettlingTime, MaxControlVoltage
performance_Q44 = zeros(length(Q44s),3);

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

figure(2); clf;
hold on; grid on;
ylabel('$u$ [V]','Interpreter','latex')
xlabel('Time [s]')

for i = 1:length(Q44s)

    Q44 = Q44s(i);

    Q = diag([0.25 4 0 Q44]);
    R = 0.003;

    K = lqr(A, B, Q, R);

    Acl = A - B*K;
    Bcl = B*K(1);

    sys_cl = ss(Acl, Bcl, eye(4), zeros(4,1));

    [~,t,x] = lsim(sys_cl, x_ref, t);

    x_des = [x_ref; zeros(3,length(t))];
    u = -K*(x' - x_des);
    u = u';

    info = stepinfo(x(:,1), t, x_ref(end));

    performance_Q44(i,:) = [ ...
        Q44, ...
        info.SettlingTime, ...
        max(abs(u)) ...
    ];

    figure(1)
    plot(ax1,t,x(:,1))
    plot(ax2,t,x(:,2))
    plot(ax3,t,x(:,3))
    plot(ax4,t,x(:,4))

    figure(2)
    plot(t,u)
end

figure(1)
legend(ax1,string(Q44s),'Location','bestoutside')
title(tl,'Effect of changing $Q_{44}$','Interpreter','latex')

figure(2)
legend(string(Q44s),'Location','best')
title('Control action for different $Q_{44}$','Interpreter','latex')

performance_Q44


% Poles
figure(3); clf;
hold on; grid on;

xlabel('Real axis')
ylabel('Imaginary axis')

title('Closed-loop poles for different $Q_{44}$','Interpreter','latex')

colors = lines(length(Q44s));

for i = 1:length(Q44s)

    Q44 = Q44s(i);

    Q = diag([0.25 4 0 Q44]);
    R = 0.003;

    K = lqr(A,B,Q,R);

    Acl = A - B*K;

    p = eig(Acl);

    plot(real(p),imag(p), ...
        'x', ...
        'Color',colors(i,:), ...
        'MarkerSize',10, ...
        'LineWidth',2);
end

legend("Q44 = " + string(Q44s),'Location','bestoutside')

saveas(figure(1), 'figures/q44_states.png')
saveas(figure(2), 'figures/q44_voltage.png')
saveas(figure(3), 'figures/q44_poles.png')

%% Effect of R

t = 0:0.001:5;
x_ref = 0.1*ones(size(t));

Rs = [0.001 0.003 0.005 0.01 0.03 0.05];

% Columns: R, SettlingTime, MaxControlVoltage
performance_R = zeros(length(Rs),3);

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

figure(2); clf;
hold on; grid on;
ylabel('$u$ [V]','Interpreter','latex')
xlabel('Time [s]')

for i = 1:length(Rs)

    R = Rs(i);

    Q = diag([0.25 4 0 0]);

    K = lqr(A, B, Q, R);

    Acl = A - B*K;
    Bcl = B*K(1);

    sys_cl = ss(Acl, Bcl, eye(4), zeros(4,1));

    [~,t,x] = lsim(sys_cl, x_ref, t);

    x_des = [x_ref; zeros(3,length(t))];
    u = -K*(x' - x_des);
    u = u';

    info = stepinfo(x(:,1), t, x_ref(end));

    performance_R(i,:) = [ ...
        R, ...
        info.SettlingTime, ...
        max(abs(u)) ...
    ];

    figure(1)
    plot(ax1,t,x(:,1))
    plot(ax2,t,x(:,2))
    plot(ax3,t,x(:,3))
    plot(ax4,t,x(:,4))

    figure(2)
    plot(t,u)
end

figure(1)
legend(ax1,string(Rs),'Location','bestoutside')
title(tl,'Effect of changing $R$','Interpreter','latex')

figure(2)
legend(string(Rs),'Location','best')
title('Control action for different $R$','Interpreter','latex')

performance_R


% Poles
figure(3); clf;
hold on; grid on;

xlabel('Real axis')
ylabel('Imaginary axis')

title('Closed-loop poles for different $R$','Interpreter','latex')

colors = lines(length(Rs));

for i = 1:length(Rs)

    R = Rs(i);

    Q = diag([0.25 4 0 0]);

    K = lqr(A,B,Q,R);

    Acl = A - B*K;

    p = eig(Acl);

    plot(real(p),imag(p), ...
        'x', ...
        'Color',colors(i,:), ...
        'MarkerSize',10, ...
        'LineWidth',2);
end

legend("R = " + string(Rs),'Location','bestoutside')

saveas(figure(1), 'figures/r_states.png')
saveas(figure(2), 'figures/r_voltage.png')
saveas(figure(3), 'figures/r_poles.png')