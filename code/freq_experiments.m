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

%% Cutoff frequency 

Q = diag([5 8 1 3]); R = 0.03;
K = lqr(A, B, Q, R);


f = 50;           % 0.1 0.3 0.5 1 2 5 10 20 50
Wc = f*2*pi;            % Cut-off frequency [rad/s]

%%


%states = out.real_states;
%save('real_freq_data/real_freq_50_states.mat', 'states')
%voltage = out.real_volt;
%save('real_freq_data/real_freq_50_volt.mat', 'voltage')

%%

states = open("real_freq_data/real_freq_0_1_states.mat"); 
voltage = open("real_freq_data/real_freq_0_1_volt.mat"); 
states = states.states;
voltage = voltage.voltage;

real_x = reshape(states.Data, 4, numel(t));

real_u = reshape(voltage.Data, 1, numel(t));

t = 0:Ts:5;
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

%saveas(figure(1), '../report/figures/real_freq_states_01.png')
%saveas(figure(2), '../report/figures/real_freq_volt_01.png');

%%

frequencies = {'0_3', '0_5', '1', '2', '5', '10'};
t = 0:Ts:5;
num_states = 4;
labels = {'$x$ [m]', '$\alpha$ [rad]', '$\dot{x}$ [m/s]', '$\dot{\alpha}$ [rad/s]'};

fig1 = figure(1); clf; 
tl = tiledlayout(4,1); 

fig2 = figure(2); clf; 
hold on; grid on;
ylabel('$u$ [V]','Interpreter','latex')
xlabel('Time [s]')

ax = [];
for i = 1:num_states
    ax(i) = nexttile(tl);
    hold(ax(i), 'on'); grid(ax(i), 'on');
    ylabel(ax(i), labels{i}, 'Interpreter', 'latex')
end
xlabel(tl, 'Time [s]')

for f = 1:length(frequencies)
    freq_str = frequencies{f};
    
    state_file = sprintf('real_freq_data/real_freq_%s_states.mat', freq_str);
    volt_file  = sprintf('real_freq_data/real_freq_%s_volt.mat', freq_str);
    
    s_data = load(state_file);
    v_data = load(volt_file);
    
    current_x = reshape(s_data.states.Data, num_states, []);
    current_u = reshape(v_data.voltage.Data, 1, []);
    
    display_name = [strrep(freq_str, '_', '.'), ' Hz'];
    
    for i = 1:num_states
        plot(ax(i), t, current_x(i, :), 'DisplayName', display_name);
    end
    
    figure(2);
    plot(t, current_u, 'DisplayName', display_name);

end

figure(1)
legend(ax(1), 'show', 'Location', 'northeastoutside');

figure(2);
legend('show', 'Location', 'northeast');

%saveas(figure(1), '../report/figures/real_freq_states.png')
%saveas(figure(2), '../report/figures/real_freq_volt.png');

%%

states = open("real_freq_data/real_freq_50_states.mat"); 
voltage = open("real_freq_data/real_freq_50_volt.mat"); 
states = states.states;
voltage = voltage.voltage;

real_x = reshape(states.Data, 4, numel(t));

real_u = reshape(voltage.Data, 1, numel(t));

t = 0:Ts:5;
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

saveas(figure(1), '../report/figures/real_freq_states_50.png')
saveas(figure(2), '../report/figures/real_freq_volt_50.png');