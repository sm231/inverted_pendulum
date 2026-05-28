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


%% Individual plot

path = "data/Q3_cfreq_10_";

% Load Position Data
p_data = load(path + "p.mat");
tp = p_data.data.Time;
pos = reshape(p_data.data.Data, 1, []);

% Load Angle Data
a_data = load(path + "a.mat");
ta = a_data.data.Time;
ang = reshape(a_data.data.Data, 1, []);

figure(1); clf;
tiledlayout(2,1);

% Position
nexttile;
stairs(tp, pos, 'LineWidth', 1.5);
grid on;
ylabel('$x$ [m]');

% Angle
nexttile;
stairs(ta, ang, 'LineWidth', 1.5);
grid on;
ylabel('$\alpha$ [rad]');
xlabel('Time [s]');


% Load Voltage
v_data = load(path + "v.mat");
tv = v_data.data.Time;
volts = reshape(v_data.data.Data, 1, []);

figure(2); clf;
stairs(tv, volts, 'r', 'LineWidth', 1.5);
grid on;
ylabel('$u$ [V]');
xlabel('Time [s]');



%% cutted plot

cut_start = 6; 
cut_end   = 27.5;

p_data = load(path + "p.mat");
tp = p_data.data.Time;
pos = reshape(p_data.data.Data, 1, []);

a_data = load(path + "a.mat");
ta = a_data.data.Time;
ang = reshape(a_data.data.Data, 1, []);

mask_p = (tp >= cut_start) & (tp <= cut_end);
mask_a = (ta >= cut_start) & (ta <= cut_end);

figure(1); clf;
tiledlayout(2,1);

% Tile 1: Position
nexttile;
stairs(tp(mask_p), pos(mask_p), 'LineWidth', 1.5);
grid on; ylabel('$x$ [m]');
xlim([cut_start cut_end]); % Forces the axis to match exactly

% Tile 2: Angle
nexttile;
stairs(ta(mask_a), ang(mask_a), 'LineWidth', 1.5);
grid on; ylabel('$\alpha$ [rad]'); xlabel('Time [s]');
xlim([cut_start cut_end]);


v_data = load(path + "v.mat");
tv = v_data.data.Time;
volts = reshape(v_data.data.Data, 1, []);

% Create mask for voltage time
mask_v = (tv >= cut_start) & (tv <= cut_end);

figure(2); clf;
stairs(tv(mask_v), volts(mask_v), 'r', 'LineWidth', 1.5);
grid on; ylabel('$u$ [V]'); xlabel('Time [s]');
xlim([cut_start cut_end]);


saveas(figure(1), '../report/figures/setup_cfreq_10_states.png')
saveas(figure(2), '../report/figures/setup_cfreq_10_volt.png');