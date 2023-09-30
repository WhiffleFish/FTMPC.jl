%% Plot state trajectories
clc; clear; close;

% Space
x = [-2, 2];
y = [-2, 2];
z = [-9.5, 5];
reference = [-.4 .4 -5];

% Load mode 7 data
load('nominal7.mat')
load('nonrobust7.mat')
load('unit7.mat')

% Define plot parameters
rotor_fail_7 = dft_nominal_7(floor(length(dft_nominal_7)/4));
imm_delay_7 = dft_nominal_7(floor(length(dft_nominal_7)/4) + 3);
reference_time = dft_nominal_7;

% Create figure
figure(1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First subfigure
subplot(3, 1, 1);
h1 = plot(dft_nonrobust_7(1:length(dfx_nonrobust_7)), dfx_nonrobust_7(1,:),'-', ...
    'LineWidth', 2, 'Color', "blue");
hold on
h2 = plot(dft_unit_7(1:length(dfx_unit_7)), dfx_unit_7(1,:), '-', ...
    'LineWidth', 2, 'Color', "magenta");

h3 = plot(dft_nominal_7(1:length(dfx_nominal_7)), dfx_nominal_7(1,:), '-', ...
    'LineWidth', 2, 'Color', "cyan");
ylabel("x [m]")

plot([rotor_fail_7,rotor_fail_7], x, '--', 'LineWidth', 1, 'Color', "black"); 
hold on
plot([imm_delay_7, imm_delay_7], x, '--', 'LineWidth', 1, 'Color', "green");
plot(reference_time , reference(1)*ones(1, length(reference_time )), '--', ...
    'LineWidth', 1, 'Color', "red");
grid on

% Rectangle box
x_box = [0, reference_time(end)];
y_box = [-2, -1];
rectangle('Position', [x_box(1), y_box(1), diff(x_box), diff(y_box)],...
    'FaceColor', [0.85, 0.85, 0.85], 'EdgeColor', 'k');
y_box = [1, 2];
rectangle('Position', [x_box(1), y_box(1), diff(x_box), diff(y_box)],...
    'FaceColor', [0.85, 0.85, 0.85], 'EdgeColor', 'k');

plot(dft_unit_7(1:length(dfx_unit_7)), dfx_unit_7(1,:), '-', ...
    'LineWidth', 2, 'Color', "magenta")
plot(dft_nonrobust_7(1:length(dfx_nonrobust_7)), dfx_nonrobust_7(1,:),'-', ...
    'LineWidth', 2, 'Color', "blue")
plot(dft_nominal_7(1:length(dfx_nominal_7)), dfx_nominal_7(1,:), '-', ...
    'LineWidth', 2, 'Color', "cyan")
ylabel("x [m]")

plot([rotor_fail_7,rotor_fail_7], x, '--', 'LineWidth', 1, 'Color', "black"); 
plot([imm_delay_7, imm_delay_7], x, '--', 'LineWidth', 1, 'Color', "green");

plot(dft_nonrobust_7(length(dfx_nonrobust_7)), dfx_nonrobust_7(1, end), ...
    'rx', 'MarkerSize', 10, 'LineWidth', 2);
plot(dft_unit_7(length(dfx_unit_7)), dfx_unit_7(1, end), ...
    'rx', 'MarkerSize', 10, 'LineWidth', 2);
xlim(x_box)
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3, 1, 2);
plot([rotor_fail_7,rotor_fail_7], y, '--', 'LineWidth', 1, 'Color', "black"); 
hold on
plot([imm_delay_7, imm_delay_7], y, '--', 'LineWidth', 1, 'Color', "green");
plot(reference_time , reference(2)*ones(1, length(reference_time)), '--', ...
    'LineWidth', 1, 'Color', "red");
grid on
plot(dft_unit_7(1:length(dfx_unit_7)), dfx_unit_7(2,:), '-', ...
    'LineWidth', 2, 'Color', "magenta")
plot(dft_nonrobust_7(1:length(dfx_nonrobust_7)), dfx_nonrobust_7(2,:),'-', ...
    'LineWidth', 2, 'Color', "blue")
plot(dft_nominal_7(1:length(dfx_nominal_7)), dfx_nominal_7(2,:), '-', ...
    'LineWidth', 2, 'Color', "cyan")
ylim(y)
ylabel("y [m]")

x_box = [0, reference_time(end)];
y_box = [-2, -1];
rectangle('Position', [x_box(1), y_box(1), diff(x_box), diff(y_box)],...
    'FaceColor', [0.85, 0.85, 0.85], 'EdgeColor', 'k');
y_box = [1, 7];
rectangle('Position', [x_box(1), y_box(1), diff(x_box), diff(y_box)],...
    'FaceColor', [0.85, 0.85, 0.85], 'EdgeColor', 'k');

plot([rotor_fail_7,rotor_fail_7], y, '--', 'LineWidth', 1, 'Color', "black"); 
plot([imm_delay_7, imm_delay_7], y, '--', 'LineWidth', 1, 'Color', "green");

plot(dft_nonrobust_7(length(dfx_nonrobust_7)), dfx_nonrobust_7(2, end), ...
    'rx', 'MarkerSize', 10, 'LineWidth', 2);
plot(dft_unit_7(length(dfx_unit_7)), dfx_unit_7(2, end), ...
    'rx', 'MarkerSize', 10, 'LineWidth', 2);
xlim(x_box)
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3, 1, 3);
plot([rotor_fail_7,rotor_fail_7], z, '--', 'LineWidth', 1, 'Color', "black"); 
hold on
plot([imm_delay_7, imm_delay_7], z, '--', 'LineWidth', 1, 'Color', "green");
plot(reference_time , reference(3)*ones(1, length(reference_time )), '--', ...
    'LineWidth', 1, 'Color', "red");
grid on
plot(dft_unit_7(1:length(dfx_unit_7)), -dfx_unit_7(3,:), '-', ...
    'LineWidth', 2, 'Color', "magenta")
plot(dft_nonrobust_7(1:length(dfx_nonrobust_7)), -dfx_nonrobust_7(3,:),'-', ...
    'LineWidth', 2, 'Color', "blue")
plot(dft_nominal_7(1:length(dfx_nominal_7)), -dfx_nominal_7(3,:), '-', ...
    'LineWidth', 2, 'Color', "cyan")
ylim(z)
ylabel("z [m]")

x_box = [0, reference_time(end)];
y_box = [-9.5, -6];
rectangle('Position', [x_box(1), y_box(1), diff(x_box), diff(y_box)],...
    'FaceColor', [0.85, 0.85, 0.85], 'EdgeColor', 'k');
plot([rotor_fail_7,rotor_fail_7], z, '--', 'LineWidth', 1, 'Color', "black"); 
plot([imm_delay_7, imm_delay_7], z, '--', 'LineWidth', 1, 'Color', "green");

plot(dft_nonrobust_7(length(dfx_nonrobust_7)), -dfx_nonrobust_7(3, end), ...
    'rx', 'MarkerSize', 10, 'LineWidth', 2);
plot(dft_unit_7(length(dfx_unit_7)), -dfx_unit_7(3, end), ...
    'rx', 'MarkerSize', 10, 'LineWidth', 2);
xlim(x_box)
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlabel("Time (s) ")
labels = {'\color{black} Non-Robust', ...
          '\color{black} First-Step Consensus', ...
          '\color{black} FGMPC'};
legend([h1, h2, h3], labels, 'FontSize', 8, 'FontName', 'cmr12',...
    'TextColor', 'black', "Position", [0.8, 0.9, 0, 0]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



print(gcf,'sevenmode.png','-dpng','-r500')



