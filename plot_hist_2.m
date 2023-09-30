%% Plot state trajectories
clc; clear; close;

% Space
x = [-2, 2];
y = [-2, 2];
z = [-9.5, 5];
reference = [-.4 .4 -5];

% Load mode 2 data
load('nominal2.mat')
load('nonrobust2.mat')
load('unit2.mat')

% Define plot parameters
rotor_fail_2 = dft_nominal_2(floor(length(dft_nominal_2)/4));
imm_delay_2 = dft_nominal_2(floor(length(dft_nominal_2)/4) + 3);
reference_time = dft_nominal_2;

% Create figure
figure(1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First subfigure
subplot(3, 1, 1);
h1 = plot([rotor_fail_2,rotor_fail_2], x, '--', 'LineWidth', 1, 'Color', "black"); 
hold on
h2 = plot([imm_delay_2, imm_delay_2], x, '--', 'LineWidth', 1, 'Color', "green");
h3 = plot(reference_time , reference(1)*ones(1, length(reference_time )), '--', ...
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

plot(dft_unit_2(1:length(dfx_unit_2)), dfx_unit_2(1,:), '-', ...
    'LineWidth', 2, 'Color', "magenta")
plot(dft_nonrobust_2(1:length(dfx_nonrobust_2)), dfx_nonrobust_2(1,:),'-', ...
    'LineWidth', 2, 'Color', "blue")
plot(dft_nominal_2(1:length(dfx_nominal_2)), dfx_nominal_2(1,:), '-', ...
    'LineWidth', 2, 'Color', "cyan")
ylabel("x [m]")

plot([rotor_fail_2,rotor_fail_2], x, '--', 'LineWidth', 1, 'Color', "black"); 
plot([imm_delay_2, imm_delay_2], x, '--', 'LineWidth', 1, 'Color', "green");

plot(dft_nonrobust_2(length(dfx_nonrobust_2)), dfx_nonrobust_2(1, end), ...
    'rx', 'MarkerSize', 10, 'LineWidth', 2);
plot(dft_unit_2(length(dfx_unit_2)), dfx_unit_2(1, end), ...
    'rx', 'MarkerSize', 10, 'LineWidth', 2);
xlim(x_box)
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3, 1, 2);
plot([rotor_fail_2,rotor_fail_2], y, '--', 'LineWidth', 1, 'Color', "black"); 
hold on
plot([imm_delay_2, imm_delay_2], y, '--', 'LineWidth', 1, 'Color', "green");
plot(reference_time , reference(2)*ones(1, length(reference_time )), '--', ...
    'LineWidth', 1, 'Color', "red");
grid on
plot(dft_unit_2(1:length(dfx_unit_2)), dfx_unit_2(2,:), '-', ...
    'LineWidth', 2, 'Color', "magenta")
plot(dft_nonrobust_2(1:length(dfx_nonrobust_2)), dfx_nonrobust_2(2,:),'-', ...
    'LineWidth', 2, 'Color', "blue")
plot(dft_nominal_2(1:length(dfx_nominal_2)), dfx_nominal_2(2,:), '-', ...
    'LineWidth', 2, 'Color', "cyan")
xlim(x_box)
ylim(y)
ylabel("y [m]")

x_box = [0, reference_time(end)];
y_box = [-2, -1];
rectangle('Position', [x_box(1), y_box(1), diff(x_box), diff(y_box)],...
    'FaceColor', [0.85, 0.85, 0.85], 'EdgeColor', 'k');
y_box = [1, 2];
rectangle('Position', [x_box(1), y_box(1), diff(x_box), diff(y_box)],...
    'FaceColor', [0.85, 0.85, 0.85], 'EdgeColor', 'k');

plot([rotor_fail_2,rotor_fail_2], y, '--', 'LineWidth', 1, 'Color', "black"); 
plot([imm_delay_2, imm_delay_2], y, '--', 'LineWidth', 1, 'Color', "green");

plot(dft_nonrobust_2(length(dfx_nonrobust_2)), dfx_nonrobust_2(2, end), ...
    'rx', 'MarkerSize', 10, 'LineWidth', 2);
plot(dft_unit_2(length(dfx_unit_2)), dfx_unit_2(2, end), ...
    'rx', 'MarkerSize', 10, 'LineWidth', 2);

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3, 1, 3);
plot([rotor_fail_2,rotor_fail_2], z, '--', 'LineWidth', 1, 'Color', "black"); 
hold on
plot([imm_delay_2, imm_delay_2], z, '--', 'LineWidth', 1, 'Color', "green");
plot(reference_time , reference(3)*ones(1, length(reference_time )), '--', ...
    'LineWidth', 1, 'Color', "red");
grid on
plot(dft_unit_2(1:length(dfx_unit_2)), -dfx_unit_2(3,:), '-', ...
    'LineWidth', 2, 'Color', "magenta")
plot(dft_nonrobust_2(1:length(dfx_nonrobust_2)), -dfx_nonrobust_2(3,:),'-', ...
    'LineWidth', 2, 'Color', "blue")
plot(dft_nominal_2(1:length(dfx_nominal_2)), -dfx_nominal_2(3,:), '-', ...
    'LineWidth', 2, 'Color', "cyan")
ylim(z)

ylabel("z [m]")

x_box = [0, reference_time(end)];
y_box = [-9.5, -6];
rectangle('Position', [x_box(1), y_box(1), diff(x_box), diff(y_box)],...
    'FaceColor', [0.85, 0.85, 0.85], 'EdgeColor', 'k');
plot([rotor_fail_2,rotor_fail_2], z, '--', 'LineWidth', 1, 'Color', "black"); 
plot([imm_delay_2, imm_delay_2], z, '--', 'LineWidth', 1, 'Color', "green");

xlim(x_box)
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlabel("Time (s) ")
labels = {'\color{black} Rotor Fail', ...
          '\color{black} IMM Delay', ...
          '\color{black} Reference'};
legend([h1, h2, h3], labels, 'FontSize', 8, 'FontName', 'cmr12',...
    'TextColor', 'black', "Position", [0.8, 0.9, 0, 0]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print(gcf,'twomode.png','-dpng','-r500')

