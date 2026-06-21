clear; clc; close all;

% ── Start total timer ──────────────────────────────────────────────────
t_total_start = tic;                                                       % <-- ADDED

% PARAMETERS
wn = 2*pi*50;
wf = 4*pi;
Rf = 0.07; Lf = 5.2e-3; Cf = 100e-6;
Rc = 0.05; Lc = 1.5e-4;
Pref = 50000;
Qref = 0;
mq = 1e-5;
Kpv = 3; Kiv = 400;
Kpc = 5; Kic = 10;
J = 0.05; D = 10;
Vdref = sqrt(2/3)*400;
Vqref = 0;
vbD = sqrt(2/3)*400;
vbQ = 0;

% STEP 1: Run nonlinear ODE to find true steady state
x0 = [Vdref; 0; 0; 0; 0; 0; Pref; 0; 0; 0; 0; 0; 0; 0; 0];
tspan_nl = [0 30];
options_nl = odeset('RelTol',1e-6,'AbsTol',1e-9,'MaxStep',5e-6);
t_ode_start = tic;                                                         % <-- ADDED
[t_nl, X_nl] = ode45(@(t,x) inverter_dynamics(x,Pref,Qref), tspan_nl, x0, options_nl);
fprintf('ODE solve time (steady-state): %.3f s\n', toc(t_ode_start))      % <-- ADDED

% Use final state as true operating point
x_eq = X_nl(end, :)';

fprintf('True steady-state operating point:\n')
fprintf('  Vod = %.4f V\n',  x_eq(1))
fprintf('  Voq = %.4f V\n',  x_eq(2))
fprintf('  Ild = %.4f A\n',  x_eq(3))
fprintf('  Ilq = %.4f A\n',  x_eq(4))
fprintf('  Iod = %.4f A\n',  x_eq(5))
fprintf('  Ioq = %.4f A\n',  x_eq(6))
fprintf('  Pinv = %.4f W\n', x_eq(7))
fprintf('  Qinv = %.4f VAR\n',x_eq(8))
fprintf('  dw   = %.6f rad/s\n', x_eq(9))

n = 15;
m = 2;

% STEP 2: Linearize around true steady state
t_lin_start = tic;                                                         % <-- ADDED
f = @(x,Pref,Qref) inverter_dynamics(x,Pref,Qref);
A = zeros(n,n);
B = zeros(n,m);
eps_fd = 1e-6;

for k = 1:n
    x_p = x_eq; x_p(k) = x_p(k) + eps_fd;
    x_m = x_eq; x_m(k) = x_m(k) - eps_fd;
    A(:,k) = (f(x_p,Pref,Qref) - f(x_m,Pref,Qref))/(2*eps_fd);
end

for k = 1:m
    if k == 1
        B(:,k) = (f(x_eq,Pref+eps_fd,Qref) - f(x_eq,Pref-eps_fd,Qref))/(2*eps_fd);
    else
        B(:,k) = (f(x_eq,Pref,Qref+eps_fd) - f(x_eq,Pref,Qref-eps_fd))/(2*eps_fd);
    end
end
fprintf('Linearization time: %.3f s\n', toc(t_lin_start))                 % <-- ADDED

% STEP 3: Simulate linearized model
t_linsim_start = tic;                                                      % <-- ADDED
u_step = [1000; 0];
ode_lin = @(t, x) A*(x - x_eq) + B*u_step;

tspan   = [0 10];
options = odeset('RelTol',1e-6,'AbsTol',1e-9,'MaxStep',1e-3);
[t_lin, x_lin] = ode45(ode_lin, tspan, x_eq, options);
fprintf('Linearized ODE simulation time: %.3f s\n', toc(t_linsim_start))  % <-- ADDED

% -----------------------------------------------------------------------
% Figure 1: Vod & Voq
% -----------------------------------------------------------------------
t_plot_start = tic;                                                        % <-- ADDED
figure('Name','Output Voltages','NumberTitle','off', ...
       'Units','normalized','Position',[0.05 0.55 0.42 0.35])

subplot(2,1,1)
plot(t_lin, x_lin(:,1), 'Color',[0.00 0.45 0.74], 'LineWidth',1.6, ...
    'DisplayName','Output voltage d-axis, V_{od}')
ylabel('V_{od} (V)', 'FontSize',9)
legend('Location','best','FontSize',8)
grid on

subplot(2,1,2)
plot(t_lin, x_lin(:,2), 'Color',[0.85 0.33 0.10], 'LineWidth',1.6, ...
    'DisplayName','Output voltage q-axis, V_{oq}')
xlabel('Time (s)', 'FontSize',9)
ylabel('V_{oq} (V)', 'FontSize',9)
legend('Location','best','FontSize',8)
grid on

% -----------------------------------------------------------------------
% Figure 2: Iod & Ioq
% -----------------------------------------------------------------------
figure('Name','Output Currents','NumberTitle','off', ...
       'Units','normalized','Position',[0.52 0.55 0.42 0.35])

subplot(2,1,1)
plot(t_lin, x_lin(:,5), 'Color',[0.93 0.69 0.13], 'LineWidth',1.6, ...
    'DisplayName','Output current d-axis, I_{od}')
ylabel('I_{od} (A)', 'FontSize',9)
legend('Location','best','FontSize',8)
grid on

subplot(2,1,2)
plot(t_lin, x_lin(:,6), 'Color',[0.30 0.75 0.93], 'LineWidth',1.6, ...
    'DisplayName','Output current q-axis, I_{oq}')
xlabel('Time (s)', 'FontSize',9)
ylabel('I_{oq} (A)', 'FontSize',9)
legend('Location','best','FontSize',8)
grid on

% -----------------------------------------------------------------------
% Figure 3: Pinv & Qinv
% -----------------------------------------------------------------------
figure('Name','Active and Reactive Power','NumberTitle','off', ...
       'Units','normalized','Position',[0.05 0.10 0.42 0.35])

subplot(2,1,1)
plot(t_lin, x_lin(:,7), 'Color',[0.64 0.08 0.18], 'LineWidth',1.6, ...
    'DisplayName','Active power, P_{inv}')
ylabel('P_{inv} (W)', 'FontSize',9)
legend('Location','best','FontSize',8)
grid on

subplot(2,1,2)
plot(t_lin, x_lin(:,8), 'Color',[0.00 0.60 0.00], 'LineWidth',1.6, ...
    'DisplayName','Reactive power, Q_{inv}')
xlabel('Time (s)', 'FontSize',9)
ylabel('Q_{inv} (VAR)', 'FontSize',9)
legend('Location','best','FontSize',8)
grid on

% -----------------------------------------------------------------------
% Figure 4: Delta-omega VSM
% -----------------------------------------------------------------------
figure('Name','VSM Frequency Deviation','NumberTitle','off', ...
       'Units','normalized','Position',[0.52 0.10 0.42 0.35])

plot(t_lin, x_lin(:,9), 'Color',[0.10 0.10 0.80], 'LineWidth',1.6, ...
    'DisplayName','VSM frequency deviation, \Delta\omega')
xlabel('Time (s)', 'FontSize',9)
ylabel('\Delta\omega (rad/s)', 'FontSize',9)
legend('Location','best','FontSize',8)
grid on

% -----------------------------------------------------------------------
% Figure 5: Eigenvalue Plot
% -----------------------------------------------------------------------
t_eig_start = tic;                                                         % <-- ADDED
eig_vals = eig(A);

colors = [
    0.00  0.45  0.74;
    0.85  0.33  0.10;
    0.93  0.69  0.13;
    0.49  0.18  0.56;
    0.47  0.67  0.19;
    0.30  0.75  0.93;
    0.64  0.08  0.18;
    1.00  0.60  0.00;
    0.00  0.60  0.00;
    0.60  0.00  0.60;
    0.00  0.45  0.45;
    0.80  0.80  0.00;
    0.10  0.10  0.80;
    0.80  0.40  0.00;
    0.40  0.40  0.40;
];

state_labels = {
    '\lambda_1  –  V_{od}',
    '\lambda_2  –  V_{oq}',
    '\lambda_3  –  I_{Ld}',
    '\lambda_4  –  I_{Lq}',
    '\lambda_5  –  I_{od}',
    '\lambda_6  –  I_{oq}',
    '\lambda_7  –  P_{inv}',
    '\lambda_8  –  Q_{inv}',
    '\lambda_9  –  \Delta\omega',
    '\lambda_{10} –  \alpha_{VSM}',
    '\lambda_{11} –  \delta_{VSM}',
    '\lambda_{12} –  \gamma_d',
    '\lambda_{13} –  \gamma_q',
    '\lambda_{14} –  \zeta_d',
    '\lambda_{15} –  \zeta_q'};

figure('Name','Eigenvalue Plot','NumberTitle','off', ...
       'Units','normalized','Position',[0.15 0.10 0.70 0.75])
hold on
for i = 1:n
    plot(real(eig_vals(i)), imag(eig_vals(i)), 'x', ...
         'Color', colors(i,:), 'MarkerSize', 14, 'LineWidth', 2.5, ...
         'DisplayName', state_labels{i})
end
xline(0, 'k--', 'LineWidth', 1.2, 'HandleVisibility','off')
xlabel('Real Part  (s^{-1})',     'FontSize', 11)
ylabel('Imaginary Part  (rad/s)', 'FontSize', 11)
legend('Location', 'bestoutside', 'FontSize', 9)
grid on
hold off
fprintf('Eigenvalue analysis time: %.3f s\n', toc(t_eig_start))           % <-- ADDED
fprintf('Figure plotting time: %.3f s\n', toc(t_plot_start))              % <-- ADDED

% PRINT TABLE
fprintf('\n===== Eigenvalue Table =====\n')
for i = 1:n
    fprintf('%2d (%s): %10.4f + j%10.4f\n', ...
        i, state_labels{i}, real(eig_vals(i)), imag(eig_vals(i)));
end

% ── Final Timing Summary ───────────────────────────────────────────────
fprintf('\n========== Execution Time Summary ==========\n')               % <-- ADDED
fprintf('  ODE solve (steady-state) : %.3f s\n', toc(t_ode_start))       % <-- ADDED
fprintf('  Linearization            : %.3f s\n', toc(t_lin_start))       % <-- ADDED
fprintf('  Linearized ODE sim       : %.3f s\n', toc(t_linsim_start))    % <-- ADDED
fprintf('  Eigenvalue analysis      : %.3f s\n', toc(t_eig_start))       % <-- ADDED
fprintf('  Figure plotting          : %.3f s\n', toc(t_plot_start))      % <-- ADDED
fprintf('  TOTAL                    : %.3f s\n', toc(t_total_start))     % <-- ADDED
fprintf('=============================================\n')                 % <-- ADDED


% ── DYNAMICS FUNCTION (unchanged) ────────────────────────────────────────────
function dx = inverter_dynamics(x,Pref,Qref)
    Rf = 0.07; Lf = 5.2e-3; Cf = 100e-6;
    Rc = 0.05; Lc = 1.5e-4;
    Kpv = 3; Kiv = 400;
    Kpc = 5; Kic = 10;
    J = 0.05; D = 10;
    mq = 1e-5;
    wn = 2*pi*50;
    wf = 4*pi;
    Vdref = sqrt(2/3)*400;
    Vqref = 0;
    vbD = sqrt(2/3)*400;
    vbQ = 0;

    vod=x(1); voq=x(2); ild=x(3); ilq=x(4);
    iod=x(5); ioq=x(6); Pinv=x(7); Qinv=x(8);
    dw=x(9); a=x(10); delta=x(11);
    gd=x(12); gq=x(13); zd=x(14); zq=x(15);

    vbd = cos(delta)*vbD + sin(delta)*vbQ;
    vbq = -sin(delta)*vbD + cos(delta)*vbQ;

    pinv = (3/2)*(vod*iod + voq*ioq);
    qinv = (3/2)*(voq*iod - vod*ioq);

    Pinv_dot = -wf*Pinv + wf*pinv;
    Qinv_dot = -wf*Qinv + wf*qinv;

    dw_dot = ((Pref - Pinv)/wn - D*dw)/J;
    w = wn + dw;
    a_dot = w;
    delta_dot = dw;

    vodref = Vdref + (Qref - Qinv)*mq - ild*Rf + ilq*w*Lf;
    voqref = Vqref - ilq*Rf - ild*w*Lf;

    gd_dot = vodref - vod;
    gq_dot = voqref - voq;
    ifd = Kpv*(vodref - vod) + Kiv*gd + iod - w*Cf*voq;
    ifq = Kpv*(voqref - voq) + Kiv*gq + ioq + w*Cf*vod;

    zd_dot = ifd - ild;
    zq_dot = ifq - ilq;
    vfd = Kpc*(ifd - ild) + Kic*zd + vod - w*Lf*ilq;
    vfq = Kpc*(ifq - ilq) + Kic*zq + voq + w*Lf*ild;

    ild_dot = (-Rf/Lf)*ild + w*ilq + (vfd - vod)/Lf;
    ilq_dot = -w*ild + (-Rf/Lf)*ilq + (vfq - voq)/Lf;
    vod_dot = w*voq + (ild - iod)/Cf;
    voq_dot = -w*vod + (ilq - ioq)/Cf;

    iod_dot = (-Rc/Lc)*iod + w*ioq + (vod - vbd)/Lc;
    ioq_dot = -w*iod + (-Rc/Lc)*ioq + (voq - vbq)/Lc;

    dx = [vod_dot; voq_dot; ild_dot; ilq_dot; iod_dot; ioq_dot;
          Pinv_dot; Qinv_dot; dw_dot; a_dot; delta_dot;
          gd_dot; gq_dot; zd_dot; zq_dot];
end