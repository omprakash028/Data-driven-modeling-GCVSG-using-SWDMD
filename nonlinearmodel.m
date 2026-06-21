clear all
clc

% ── Start total timer ──────────────────────────────────────────────────
t_total_start = tic;

% Parameters
Vdc  = 1000;
Vg   = 400;
Vref = sqrt(2/3)*400;
Vqref = 0;
f = 50;
wf = 4*pi;
fc = 10e3;
wn = 2*pi*f;
Rf = 0.07;
Lf = 5.2e-3;
Cf = 100e-6;
Rc   = 0.05;
Lc   = 1.5e-4;
Rv   = Rf;
Lv   = Lf;
Pref = 50000;
Qref = 0;
ka   = 3.14e-3;
mq   = 1e-5;
Kpv  = 3;
Kiv  = 400;
Kpc  = 5;
Kic  = 10;
J    = 0.05;
D    = 10;

% Pack parameters into struct
params.wf = wf;
params.wn = wn;
params.Rf = Rf;
params.Lf = Lf;
params.Cf = Cf;
params.Rc = Rc;
params.Lc = Lc;
params.Pref  = Pref;
params.Qref  = Qref;
params.mq = mq;
params.J = J;
params.D = D;
params.Kpv = Kpv;
params.Kiv = Kiv;
params.Kpc = Kpc;
params.Kic = Kic;
params.Vdref = sqrt(2/3)*400;
params.Vqref = 0;
params.vbD = sqrt(2/3)*400;
params.vbQ = 0;

% Initial conditions
x0 = [params.Vdref; 0; 0; 0; 0; 0; Pref; 0; 0; 0; 0; 0; 0; 0; 0];

% Solve ODE – 10 seconds only
t_ode_start = tic;                                                         % <-- ADDED
tspan   = [0 10];
options = odeset('RelTol',1e-6,'AbsTol',1e-9,'MaxStep',5e-6);
[t, X]  = ode45(@(t,x) inverter_model(t, x, params), tspan, x0, options);
fprintf('ODE solve time: %.3f s\n', toc(t_ode_start))                     % <-- ADDED

% -----------------------------------------------------------------------
% Figure 1: Output Voltages – Vod & Voq  (2 subplots stacked)
% -----------------------------------------------------------------------
t_plot_start = tic;                                                        % <-- ADDED
figure('Name','Output Voltages','NumberTitle','off', ...
       'Units','normalized','Position',[0.05 0.55 0.42 0.35])

subplot(2,1,1)
plot(t, X(:,1), 'Color',[0.00 0.45 0.74], 'LineWidth',1.8)
ylabel('V_{od} (V)', 'FontSize',10)
grid on;  xlim([0 10])

subplot(2,1,2)
plot(t, X(:,2), 'Color',[0.85 0.33 0.10], 'LineWidth',1.8)
xlabel('Time (s)', 'FontSize',10)
ylabel('V_{oq} (V)', 'FontSize',10)
grid on;  xlim([0 10])

% -----------------------------------------------------------------------
% Figure 2: Output Currents – Iod & Ioq  (2 subplots stacked)
% -----------------------------------------------------------------------
figure('Name','Output Currents','NumberTitle','off', ...
       'Units','normalized','Position',[0.52 0.55 0.42 0.35])

subplot(2,1,1)
plot(t, X(:,5), 'Color',[0.47 0.67 0.19], 'LineWidth',1.8)
ylabel('I_{od} (A)', 'FontSize',10)
grid on;  xlim([0 10])

subplot(2,1,2)
plot(t, X(:,6), 'Color',[0.49 0.18 0.56], 'LineWidth',1.8)
xlabel('Time (s)', 'FontSize',10)
ylabel('I_{oq} (A)', 'FontSize',10)
grid on;  xlim([0 10])

% -----------------------------------------------------------------------
% Figure 3: Active & Reactive Power – Pinv & Qinv  (2 subplots stacked)
% -----------------------------------------------------------------------
figure('Name','Active and Reactive Power','NumberTitle','off', ...
       'Units','normalized','Position',[0.05 0.10 0.42 0.35])

subplot(2,1,1)
plot(t, X(:,7), 'Color',[0.64 0.08 0.18], 'LineWidth',1.8)
ylabel('P_{inv} (W)', 'FontSize',10)
grid on;  xlim([0 10])

subplot(2,1,2)
plot(t, X(:,8), 'Color',[0.00 0.60 0.00], 'LineWidth',1.8)
xlabel('Time (s)', 'FontSize',10)
ylabel('Q_{inv} (VAR)', 'FontSize',10)
grid on;  xlim([0 10])

% -----------------------------------------------------------------------
% Figure 4: VSM Frequency Deviation – Delta-omega  (single plot)
% -----------------------------------------------------------------------
figure('Name','VSM Frequency Deviation','NumberTitle','off', ...
       'Units','normalized','Position',[0.52 0.10 0.42 0.35])

plot(t, X(:,9), 'Color',[0.10 0.10 0.80], 'LineWidth',1.8)
xlabel('Time (s)', 'FontSize',10)
ylabel('\Delta\omega_{VSM} (rad/s)', 'FontSize',10)
grid on;  xlim([0 10])
fprintf('Plotting time: %.3f s\n', toc(t_plot_start))                     % <-- ADDED

% -----------------------------------------------------------------------
% Eigenvalue Analysis (at t = 10 s)
% -----------------------------------------------------------------------
t_eig_start = tic;                                                         % <-- ADDED
x_eq = X(end, :)';
n_states = length(x_eq);
A = zeros(n_states, n_states);
eps_fd = 1e-6;

for k = 1:n_states
    x_p    = x_eq;  x_p(k) = x_p(k) + eps_fd;
    x_m    = x_eq;  x_m(k) = x_m(k) - eps_fd;
    A(:,k) = (inverter_model(t(end),x_p,params) - ...
              inverter_model(t(end),x_m,params)) / (2*eps_fd);
end

eigs_val = eig(A);

colors_eig = [
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
    '\lambda_1 – V_{od}',
    '\lambda_2 – V_{oq}',
    '\lambda_3 – I_{Ld}',
    '\lambda_4 – I_{Lq}',
    '\lambda_5 – I_{od}',
    '\lambda_6 – I_{oq}',
    '\lambda_7 – P_{inv}',
    '\lambda_8 – Q_{inv}',
    '\lambda_9 – \Delta\omega_{VSM}',
    '\lambda_{10} – \alpha_{VSM}',
    '\lambda_{11} – \delta_{VSM}',
    '\lambda_{12} – \gamma_d',
    '\lambda_{13} – \gamma_q',
    '\lambda_{14} – \zeta_d',
    '\lambda_{15} – \zeta_q'};

figure('Name','Eigenvalue Plot','NumberTitle','off', ...
       'Units','normalized','Position',[0.15 0.10 0.70 0.75])
hold on
for k = 1:n_states
    plot(real(eigs_val(k)), imag(eigs_val(k)), 'x', ...
         'Color', colors_eig(k,:), 'MarkerSize', 14, 'LineWidth', 2.5, ...
         'DisplayName', state_labels{k})
end
xline(0, 'k--', 'LineWidth', 1.2, 'HandleVisibility','off')
xlabel('Real Part  (s^{-1})', 'FontSize', 11)
ylabel('Imaginary Part  (rad/s)', 'FontSize', 11)
legend('Location', 'bestoutside', 'FontSize', 9)
grid on
hold off

% Eigenvalue summary table
fprintf('\n===== Eigenvalue Summary =====\n')
fprintf('  #      Real          Imag         |lambda|   Stable?\n')
fprintf('------------------------------------------------------------\n')
for k = 1:n_states
    flag = 'YES';
    if real(eigs_val(k)) >= 0, flag = 'NO '; end
    fprintf('  %2d  %12.4f  %12.4f  %10.4f    %s\n', ...
            k, real(eigs_val(k)), imag(eigs_val(k)), abs(eigs_val(k)), flag)
end
fprintf('Eigenvalue analysis time: %.3f s\n', toc(t_eig_start))          % <-- ADDED

% ── Final Timing Summary ───────────────────────────────────────────────% <-- ADDED
fprintf('\n========== Execution Time Summary ==========\n')               % <-- ADDED
fprintf('  ODE simulation   : %.3f s\n', toc(t_ode_start))               % <-- ADDED
fprintf('  Figure plotting  : %.3f s\n', toc(t_plot_start))              % <-- ADDED
fprintf('  Eigenvalue anal. : %.3f s\n', toc(t_eig_start))               % <-- ADDED
fprintf('  TOTAL            : %.3f s\n', toc(t_total_start))             % <-- ADDED
fprintf('=============================================\n')                 % <-- ADDED


% -----------------------------------------------------------------------
% ODE function
% -----------------------------------------------------------------------
function dxdt = inverter_model(t, x, p)

    vod      = x(1);
    voq      = x(2);
    ild      = x(3);
    ilq      = x(4);
    iod      = x(5);
    ioq      = x(6);
    Pinv     = x(7);
    Qinv     = x(8);
    dwVSM    = x(9);
    aVSM     = x(10);
    deltaVSM = x(11);
    gammad   = x(12);
    gammaq   = x(13);
    zetad    = x(14);
    zetaq    = x(15);

    % Grid voltage in inverter dq frame
    vbd = cos(deltaVSM)*p.vbD + sin(deltaVSM)*p.vbQ;
    vbq = -sin(deltaVSM)*p.vbD + cos(deltaVSM)*p.vbQ;

    % Instantaneous power
    pinv = (3/2)*(vod*iod + voq*ioq);
    qinv = (3/2)*(voq*iod - vod*ioq);

    % Low-pass filtered power
    Pinv_dot = -p.wf*Pinv + p.wf*pinv;
    Qinv_dot = -p.wf*Qinv + p.wf*qinv;

    % VSM frequency dynamics
    dwVSM_dot = ((p.Pref - Pinv)/p.wn - p.D*dwVSM) / p.J;
    wVSM = p.wn + dwVSM;
    aVSM_dot = wVSM;
    deltaVSM_dot = dwVSM;

    % Voltage reference with Q-droop
    vodrefin = p.Vdref + (p.Qref - Qinv)*p.mq;
    voqrefin = p.Vqref;

    % Virtual impedance
    vodref = vodrefin - ild*p.Rf + ilq*wVSM*p.Lf;
    voqref = voqrefin - ilq*p.Rf - ild*wVSM*p.Lf;

    % Voltage controller
    gammad_dot = vodref - vod;
    gammaq_dot = voqref - voq;
    ifdref = p.Kpv*(vodref - vod) + p.Kiv*gammad + iod - wVSM*p.Cf*voq;
    ifqref = p.Kpv*(voqref - voq) + p.Kiv*gammaq + ioq + wVSM*p.Cf*vod;

    % Current controller
    zetad_dot = ifdref - ild;
    zetaq_dot = ifqref - ilq;
    vfdref = p.Kpc*(ifdref - ild) + p.Kic*zetad + vod - wVSM*p.Lf*ilq;
    vfqref = p.Kpc*(ifqref - ilq) + p.Kic*zetaq + voq + wVSM*p.Lf*ild;

    % LC filter dynamics
    ild_dot = (-p.Rf/p.Lf)*ild + wVSM*ilq + (vfdref - vod)/p.Lf;
    ilq_dot = -wVSM*ild + (-p.Rf/p.Lf)*ilq + (vfqref - voq)/p.Lf;
    vod_dot =  wVSM*voq + (ild - iod)/p.Cf;
    voq_dot = -wVSM*vod + (ilq - ioq)/p.Cf;

    % Coupling line dynamics
    iod_dot = (-p.Rc/p.Lc)*iod + wVSM*ioq + (vod - vbd)/p.Lc;
    ioq_dot = -wVSM*iod + (-p.Rc/p.Lc)*ioq + (voq - vbq)/p.Lc;

    % Output derivative vector
    dxdt = [vod_dot; voq_dot; ild_dot; ilq_dot; iod_dot; ioq_dot; ...
            Pinv_dot; Qinv_dot; dwVSM_dot; aVSM_dot; deltaVSM_dot; ...
            gammad_dot; gammaq_dot; zetad_dot; zetaq_dot];
end