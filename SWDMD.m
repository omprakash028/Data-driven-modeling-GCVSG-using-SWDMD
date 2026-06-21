clc; clear; close all;

% ── Start total timer ──────────────────────────────────────────────────
t_total_start = tic;

%  OUTPUT FOLDER
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end
img_dir = fullfile(script_dir, 'Images');
if ~exist(img_dir,'dir')
    [ok, msg] = mkdir(img_dir);
    if ~ok, error('Cannot create Images folder:\n%s\n%s', img_dir, msg); end
end

%  GLOBAL STYLE
FS = 12;
FN = 'Times New Roman';
set(groot,'DefaultAxesFontSize',        FS);
set(groot,'DefaultAxesFontName',        FN);
set(groot,'DefaultAxesFontWeight',      'bold');
set(groot,'DefaultTextFontSize',        FS);
set(groot,'DefaultTextFontName',        FN);
set(groot,'DefaultTextFontWeight',      'bold');
set(groot,'DefaultLegendFontSize',      FS);
set(groot,'DefaultLegendFontName',      FN);
set(groot,'DefaultAxesLineWidth',       1.0);
set(groot,'DefaultLineLineWidth',       2.0);
set(groot,'DefaultAxesGridColor',       [0.5 0.5 0.5]);
set(groot,'DefaultAxesGridAlpha',       0.7);
set(groot,'DefaultAxesGridLineStyle',   '--');
set(groot,'DefaultFigureColor',         'w');
set(groot,'DefaultAxesBox',             'on');

%  STATE NAME DEFINITIONS
ylabels = {
    'v_{od} (V)',          'v_{oq} (V)', ...
    'i_{ld} (A)',          'i_{lq} (A)', ...
    'i_{od} (A)',          'i_{oq} (A)', ...
    'P_{inv} (W)',         'Q_{inv} (VAR)', ...
    '\Delta\omega_{VSM} (rad/s)'};
short_tex = {
    'V_{od}','V_{oq}','I_{ld}','I_{lq}', ...
    'I_{od}','I_{oq}','P_{inv}','Q_{inv}','\Delta\omega'};
short_csv = {'Vod','Voq','Ild','Ilq','Iod','Ioq','Pinv','Qinv','dw'};
eig_state_labels_full = {'V_{od}','V_{oq}','I_{Ld}','I_{Lq}', ...
                         'I_{od}','I_{oq}','P_{inv}','Q_{inv}', ...
                         '\Delta\omega_{VSM}'};

%  1. LOAD DATA
t_load_start = tic;
filename  = 'Out_Data1.xlsx';
raw       = readtable(filename);
time_full = raw{:,1};
X_full    = raw{:,2:10};
if size(X_full,2) ~= 9
    error('Data must have exactly 9 state columns (cols 2-10).');
end
dt = mean(diff(time_full));
fprintf('Full dataset: %d samples,  dt = %.6f s,  T_end = %.2f s\n', ...
        numel(time_full), dt, time_full(end));
fprintf('Data load time: %.3f s\n', toc(t_load_start))

%  2. TRIM TO 10 SECONDS
t_end = 10.0;
mask  = time_full <= t_end;
time  = time_full(mask);
X_all = X_full(mask, :);
N     = numel(time);
fprintf('Trimmed to %.0f s: %d samples\n', t_end, N);

%  3. SW-DMD PARAMETERS
r = 9;
w = max(2*r + 1, 20);
if w >= N
    error('Window size %d too large for %d samples in 10 s window.', w, N);
end
n_windows = N - w;
t_sw      = time(w : N-1);
fprintf('SW-DMD:  r=%d,  w=%d,  n_windows=%d\n', r, w, n_windows);

%  4. SW-DMD CORE LOOP
t_dmd_start = tic;
X_mat      = X_all';
lambda_all = zeros(r, n_windows);
x_recon    = zeros(9, n_windows);

for k = 1:n_windows
    Xw  = X_mat(:, k : k+w-1);
    Xwp = X_mat(:, k+1 : k+w);

    [U,S,V]    = svd(Xw,'econ');
    Ur = U(:,1:r);  Sr = S(1:r,1:r);  Vr = V(:,1:r);

    Atilde     = Ur' * Xwp * Vr / Sr;
    [W_e, D_e] = eig(Atilde);
    lambda     = diag(D_e);

    Phi  = Xwp * Vr / Sr * W_e;
    b    = Phi \ Xw(:,end);

    x_recon(:,k)    = real(Phi * (b .* lambda));
    lambda_all(:,k) = lambda;
end
fprintf('SW-DMD core loop time: %.3f s\n', toc(t_dmd_start))

%  5. CONTINUOUS EIGENVALUES
t_eig_start = tic;
omega_all   = log(lambda_all) / dt;

omega_last  = [];
valid_mask  = false(r,1);
best_window = 0;

for win_idx = n_windows : -1 : 1
    lambda_candidate = lambda_all(:, win_idx);
    omega_candidate  = log(lambda_candidate) / dt;
    vmask = isfinite(real(omega_candidate)) & isfinite(imag(omega_candidate));
    if sum(vmask) > sum(valid_mask)
        valid_mask  = vmask;
        omega_last  = omega_candidate(vmask);
        best_window = win_idx;
    end
    if all(vmask)
        break
    end
end

if isempty(omega_last)
    warning('SW-DMD: no finite eigenvalues found in any window.');
else
    fprintf('Eigenvalues: %d finite modes from window %d / %d\n', ...
            numel(omega_last), best_window, n_windows);
end
fprintf('Eigenvalue computation time: %.3f s\n', toc(t_eig_start))

%  6. ALIGN ORIGINAL DATA ONTO t_sw GRID
X_orig_sw = zeros(9, n_windows);
for s = 1:9
    X_orig_sw(s,:) = interp1(time, X_all(:,s), t_sw, 'linear');
end

%  7. MAX-NORMALIZATION
norm_max = @(Xin) Xin ./ max(max(abs(Xin),[],2), 1e-10);
dn_orig  = norm_max(X_orig_sw);
dn_recon = norm_max(x_recon);

%  8. ERROR METRICS
t_metrics_start = tic;
NMSE  = mean((dn_orig - dn_recon).^2, 2);
NRMSE = sqrt(NMSE);
R2    = zeros(9,1);
for s = 1:9
    yt = dn_orig(s,:);  yp = dn_recon(s,:);
    ss_res = sum((yt - yp).^2);
    ss_tot = sum((yt - mean(yt)).^2);
    if ss_tot < eps,  R2(s) = 0;
    else,             R2(s) = 1 - ss_res/ss_tot;
    end
end
yt_g      = dn_orig(:);  yp_g = dn_recon(:);
R2_global = 1 - sum((yt_g-yp_g).^2)/sum((yt_g-mean(yt_g)).^2);
fprintf('Error metrics time: %.3f s\n', toc(t_metrics_start))

%  9. PRINT ERROR TABLE
fprintf('\n%s\n', repmat('=',1,66));
fprintf('  SW-DMD Error Metrics  (0 to %.0f s,  max-normalised)\n', t_end);
fprintf('%s\n', repmat('-',1,66));
fprintf('  %-10s   %12s   %12s   %10s\n','State','NMSE','NRMSE','R^2');
fprintf('%s\n', repmat('-',1,66));
for s = 1:9
    fprintf('  %-10s   %12.6f   %12.6f   %10.6f\n', ...
            short_csv{s}, NMSE(s), NRMSE(s), R2(s));
end
fprintf('%s\n', repmat('-',1,66));
fprintf('  %-10s   %12.6f   %12.6f   %10.6f  (global)\n', ...
        'OVERALL', mean(NMSE), mean(NRMSE), R2_global);
fprintf('%s\n\n', repmat('=',1,66));

%  10. COLOURS
c_orig  = [0.00  0.447 0.741];
c_dmd   = [0.85  0.325 0.098];
c_nmse  = [0.22  0.62  0.88];
c_nrmse = [0.10  0.44  0.70];
c_r2    = [0.47  0.67  0.19];

%  11. X-AXIS HELPER
x_pos       = 1:9;
xt_vals     = 0 : 0.5 : t_end;
apply_xgrid = @(ax) set(ax, ...
    'XTick',              xt_vals, ...
    'XLim',               [0, t_end], ...
    'XTickLabelRotation', 45, ...
    'FontSize',           FS, 'FontName', FN);

%  12. STATE COMPARISON FIGURES
t_plot_start = tic;
pairs      = [1 2; 3 4; 5 6; 7 8];
fig_fnames = {'swdmd_Vod_Voq','swdmd_Ild_Ilq', ...
               'swdmd_Iod_Ioq','swdmd_Pinv_Qinv'};
leg_str    = {'Detail Switching Model (Original)', 'SW-DMD Reconstruction'};

for p = 1:4
    fig_p = figure('Color','w','Position',[80 80 1000 620]);
    set(fig_p,'Units','centimeters','Position',[2 2 16 10]);
    for k = 1:2
        idx = pairs(p,k);
        ax  = subplot(2,1,k);
        hold(ax,'on');
        h1 = plot(ax, time, X_all(:,idx), ...
                  'Color',c_orig,'LineStyle','-','LineWidth',2.5, ...
                  'DisplayName',leg_str{1});
        h2 = plot(ax, t_sw, x_recon(idx,:)', ...
                  'Color',c_dmd,'LineStyle','--','LineWidth',2.5, ...
                  'DisplayName',leg_str{2});
        hold(ax,'off');
        apply_xgrid(ax);
        grid(ax,'on');
        set(ax,'Box','on','FontSize',FS,'FontName',FN);
        yl = ylabel(ax, ylabels{idx}, 'Interpreter','tex','FontSize',FS);
        set(yl,'FontName',FN,'FontWeight','bold');
        if k == 2
            xl = xlabel(ax,'Time (s)','FontSize',FS);
            set(xl,'FontName',FN,'FontWeight','bold');
        end
        legend(ax,[h1 h2], leg_str, ...
               'Location','northeast','FontSize',FS,'FontName',FN,'Box','on');
    end
    set(fig_p,'PaperUnits','centimeters','PaperSize',[16 10], ...
              'PaperPosition',[0 0 16 10]);
    exportgraphics(fig_p, fullfile(img_dir,[fig_fnames{p} '.pdf']), ...
                   'ContentType','vector','Resolution',300);
    fprintf('  Saved: %s.pdf\n', fig_fnames{p});
end

%  13. DELTA OMEGA VSM PLOT
fig_dw = figure('Color','w','Position',[150 150 1000 420]);
set(fig_dw,'Units','centimeters','Position',[2 2 16 7]);
axDW = axes(fig_dw);
hold(axDW,'on');
h1 = plot(axDW, time, X_all(:,9), ...
          'Color',c_orig,'LineStyle','-','LineWidth',2.5,'DisplayName',leg_str{1});
h2 = plot(axDW, t_sw, x_recon(9,:)', ...
          'Color',c_dmd,'LineStyle','--','LineWidth',2.5,'DisplayName',leg_str{2});
hold(axDW,'off');
apply_xgrid(axDW);
grid(axDW,'on');
set(axDW,'Box','on','FontSize',FS,'FontName',FN);
xl_dw = xlabel(axDW,'Time (s)','FontSize',FS);
set(xl_dw,'FontName',FN,'FontWeight','bold');
yl_dw = ylabel(axDW,'\Delta\omega_{VSM} (rad/s)','Interpreter','tex','FontSize',FS);
set(yl_dw,'FontName',FN,'FontWeight','bold');
legend(axDW,[h1 h2],leg_str,'Location','northeast','FontSize',FS,'FontName',FN,'Box','on');
set(fig_dw,'PaperUnits','centimeters','PaperSize',[16 7],'PaperPosition',[0 0 16 7]);
exportgraphics(fig_dw, fullfile(img_dir,'swdmd_dw.pdf'), ...
               'ContentType','vector','Resolution',300);
fprintf('  Saved: swdmd_dw.pdf\n');

%  14. COMBINED ERROR BAR CHART (original — unchanged)
fig_err = figure('Color','w','Position',[80 60 1350 430]);
set(fig_err,'Units','centimeters','Position',[2 2 16 6]);

ax_n = subplot(1,3,1);
b1 = bar(ax_n, x_pos, NMSE, 0.62);
b1.FaceColor = c_nmse;  b1.EdgeColor = 'none';
set(ax_n,'XTick',x_pos,'XTickLabel',short_tex,'TickLabelInterpreter','tex', ...
         'XTickLabelRotation',40,'FontSize',FS,'FontName',FN,'Box','off');
yl_tmp = ylabel(ax_n,'NMSE','FontSize',FS);
set(yl_tmp,'FontName',FN,'FontWeight','bold');
xl_tmp = xlabel(ax_n,'State Variable','FontSize',FS);
set(xl_tmp,'FontName',FN,'FontWeight','bold');
grid(ax_n,'on');  ylim(ax_n,[0, max(NMSE)*1.30]);
for s = 1:9
    text(ax_n,s,NMSE(s)+max(NMSE)*0.04,sprintf('%.4f',NMSE(s)), ...
         'HorizontalAlignment','center','FontSize',9,'FontName',FN,'Rotation',90);
end

ax_nr = subplot(1,3,2);
b2 = bar(ax_nr, x_pos, NRMSE, 0.62);
b2.FaceColor = c_nrmse;  b2.EdgeColor = 'none';
set(ax_nr,'XTick',x_pos,'XTickLabel',short_tex,'TickLabelInterpreter','tex', ...
          'XTickLabelRotation',40,'FontSize',FS,'FontName',FN,'Box','off');
yl_tmp = ylabel(ax_nr,'NRMSE','FontSize',FS);
set(yl_tmp,'FontName',FN,'FontWeight','bold');
xl_tmp = xlabel(ax_nr,'State Variable','FontSize',FS);
set(xl_tmp,'FontName',FN,'FontWeight','bold');
grid(ax_nr,'on');  ylim(ax_nr,[0, max(NRMSE)*1.30]);
for s = 1:9
    text(ax_nr,s,NRMSE(s)+max(NRMSE)*0.04,sprintf('%.4f',NRMSE(s)), ...
         'HorizontalAlignment','center','FontSize',9,'FontName',FN,'Rotation',90);
end

ax_r2 = subplot(1,3,3);
b3 = bar(ax_r2, x_pos, R2, 0.62);
b3.FaceColor = c_r2;  b3.EdgeColor = 'none';
set(ax_r2,'XTick',x_pos,'XTickLabel',short_tex,'TickLabelInterpreter','tex', ...
          'XTickLabelRotation',40,'FontSize',FS,'FontName',FN,'Box','off');
yl_tmp = ylabel(ax_r2,'R^2 Score','FontSize',FS);
set(yl_tmp,'FontName',FN,'FontWeight','bold');
xl_tmp = xlabel(ax_r2,'State Variable','FontSize',FS);
set(xl_tmp,'FontName',FN,'FontWeight','bold');
ylim(ax_r2,[min(0,min(R2)-0.06), 1.15]);
yline(ax_r2,1.0,'k--','LineWidth',1.2,'HandleVisibility','off');
grid(ax_r2,'on');
for s = 1:9
    text(ax_r2,s,R2(s)+0.03,sprintf('%.4f',R2(s)), ...
         'HorizontalAlignment','center','FontSize',9,'FontName',FN,'Rotation',90);
end

set(fig_err,'PaperUnits','centimeters','PaperSize',[16 6],'PaperPosition',[0 0 16 6]);
exportgraphics(fig_err, fullfile(img_dir,'swdmd_error_bars.pdf'), ...
               'ContentType','vector','Resolution',300);
fprintf('  Saved: swdmd_error_bars.pdf\n');

%  14b. NMSE — INDIVIDUAL FIGURE
fig_nmse = figure('Color','w');
set(fig_nmse,'Units','centimeters','Position',[2 2 14 9]);
ax_nmse = axes(fig_nmse);
b_nmse = bar(ax_nmse, x_pos, NMSE, 0.62);
b_nmse.FaceColor = c_nmse;  b_nmse.EdgeColor = 'none';
set(ax_nmse,'XTick',x_pos,'XTickLabel',short_tex,'TickLabelInterpreter','tex', ...
    'XTickLabelRotation',40,'FontSize',FS,'FontName',FN,'Box','off');
yl_tmp = ylabel(ax_nmse,'NMSE','FontSize',FS);
set(yl_tmp,'FontName',FN,'FontWeight','bold');
xl_tmp = xlabel(ax_nmse,'State Variable','FontSize',FS);
set(xl_tmp,'FontName',FN,'FontWeight','bold');
grid(ax_nmse,'on');
ylim(ax_nmse,[0, max(NMSE)*1.35]);
for s = 1:9
    text(ax_nmse, s, NMSE(s)+max(NMSE)*0.04, sprintf('%.4f',NMSE(s)), ...
        'HorizontalAlignment','center','FontSize',9,'FontName',FN,'Rotation',90);
end
set(fig_nmse,'PaperUnits','centimeters','PaperSize',[14 9],'PaperPosition',[0 0 14 9]);
exportgraphics(fig_nmse, fullfile(img_dir,'swdmd_nmse.pdf'), ...
    'ContentType','vector','Resolution',300);
fprintf('  Saved: swdmd_nmse.pdf\n');

%  14c. NRMSE — INDIVIDUAL FIGURE
fig_nrmse = figure('Color','w');
set(fig_nrmse,'Units','centimeters','Position',[2 2 14 9]);
ax_nrmse = axes(fig_nrmse);
b_nrmse = bar(ax_nrmse, x_pos, NRMSE, 0.62);
b_nrmse.FaceColor = c_nrmse;  b_nrmse.EdgeColor = 'none';
set(ax_nrmse,'XTick',x_pos,'XTickLabel',short_tex,'TickLabelInterpreter','tex', ...
    'XTickLabelRotation',40,'FontSize',FS,'FontName',FN,'Box','off');
yl_tmp = ylabel(ax_nrmse,'NRMSE','FontSize',FS);
set(yl_tmp,'FontName',FN,'FontWeight','bold');
xl_tmp = xlabel(ax_nrmse,'State Variable','FontSize',FS);
set(xl_tmp,'FontName',FN,'FontWeight','bold');
grid(ax_nrmse,'on');
ylim(ax_nrmse,[0, max(NRMSE)*1.35]);
for s = 1:9
    text(ax_nrmse, s, NRMSE(s)+max(NRMSE)*0.04, sprintf('%.4f',NRMSE(s)), ...
        'HorizontalAlignment','center','FontSize',9,'FontName',FN,'Rotation',90);
end
set(fig_nrmse,'PaperUnits','centimeters','PaperSize',[14 9],'PaperPosition',[0 0 14 9]);
exportgraphics(fig_nrmse, fullfile(img_dir,'swdmd_nrmse.pdf'), ...
    'ContentType','vector','Resolution',300);
fprintf('  Saved: swdmd_nrmse.pdf\n');

%  14d. R2 — INDIVIDUAL FIGURE
fig_r2_solo = figure('Color','w');
set(fig_r2_solo,'Units','centimeters','Position',[2 2 14 9]);
ax_r2_solo = axes(fig_r2_solo);
b_r2_solo = bar(ax_r2_solo, x_pos, R2, 0.62);
b_r2_solo.FaceColor = c_r2;  b_r2_solo.EdgeColor = 'none';
set(ax_r2_solo,'XTick',x_pos,'XTickLabel',short_tex,'TickLabelInterpreter','tex', ...
    'XTickLabelRotation',40,'FontSize',FS,'FontName',FN,'Box','off');
yl_tmp = ylabel(ax_r2_solo,'R^2 Score','FontSize',FS);
set(yl_tmp,'FontName',FN,'FontWeight','bold');
xl_tmp = xlabel(ax_r2_solo,'State Variable','FontSize',FS);
set(xl_tmp,'FontName',FN,'FontWeight','bold');
ylim(ax_r2_solo,[min(0,min(R2)-0.08), 1.18]);
yline(ax_r2_solo, 1.0,'k--','LineWidth',1.2,'HandleVisibility','off');
grid(ax_r2_solo,'on');
for s = 1:9
    text(ax_r2_solo, s, R2(s)+0.03, sprintf('%.4f',R2(s)), ...
        'HorizontalAlignment','center','FontSize',9,'FontName',FN,'Rotation',90);
end
set(fig_r2_solo,'PaperUnits','centimeters','PaperSize',[14 9],'PaperPosition',[0 0 14 9]);
exportgraphics(fig_r2_solo, fullfile(img_dir,'swdmd_r2_solo.pdf'), ...
    'ContentType','vector','Resolution',300);
fprintf('  Saved: swdmd_r2_solo.pdf\n');

%  15. R² COMPARISON FIGURE
fig_r2 = figure('Color','w','Position',[300 180 860 500]);
set(fig_r2,'Units','centimeters','Position',[2 2 16 9]);
axR2 = axes(fig_r2);

br = bar(axR2, x_pos, R2, 0.62);
br.FaceColor = 'flat';
cmap_r2 = parula(9);
for s = 1:9,  br.CData(s,:) = cmap_r2(s,:);  end

hold(axR2,'on');
hl = yline(axR2, R2_global, 'r--', 'LineWidth', 2.2, ...
           'DisplayName', sprintf('Global R^2 = %.4f', R2_global));
yline(axR2, 1.0, 'k:', 'LineWidth', 1.0, 'HandleVisibility','off');
h_dummy = plot(axR2, NaN, NaN, 's', ...
               'MarkerFaceColor',[0.3 0.6 0.9], ...
               'MarkerEdgeColor','none','MarkerSize',10);
hold(axR2,'off');

set(axR2,'XTick',x_pos,'XTickLabel',short_tex,'TickLabelInterpreter','tex', ...
         'XTickLabelRotation',40,'FontSize',FS,'FontName',FN,'Box','off');
ylim(axR2,[min(0,min(R2)-0.08), 1.18]);
grid(axR2,'on');
yl_tmp = ylabel(axR2,'R^2 Score','FontSize',FS);
set(yl_tmp,'FontName',FN,'FontWeight','bold');
xl_tmp = xlabel(axR2,'State Variable','FontSize',FS);
set(xl_tmp,'FontName',FN,'FontWeight','bold');
legend(axR2, [h_dummy, hl], ...
       {'SW-DMD vs Detail Switching Model', ...
        sprintf('Global R^2 = %.4f', R2_global)}, ...
       'Location','south','FontSize',FS,'FontName',FN,'Box','on');

for s = 1:9
    text(axR2,s,R2(s)+0.025,sprintf('%.4f',R2(s)), ...
         'HorizontalAlignment','center','FontSize',11,'FontName',FN,'FontWeight','bold');
end
set(fig_r2,'PaperUnits','centimeters','PaperSize',[16 9],'PaperPosition',[0 0 16 9]);
exportgraphics(fig_r2, fullfile(img_dir,'swdmd_r2_comparison.pdf'), ...
               'ContentType','vector','Resolution',300);
fprintf('  Saved: swdmd_r2_comparison.pdf\n');

%  16. EXPORT CSV
writetable(table(short_csv', NMSE, NRMSE, R2, ...
    'VariableNames',{'State','NMSE','NRMSE','R2'}), ...
    fullfile(img_dir,'swdmd_error_table.csv'));
fprintf('  Saved: swdmd_error_table.csv\n');
fprintf('Figure plotting & export time: %.3f s\n', toc(t_plot_start))

%  17. EIGENVALUE MAP
t_eigmap_start  = tic;
eig_labels_plot = eig_state_labels_full(valid_mask);
n_valid         = numel(omega_last);
colors_plot     = lines(max(n_valid,1));

fig_eig = figure('Color','w','Units','centimeters','Position',[2 2 20 16]);
axE = axes(fig_eig);
hold(axE,'on');
xline(axE, 0,'k--','LineWidth',1.5,'HandleVisibility','off');

if n_valid > 0
    for k = 1:n_valid
        plot(axE, real(omega_last(k)), imag(omega_last(k)), 'x', ...
             'Color',      colors_plot(k,:), ...
             'MarkerSize', 14, ...
             'LineWidth',  2.5, ...
             'DisplayName', eig_labels_plot{k});
    end
    for k = 1:n_valid
        text(axE, real(omega_last(k)), imag(omega_last(k)), ...
             sprintf('  %s', eig_labels_plot{k}), ...
             'Interpreter','tex','FontSize',9,'FontName',FN,'Color',colors_plot(k,:));
    end
else
    text(0.5,0.5,'No finite eigenvalues found', ...
         'Units','normalized','HorizontalAlignment','center','FontSize',12);
end
hold(axE,'off');
drawnow;

if n_valid > 0
    re_all = real(omega_last);
    im_all = imag(omega_last);
    re_pad = max(50,  0.20*(max(re_all)-min(re_all)+eps));
    im_pad = max(500, 0.20*(max(im_all)-min(im_all)+eps));
    xlim(axE,[min(re_all)-re_pad,  max(re_all)+re_pad]);
    ylim(axE,[min(im_all)-im_pad,  max(im_all)+im_pad]);
end

xl_e = xlabel(axE,'Real Part (s^{-1})','FontSize',FS);
set(xl_e,'FontName',FN,'FontWeight','bold');
yl_e = ylabel(axE,'Imaginary Part (rad/s)','FontSize',FS);
set(yl_e,'FontName',FN,'FontWeight','bold');
if n_valid > 0
    legend(axE,'Location','bestoutside','FontSize',FS,'FontName',FN);
end
grid(axE,'on');  box(axE,'on');
set(axE,'FontSize',FS,'FontName',FN);

set(fig_eig,'PaperUnits','centimeters','PaperSize',[20 16],'PaperPosition',[0 0 20 16]);
exportgraphics(fig_eig, fullfile(img_dir,'swdmd_eigenvalue_map.pdf'), ...
               'ContentType','vector','Resolution',300);
fprintf('  Saved: swdmd_eigenvalue_map.pdf\n');
fprintf('Eigenvalue map time: %.3f s\n', toc(t_eigmap_start))

%  18. EIGENVALUE TABLE
fprintf('\n%s\n', repmat('=',1,66));
fprintf('  SW-DMD Eigenvalue Summary  (window %d / %d)\n', best_window, n_windows);
fprintf('  %-4s  %-22s  %10s  %12s  %8s  %s\n', ...
        '#','Mode','Real','Imag','|omega|','Stable?');
fprintf('%s\n', repmat('-',1,66));
for k = 1:n_valid
    flag = 'YES';
    if real(omega_last(k)) >= 0, flag = 'NO '; end
    fprintf('  %2d   %-22s  %8.4f  %+10.4fj  %6.4f    %s\n', ...
            k, eig_labels_plot{k}, ...
            real(omega_last(k)), imag(omega_last(k)), ...
            abs(omega_last(k)), flag);
end
fprintf('%s\n\n', repmat('=',1,66));
fprintf('All outputs saved to:  %s\n', img_dir);

% ── Final Timing Summary ───────────────────────────────────────────────
fprintf('\n========== Execution Time Summary ==========\n')
fprintf('  Data load              : %.3f s\n', toc(t_load_start))
fprintf('  SW-DMD core loop       : %.3f s\n', toc(t_dmd_start))
fprintf('  Eigenvalue computation : %.3f s\n', toc(t_eig_start))
fprintf('  Error metrics          : %.3f s\n', toc(t_metrics_start))
fprintf('  Figure plotting/export : %.3f s\n', toc(t_plot_start))
fprintf('  Eigenvalue map         : %.3f s\n', toc(t_eigmap_start))
fprintf('  TOTAL                  : %.3f s\n', toc(t_total_start))
fprintf('=============================================\n')