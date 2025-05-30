
%% WEAK MEAS ONLY, WITH RESPECT TO TIME

path = "/Users/sudirekci/Documents/MATLAB/weak_measurements/data/";

load(path+ "fisher_info_no_strong.mat");

load(path+ "mle_variance_no_strong_5_atoms.mat");
tau = tau_mle;
g = g_mle;
mle_variance_5 = mle_variance;
T_arr_5 = T_arr_mle;

load(path+ "mle_variance_no_strong_20_atoms.mat");
mle_variance_20 = mle_variance;
T_arr_20 = T_arr_mle;

load(path+ "mle_variance_no_strong_50_atoms.mat");
mle_variance_50 = mle_variance;
T_arr_50 = T_arr_mle;

% s = struct("tau_mle", tau, "g_mle" ,g, "T_arr_mle", T_arr_20, ...
%   "no_atoms", 20, "mle_variance", mle_variance_20);
% save(path+ "mle_variance_no_strong_20_atoms.mat", '-struct', 's');


figure
hold on
plot(g_fisher^2*T_arr_fisher/tau_fisher, log(total_fisher_info), Color=[10,10,255]/256, LineWidth=3)
plot(g^2*T_arr_5/tau, log(1./(mle_variance_5.*5)), LineWidth=3)
plot(g^2*T_arr_20/tau, log(1./(mle_variance_20.*20)), LineWidth=3)
plot(g^2*T_arr_50/tau, log(1./(mle_variance_50.*50)), LineWidth=3)
legend("Fisher information", "N = 5", "N = 20", "N = 50", "Location","northwest")
set(gca, 'XScale', 'log')
ax = gca;
ax.FontSize = 14; 
ax.FontName = "times";
xlabel('g^2 T/tau',FontName="times")
ylabel('ln(I/N)',FontName="times")



%% WITH STRONG MEAS, WITH RESPECT TO TIME

path = "/Users/sudirekci/Documents/MATLAB/weak_measurements/data/";

load(path+ "fisher_info_with_strong.mat");

load(path+ "mle_variance_with_strong_5_atoms.mat");
tau = tau_mle;
g = g_mle;
mle_variance_5 = mle_variance;
T_arr_5 = T_arr_mle;

load(path+ "mle_variance_with_strong_20_atoms.mat");
mle_variance_20 = mle_variance;
T_arr_20 = T_arr_mle;

load(path+ "mle_variance_with_strong_50_atoms.mat");
mle_variance_50 = mle_variance;
T_arr_50 = T_arr_mle;

% s = struct("tau_mle", tau, "g_mle", g, "T_arr_mle", T_arr_50, ...
%   "no_atoms", 50, "mle_variance", mle_variance_50);
% save(path+ "mle_variance_with_strong_50_atoms.mat", '-struct', 's');


figure
hold on
plot(T_arr_fisher, log(total_fisher_info), Color=[10,10,255]/256, LineWidth=3)
plot(T_arr_5, log(1./(mle_variance_5.*5)), LineWidth=3)
plot(T_arr_20, log(1./(mle_variance_20.*20)), LineWidth=3)
plot(T_arr_50, log(1./(mle_variance_50.*50)), LineWidth=3)
legend("Fisher information", "N = 5", "N = 20", "N = 50", "Location","northwest")
set(gca, 'XScale', 'log')
ax = gca;
ax.FontSize = 14; 
ax.FontName = "times";
xlabel('g^2 T/tau',FontName="times")
ylabel('ln(I/N)',FontName="times")



%% WEAK MEAS ONLY, WITH RESPECT TO MEASUREMENT STRENGTH

path = "/Users/sudirekci/Documents/MATLAB/weak_measurements/data/";

load(path+ "fisher_info_wrt_g_no_strong.mat");
g_fisher = g_arr_fisher;

load(path+ "mle_variance_wrt_g_no_strong_5_atoms.mat");
tau = tau_mle;
T = T_mle;
mle_variance_5 = mle_variance_mle;
g_arr_5 = g_mle;

load(path+ "mle_variance_wrt_g_no_strong_20_atoms.mat");
mle_variance_20 = mle_variance_mle;
g_arr_20 = g_mle;

load(path+ "mle_variance_wrt_g_no_strong_80_atoms.mat");
mle_variance_80 = mle_variance_mle;
g_arr_80 = g_mle;

% s = struct("tau_mle", tau, "g_mle" , g_arr_20, "T_mle", T, ... 
% "no_atoms", 20, "mle_variance_mle", mle_variance_20); 
% save(path+ "mle_variance_wrt_g_no_strong_20_atoms.mat", '-struct', 's');


figure
hold on
plot(g_fisher, total_fisher_info/(4*T^2), Color=[10,10,255]/256, LineWidth=3)
plot(g_arr_5, 1./(mle_variance_5.*5)/(4*T^2), '-o', LineWidth=3)
plot(g_arr_20, 1./(mle_variance_20.*20)/(4*T^2), '-o', LineWidth=3)
plot(g_arr_80, 1./(mle_variance_80.*80)/(4*T^2), '-o', LineWidth=3)
legend("Fisher information", "N = 5", "N = 20", "N = 80", "Location","northwest")
set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')
ax = gca;
ax.FontSize = 14; 
ax.FontName = "times";
xlabel('g^2 T/tau',FontName="times")
ylabel('I/4NT^2',FontName="times")



%% WITH STRONG MEAS, WITH RESPECT TO MEASUREMENT STRENGTH

path = "/Users/sudirekci/Documents/MATLAB/weak_measurements/data/";

load(path+ "fisher_info_wrt_g_with_strong.mat");
g_fisher = g_arr_fisher;

load(path+ "mle_variance_wrt_g_with_strong_5_atoms.mat");
tau = tau_mle;
T = T_mle;
mle_variance_5 = mle_variance_mle;
g_arr_5 = g_mle;

load(path+ "mle_variance_wrt_g_with_strong_20_atoms.mat");
mle_variance_20 = mle_variance_mle;
g_arr_20 = g_mle;

load(path+ "mle_variance_wrt_g_with_strong_80_atoms.mat");
mle_variance_80 = mle_variance_mle;
g_arr_80 = g_mle;

% s = struct("tau_mle", tau, "g_mle" , g_arr_80, "T_mle", T, ...
%   "no_atoms", 80, "mle_variance_mle", mle_variance_80);
% save(path+ "mle_variance_wrt_g_with_strong_80_atoms.mat", '-struct', 's');


figure
hold on
plot(g_fisher, total_fisher_info/(4*T^2), Color=[10,10,255]/256, LineWidth=3)
plot(g_arr_5, 1./(mle_variance_5.*5)/(4*T^2), '-o', LineWidth=3)
plot(g_arr_20, 1./(mle_variance_20.*20)/(4*T^2), '-o', LineWidth=3)
plot(g_arr_80, 1./(mle_variance_80.*80)/(4*T^2), '-o', LineWidth=3)
legend("Fisher information", "N = 5", "N = 20", "N = 80", "Location","northeast")
set(gca, 'XScale', 'log')
ax = gca;
ax.FontSize = 14; 
ax.FontName = "times";
xlabel('g^2 T/tau',FontName="times")
ylabel('I/4NT^2',FontName="times")


%%

path = "/Users/sudirekci/Documents/MATLAB/weak_measurements/data/";

load(path + "rosenband_freq_optimal_ensemble.mat")
load(path + "delta_omega_our_protocol.mat");

omega_corner = pi; % from -pi to pi
delta_omega = omega_corner/sqrt(12);
N = 64;

% s = struct("no_atoms", N, "T_arr", T_arr, ...
%     "bmse", bmse, "optimal_gs", optimal_gs);
% save(path + "delta_omega_our_protocol.mat", '-struct', 's');

figure
hold on
plot(delta_omega*T_arr, 1./(T_arr), LineWidth=3)
plot(delta_omega*rosenband_T_arr, sqrt(rosenband_bmse)*2*sqrt(N), LineWidth=3)
plot(delta_omega*T_arr, sqrt(bmse)*2*sqrt(N), '-o', LineWidth=3)
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
ax = gca;
ax.FontSize = 20; 
ax.FontName = "times";
xlabel("\delta\omega T")
lab = ylabel('$$\sqrt N \Delta\omega$$');
set(lab,'Interpreter','latex') 
legend("1/T", "Rosenband", "Our protocol")