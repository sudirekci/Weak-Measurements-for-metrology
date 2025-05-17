
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

figure
hold on
plot(g_fisher^2*T_arr_fisher/tau_fisher, log(total_fisher_info), Color=[10,10,255]/256, LineWidth=3)
plot(g^2*T_arr_5/tau, log(1./(mle_variance_5.*5)), LineWidth=3)
plot(g^2*T_arr_5/tau, log(1./(mle_variance_20.*20)), LineWidth=3)
plot(g^2*T_arr_5/tau, log(1./(mle_variance_50.*50)), LineWidth=3)
legend("Fisher information", "N = 5", "N = 20", "N = 50")
set(gca, 'XScale', 'log')
ax = gca;
ax.FontSize = 14; 
ax.FontName = "times";
xlabel('g^2 T/tau',FontName="times")
ylabel('ln(I/N)',FontName="times")

