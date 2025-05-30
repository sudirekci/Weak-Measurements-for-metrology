
close all

path = "/Users/sudirekci/Documents/MATLAB/weak_measurements/data/";

omega_corner = pi; % from -pi to pi
delta_omega = omega_corner/sqrt(12);

L = 6;

T_arr = 10.^(log10(0.1):0.02:log10(10.));
rosenband_arr = zeros([L, length(T_arr)]);

for i = 1:L

    no_ensembles = i;

    load(path + "rosenband_freq_" + num2str(no_ensembles) + "_ensemble.mat");

    rosenband_arr(i, :) = bmse;

end

rosenband_final = min(rosenband_arr, [], 1);

% s = struct("rosenband_T_arr", T_arr, "rosenband_bmse", rosenband_final);
% save(path + "rosenband_freq_optimal_ensemble.mat", '-struct', 's');

N = 64;

figure
hold on
plot(T_arr, 1./(T_arr), LineWidth=3)
plot(T_arr, sqrt(rosenband_final)*2*sqrt(N), LineWidth=3)
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
ax = gca;
ax.FontSize = 20; 
ax.FontName = "times";
xlabel("\delta\omega T")
lab = ylabel('$$\sqrt N \Delta\omega$$');
set(lab,'Interpreter','latex') 



figure
hold on
grid on
plot(delta_omega*T_arr, 1./(T_arr), LineWidth=3)
plot(delta_omega*T_arr, sqrt(rosenband_arr)*2*sqrt(N), LineWidth=3)
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
ax = gca;
ax.FontSize = 20; 
ax.FontName = "times";
xlabel("\delta\omega T")
lab = ylabel('$$\sqrt N \Delta\omega$$');
set(lab,'Interpreter','latex') 
h = legend('$$1/T$$', "M=1", "M=2", "M=3", "M=4", "M=5", "M=6");
set(h,'Interpreter','latex') 

% figure
% hold on
% plot(delta_omega*T_arr, 1.55./(T_arr), LineWidth=3)
% plot(delta_omega*T_arr, sqrt(rosenband_arr(2, :))*2*sqrt(N), LineWidth=3)
% set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')
% ax = gca;
% ax.FontSize = 20; 
% ax.FontName = "times";
% xlabel("\delta\omega T")
% lab = ylabel('$$\sqrt N \Delta\omega$$');
% set(lab,'Interpreter','latex') 