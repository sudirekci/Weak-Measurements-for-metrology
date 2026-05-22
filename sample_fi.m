
% CREATE TARGET PATH

targetPath = fullfile(pwd, 'data_folder');

if ~isfolder(targetPath)
    mkdir(targetPath);
end



%% FISHER INFO SAMPLING

%close all

T = 20;  % total time
tau = 0.1;  %weak meas time
w = 1; % frequency

no_atoms = 10;

no_trajectory_samples = 1000; % number of samples
g_array = 10.^(-2.:0.07:-0.2);

strong_meas = 0; % 1 for strong measurement at the end
classical = 0; % 1 for 

fi_fit = (8*g_array.^2.*no_atoms*T^3)/(3*tau)./(1+3/4*g_array.^2*T/tau+ ...
                                +2/3*g_array.^4*T^2/tau^2);

tic

[total_fisher_info] = avg_fisher_information_per_trajectory(g_array, ...
    no_trajectory_samples, T, tau, w, no_atoms, strong_meas, classical, ...
    0);

toc

figure
hold on
plot(1:length(g_array), total_fisher_info/(no_atoms*T^2),LineWidth=3)
%plot(1:length(g_array), 4-3*g_array.^2*T/(tau),LineWidth=3)
plot(1:length(g_array), fi_fit/(no_atoms*T^2),LineWidth=3)
%set(gca, 'XScale', 'log')
ax = gca;
ax.FontSize = 22; 
ax.FontName = "times";
xlabel('Measurement strength g',FontSize=26,FontName="times")
ylabel('FI/T^2',FontSize=26,FontName="times")

disp("FI sampling finished");


%% FISHER INFO averaged over w, wrt T

close all

no_atoms = 10;
tau = 0.1;

g = 0.1; % measurement strength

T_arr = 10.^(log10(1.):0.05:log10(100.));

no_trajectory_samples = 1e3;
w_array = 0.001:pi/(2*tau)/40:pi/(2*tau)-0.001;

total_fisher_info = zeros(length(T_arr), length(w_array));

strong_meas = 0;
classical = 0;
p_e = 0; % error probability in a weak measurement

tic

for j = 1:length(T_arr)

    T = T_arr(j);

    parfor i = 1:length(w_array)
    
        total_fisher_info(j, i) = avg_fisher_information_per_trajectory(g, ...
            no_trajectory_samples, T, tau, w_array(i), no_atoms, ...
            strong_meas, classical, p_e);
    
    end

end

toc

avg_fisher_info = mean(total_fisher_info, 2);

s = struct("no_atoms", no_atoms, "T_arr", T_arr, "tau", tau,...
    "g", g, "avg_fisher_per_atom", avg_fisher_info.'/no_atoms, "w_array",...
    w_array);
save(targetPath + "avg_fisher_info_wrt_T_no_strong.mat", '-struct', 's');

figure
hold on
plot(T_arr*g^2/tau, avg_fisher_info.'./(4*no_atoms*T_arr.^2),'-o', LineWidth=2)
ax = gca;
ax.FontSize = 22; 
ax.FontName = "times";
xlabel('T (units of tau/g^2)',FontSize=26,FontName="times")
ylabel('FI/T^2',FontSize=26,FontName="times")
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
%xticks([0.002, 0.01, 0.1, 0.5])

disp("FI sampling finished");



%% CLASSICAL CASE

close all

T = 20;  % total time
tau = 0.1;  %weak meas time
w = 1;

no_atoms = 10;

no_samples = 100;
g = 0.1;
g_array = [0.1];

strong_meas = 0;
p_e = 0;

classical = 0;

[mle_variance_1] = variance_of_mle(g_array, no_samples, T, ...
    tau, w, no_atoms, strong_meas, classical, p_e);

classical = 1;

[mle_variance_2] = variance_of_mle(g_array, no_samples, T, ...
    tau, w, no_atoms, strong_meas, classical, p_e);

disp(mle_variance_1);
disp(mle_variance_2);


%% Delta Omega wrt TIME

close all

omega_corner = pi; % from -pi to pi
delta_omega = omega_corner;

N = 64;

strong_meas = 1;
classical = 0;

no_samples = 1e3;

tau_default = pi/(2*omega_corner);

% T_1 = 10.^(log10(0.1):0.17:log10(1.6))/delta_omega;
% T_arr = cat(2, [T_1, tau_default*round(10.^(log10(2):0.2:log10(75)))]);

T_arr = tau_default*round(10.^(log10(2):0.2:log10(75)));

mle = zeros([1, length(T_arr)]);
mle_variance = zeros([1, length(T_arr)]);
w_bmse = zeros([1, no_samples]);

% tic
% 
% optimal_gs = optimize_over_g(T_arr, tau_default, omega_corner, N, strong_meas);
% 
% toc

%load(path + "delta_omega_our_protocol_optimal_gs.mat");

tic

%optimal_gs = cat(2,zeros(1,length(T_1)),...
%    [0.2622,0.2431,0.213,0.1718,0.14,0.115,0.093,0.077]);

optimal_gs = [0.2622,0.2431,0.213,0.1718,0.14,0.115,0.093,0.077];


for i = 1:length(T_arr)

    T = T_arr(i);

    if pi/(2*omega_corner)+1e-5 > T

        tau = T;

    else

        tau = tau_default;

    end

    g_opt = optimal_gs(i);

    disp("T = " + T + ", tau = " + tau + ", g = " + g_opt);

    w_samples = (rand([1, no_samples]))*omega_corner;

    parfor j=1:no_samples

        w_bmse(j) = variance_of_mle([g_opt], 1, T, ...
            tau, w_samples(j), N, strong_meas, classical, 0);

    end

    mle(i) = mean(w_bmse, "all");
    mle_variance(i) = var(w_bmse,0, "all");

    disp(sqrt(mle(i))*2*sqrt(N)*T_arr(i));


end

toc

s = struct("no_atoms", N, "T_arr", T_arr, ...
   "bmse", mle, "optimal_gs", optimal_gs, "omega_corner",...
   omega_corner);
save(targetPath + "/delta_omega_our_protocol.mat", '-struct', 's');

%load(targetPath + "/rosenband_freq_optimal_ensemble.mat"); 
% load(path + "delta_omega_our_protocol.mat");
% mle_variance = bmse;

disp(mean(sqrt(mle)*2*sqrt(N).*T_arr));


figure
hold on
plot(delta_omega*T_arr, 1./(2*T_arr)/delta_omega, LineWidth=3)
plot(delta_omega*T_arr, sqrt(mle)*sqrt(N)/delta_omega, '-o', LineWidth=3)
%plot(delta_omega*rosenband_T_arr, sqrt(rosenband_bmse)*2*sqrt(N), LineWidth=3)
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
ax = gca;
ax.FontSize = 20; 
ax.FontName = "times";
xlabel("$$\delta\omega T$$", 'Interpreter','latex')
ylabel('$$\sqrt N \Delta\omega/\delta\omega$$', 'Interpreter','latex');
h = legend('$$1/T$$', "Our protocol", Location='southwest');
set(h,'Interpreter','latex') 




%% FI vs MLE wrt MEASUREMENT STRENGTH

close all

no_atoms = 20

tau = 0.1;
w = 5.;
strong_meas = 0;
classical = 0;

no_samples = 1e3;

T = 10;

g_arr = 10.^(-2.:0.07:-0.10);

good_probs = squeeze(calculate_good_prob(T, tau, no_atoms, g_arr)).';


fi_analytical = fi_analytical_get(g_arr, T, tau, strong_meas, classical);

mse_analytical = 1./(fi_analytical*no_atoms).*good_probs+...
    (1-good_probs).*pi^2./(48*tau^2);


if isempty(gcp('nocreate'))

    parpool;

end

tic

total_fisher_info = avg_fisher_information_per_trajectory(g_arr, ...
     no_samples, T, tau, w, no_atoms, strong_meas, classical, 0)/no_atoms;

mle_variance = variance_of_mle(g_arr, no_samples, T, ...
          tau, w, no_atoms, strong_meas, classical, 0)*no_atoms;

toc


s = struct("tau_mle", tau, "g_mle" , g_arr, "T_mle", T, ...
                   "mle_variance_mle", mle_variance, "w",w);
save(targetPath+ "/mle_variance_wrt_g_" + num2str(no_atoms) + ...
     "_atoms.mat", '-struct', 's');


figure
hold on
plot(g_arr.^2*T/tau, total_fisher_info./(4*T^2), '-o', ...
    Color=[255,218,8]/256,LineWidth=3)
plot(g_arr.^2*T/tau, fi_analytical./(4*T^2), Color=[10,10,255]/256, ...
    LineWidth=3)
plot(g_arr.^2*T/tau, 1./(mle_variance)/(4*T^2), LineWidth=3, ...
     Color=[245,42,177]/256)
plot(g_arr.^2*T/tau, 1./(mse_analytical*no_atoms)/(4*T^2), LineWidth=3, ...
    Color=[50, 168, 82]/256)
set(gca, 'XScale', 'log')
ax = gca;
ax.FontSize = 14; 
ax.FontName = "times";
xlabel('g^2 T/tau',FontName="times")
ylabel('FI/(4NT^2)',FontName="times")
legend("FI","FI Fit", "1/MSE", "1/MSE Analytic Fit", FontSize=12, ...
    FontName="times", Location='northeast')
title(strcat("N=",num2str(no_atoms), " T=", num2str(T)))






%% FI vs MLE wrt TIME

close all

no_atoms = 64;

tau = 0.1;
w = 1.;
strong_meas = 1;
classical = 0;

no_samples = 1e3;

g = sqrt(0.01);

T_arr = 10.^(log10(1.):0.1:log10(100.));

total_fisher_info = zeros(1, length(T_arr));
fi_quantum = zeros(1, length(T_arr));
mle_variance = zeros(1, length(T_arr));

good_probs = calculate_good_prob(T_arr, tau, no_atoms, g(1));

fi_analytical = fi_analytical_get(g, T_arr, tau, strong_meas, classical)

mse_analytical = 1./(fi_analytical*no_atoms).*good_probs+...
    (1-good_probs).*pi^2./(48*tau^2);


if isempty(gcp('nocreate'))

    parpool;

end


tic

parfor i = 1:length(T_arr)

    total_fisher_info(i) = avg_fisher_information_per_trajectory([g], ...
     no_samples, T_arr(i), tau, w, no_atoms, strong_meas, classical, 0)/no_atoms;

    mle_variance(i) = variance_of_mle(g, no_samples, T_arr(i), ...
            tau, w, no_atoms, strong_meas, classical, 0)*no_atoms;
    
end


toc


figure
hold on
plot(g^2*T_arr/tau, total_fisher_info./(4*T_arr.^2), '-o', ...
    Color=[255,218,8]/256,LineWidth=3)
plot(g^2*T_arr/tau, fi_analytical./(4*T_arr.^2), Color=[10,10,255]/256, ...
    LineWidth=3)
plot(g^2*T_arr/tau, 1./(mle_variance)./(4*T_arr.^2), ...
    LineWidth=3, Color=[245,42,177]/256)
plot(g^2*T_arr/tau, 1./(mse_analytical)./(4*T_arr.^2*no_atoms), ...
    LineWidth=3, Color=[50, 168, 82]/256)
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
ax = gca;
ax.FontSize = 14; 
ax.FontName = "times";
xlabel('g^2 T/tau',FontName="times")
ylabel('I/(4*N*T^2)',FontName="times")
legend("FI","FI Fit", "1/MSE", "1/MSE Analytic Fit", FontSize=12, ...
    FontName="times", Location='northeast')
title(strcat("N=",num2str(no_atoms), " g=", num2str(g(1))))



%% IMPERFECT MEASUREMENTS

close all

T_arr = 10.^(log10(1.):0.1:log10(100.));  % total time
tau = 0.1;  %weak meas time
w = 5.;

no_atoms = 10;

no_trajectory_samples = 5e3;
g = 0.1;

strong_meas = 0;
classical = 0;
p_e_array = [0., 0.05, 0.1];

total_fisher_info = zeros([length(p_e_array), length(T_arr)]);

fi_analytical = fi_analytical_get(g*(1-2*p_e_array).', T_arr, tau, ...
    strong_meas, classical);

tic


for i=1:length(T_arr)

    for j = 1:length(p_e_array)

        total_fisher_info(j, i) = avg_fisher_information_per_trajectory([g], ...
            no_trajectory_samples, T_arr(i), tau, w, no_atoms, strong_meas, ...
            classical, p_e_array(j))/no_atoms;

    end

end

toc

s = struct("tau", tau, "g" ,g, "T_arr", T_arr, "p_e_array", p_e_array, ...
                   "total_fisher_info", total_fisher_info);
save(targetPath + "/fisher_info_imperf_wrt_T_no_strong.mat", ...
     '-struct', 's');


figure
hold on
plot(g.^2*T_arr/tau, total_fisher_info./(4*T_arr.^2), LineWidth=3)
plot(g.^2*T_arr/tau, fi_analytical./(4*T_arr.^2), ...
    LineWidth=3, LineStyle="--")
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
ax = gca;
ax.FontSize = 16; 
ax.FontName = "times";
xlabel('$g^2T/\tau$',FontSize=20,FontName="times",Interpreter="latex")
ylabel('FI/N',FontSize=20,FontName="times",Interpreter="latex")
legend("p_e = " + p_e_array)

disp("FI sampling finished");



%% Outlier prob approximations


no_atoms = 110;

tau = 0.1;
w = 5.;
strong_meas = 0;
classical = 1;

no_samples = 8000;

g = [0.055];

T_arr = 10.^(log10(5.):0.05:log10(40.));

disp("True prob:");

disp(calculate_good_prob(T_arr, tau, no_atoms, g(1)));

disp("Approximation:");

disp(calculate_good_prob_approx(T_arr, tau, no_atoms, g(1)));



function good_probs = calculate_good_prob(T_arr, tau, no_atoms, g)

    L = 200;
    good_probs = zeros(1, length(T_arr), length(g));

    for i = 1:length(T_arr)
        for j = 1:length(g)

            meas_g = g(j);
            T = T_arr(i);
            sigma = sqrt(tau/(8*no_atoms*T));
            xBar = meas_g/2;
    
            dx = 24*sigma/(L-1);
    
            x_array = xBar+(-12*sigma:dx:12*sigma);
    
            good_probs(1, i, j) = (dx*sqrt(1/(2*pi*sigma^2))*sum(exp(-(x_array-...
                xBar).^2/(2*sigma^2)).*(normcdf(x_array, 0, sigma)).^(T/(2*tau)-1)));

        end

    end

    %good_probs = squeeze(good_probs);

end


function good_probs = calculate_good_prob_approx(T_arr, tau, no_atoms, g)

    L = 200;
    good_probs = zeros(1, length(T_arr));

    no_atoms = no_atoms + 10;

    for i = 1:length(T_arr)

        T = T_arr(i);
        sigma = sqrt(tau/(8*no_atoms*T));
        xBar = g/2;

        dx = 24*sigma/(L-1);

        x_array = xBar+(-12*sigma:dx:12*sigma);

        good_probs(i) = 1 - 0.046*T/tau*exp(-0.5025*g^2*no_atoms*T/tau);

    end

end


% Analytical expression for fisher information per atom
function fi = fi_analytical_get(g, T, tau, strong_meas, classical)


    fi_weak_quantum = (8*g.^2.*T.^3)/(3*tau)./...
        (1+3/16*sin(2*g).^2.*T./tau +g.^2./3.*log(sec(2*g)).*T.^2/tau^2);
    fi_weak_classical = (8*g.^2.*T.^3)/(3*tau);
    
    fi_strong_quantum = 4*T.^2./(1+T/(2*tau).*log(sec(2*g)));
    fi_strong_classical = 4*T.^2+(8*g.^2*T.^3)/(3*tau);
    
    fi_analytical_quantum = strong_meas*fi_strong_quantum + (1-strong_meas)*fi_weak_quantum;
    fi_analytical_classical = strong_meas*fi_strong_classical + (1-strong_meas)*fi_weak_classical;
    
    fi = classical*fi_analytical_classical + (1-classical)*fi_analytical_quantum;

end