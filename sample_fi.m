
disp("------------------------------------------------");


%% FISHER INFO SAMPLING

%close all

T = 20;  % total time
tau = 0.1;  %weak meas time
w = 1;

no_atoms = 10;

no_trajectory_samples = 6000;
g_array = 10.^(-2.:0.07:-0.2);

strong_meas = 0;
classical = 0;

fi_fit = (8*g_array.^2.*no_atoms*T^3)/(3*tau)./(1+3/4*g_array.^2*T/tau+ ...
                                +2/3*g_array.^4*T^2/tau^2);

tic

[total_fisher_info] = avg_fisher_information_per_trajectory(g_array, ...
    no_trajectory_samples, T, tau, w, no_atoms, strong_meas, classical);

toc

% 2*tanh(-1/2*log(g_array.^2*T/tau))+2

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
%xticks([0.002, 0.01, 0.1, 0.5])

disp("FI sampling finished");


%% CLASSICAL CASE

close all

T = 20;  % total time
tau = 0.1;  %weak meas time
w = 1;

no_atoms = 10;

no_samples = 500;
g = 0.1;
g_array = [0.1];

strong_meas = 0;

classical = 0;

[mle_variance_1] = variance_of_mle(g_array, no_samples, T, ...
    tau, w, no_atoms, strong_meas, classical);

classical = 1;

[mle_variance_2] = variance_of_mle(g_array, no_samples, T, ...
    tau, w, no_atoms, strong_meas, classical);


%% Delta Omega wrt TIME

%close all

path = "/Users/sudirekci/Documents/MATLAB/weak_measurements/data/";

omega_corner = pi; % from -pi to pi
delta_omega = omega_corner/sqrt(12);

N = 64;

strong_meas = 1;
classical = 0;

no_samples = 2e4;
noRep = 4;

T_arr = 10.^(log10(0.1):0.1:log10(10.))/delta_omega;
T_arr = [0.438976];
mle_variance = zeros([1, length(T_arr)]);
w_bmse = zeros([noRep, no_samples]);

tau_default = 0.1;

% tic
% 
% optimal_gs = optimize_over_g(T_arr, tau_default, omega_corner, N, strong_meas);
% 
% s = struct("no_atoms", N, "T_arr", T_arr, ...
%     "optimal_gs", optimal_gs);
% save(path + "delta_omega_our_protocol_optimal_gs.mat", '-struct', 's');
% 
% toc

%load(path + "delta_omega_our_protocol_optimal_gs.mat");

tic

for i = 1:length(T_arr)

    T = T_arr(i);

    % pi/(2*omega_corner*T) > T

    if true

        tau = T;

    else

        tau = T/100.0001;

    end

    %g_opt = optimal_gs(i);

    g_opt = 0.058;

    disp("T = " + T + ", tau = " + tau);

    w_samples = (rand([1, no_samples])-0.5)*omega_corner+pi/(4*tau);

    for l = 1:noRep
        parfor j=1:no_samples
    
            w_bmse(l, j) = variance_of_mle([g_opt], 1, T, ...
                tau, w_samples(j), N, strong_meas, classical);
    
        end
    end

    if i > 7
        figure
        scatter(w_samples, w_bmse(1, :))
        xlabel('\omega')
        ylabel("BMSE")
    end


    mle_variance(i) = mean(w_bmse, "all");

    disp(mle_variance(i));
    disp(sqrt(mle_variance(i))*2*sqrt(N)*T_arr(i));


end

toc

% s = struct("no_atoms", N, "T_arr", T_arr, ...
%     "bmse", mle_variance, "optimal_gs", optimal_gs);
% save(path + "delta_omega_our_protocol.mat", '-struct', 's');

% load(path + "rosenband_freq_optimal_ensemble.mat");
% 
% figure
% hold on
% plot(delta_omega*T_arr, 1./(T_arr), LineWidth=3)
% plot(delta_omega*T_arr, sqrt(mle_variance)*2*sqrt(N), '-o', LineWidth=3)
% plot(delta_omega*rosenband_T_arr, sqrt(rosenband_bmse)*2*sqrt(N), LineWidth=3)
% set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')
% ax = gca;
% ax.FontSize = 20; 
% ax.FontName = "times";
% xlabel("\delta\omega T")
% lab = ylabel('2$$\sqrt N \Delta\omega$$');
% set(lab,'Interpreter','latex') 
% h = legend('$$1/T$$', "Our protocol", "Rosenband",Location='southwest');
% set(h,'Interpreter','latex') 




%% FI vs MLE wrt MEASUREMENT STRENGTH

path = "/Users/sudirekci/Documents/MATLAB/weak_measurements/data/";

%close all

no_atoms = 20

tau = 0.1;
w = 5.;
strong_meas = 0;
classical = 0;

if strong_meas
    meas = "with";
else
    meas = "no";
end

no_samples = 5e4
noRep = 10;

T = 10;

% g_arr = cat(2, 10.^(-2.:0.05:-0.15), pi/4);
g_arr = 10.^(-2.:0.07:-0.10);
g_arr = [0.112201845430196];

good_probs = squeeze(calculate_good_prob(T, tau, no_atoms, g_arr)).';


fi_weak_quantum = (8*g_arr.^2*T.^3)/(3*tau)./...
    (1+3/16*sin(2*g_arr).^2*T./tau +g_arr.^2./3.*log(sec(2*g_arr))*T.^2/tau^2);

fi_weak_classical = (8*g_arr.^2*T.^3)/(3*tau);

fi_strong_quantum = 4*T.^2./(1+T/(2*tau)*log(sec(2*g_arr)));

fi_strong_classical = 4*T.^2+(8*g_arr.^2*T.^3)/(3*tau);

fi_analytical_quantum = strong_meas*fi_strong_quantum + (1-strong_meas)*fi_weak_quantum;
fi_analytical_classical = strong_meas*fi_strong_classical + (1-strong_meas)*fi_weak_classical;

fi_analytical = classical*fi_analytical_classical + (1-classical)*fi_analytical_quantum;

mse_analytical = 1./(fi_analytical_quantum*no_atoms).*good_probs+...
    (1-good_probs).*pi^2./(48*tau^2);


if isempty(gcp('nocreate'))

    parpool;

end

tic

% total_fisher_info = avg_fisher_information_per_trajectory(g_arr, ...
%   no_samples, T, tau, w, no_atoms, strong_meas, classical)/no_atoms;


mle_variance = zeros([noRep, length(g_arr)]);

for j =1:noRep

    j

    mle_variance(j, :) = variance_of_mle(g_arr, no_samples, T, ...
          tau, w, no_atoms, strong_meas, classical);

end


toc

mle_variance = mean(mle_variance, 1);

mle_variance
1./(mle_variance.*no_atoms)/(4*T^2)


% s = struct("tau_fisher", tau, "g_arr_fisher" ,g_arr, "T_fisher", T, ...
%                   "total_fisher_info", total_fisher_info);
% save('/Users/sudirekci/Documents/MATLAB/weak_measurements/data/fisher_info_wrt_g_with_strong.mat', ...
%     '-struct', 's');


% s = struct("tau_mle", tau, "g_mle" , g_arr, "T_mle", T, ...
%                   "mle_variance_mle", mle_variance);
% save(path+ "mle_variance_wrt_g_" + meas +"_strong_" + num2str(no_atoms) + ...
%     "_atoms.mat", '-struct', 's');


% figure
% hold on
% %plot(g_arr.^2*T/tau, total_fisher_info/(4*T^2), Color=[255,218,8]/256, LineWidth=3)
% plot(g_arr.^2*T/tau, fi_analytical/(4*T^2), Color=[10,10,255]/256, LineWidth=3)
% plot(g_arr.^2*T/tau, 1./(mle_variance.*no_atoms)/(4*T^2), LineWidth=3, ...
%     Color=[245,42,177]/256)
% plot(g_arr.^2*T/tau, 1./(mse_analytical.*no_atoms)/(4*T^2), LineWidth=3, ...
%    Color=[50, 168, 82]/256)
% set(gca, 'XScale', 'log')
% ax = gca;
% ax.FontSize = 14; 
% ax.FontName = "times";
% xlabel('g^2 T/tau',FontName="times")
% ylabel('FI/(4NT^2)',FontName="times")
% legend("FI Quantum", "FI Classical", "1/MSE","1/MSE Analytic Fit", FontSize=12, FontName="times", Location='northwest')
% title(strcat("N=",num2str(no_atoms), " T=", num2str(T)))






%% FI vs MLE wrt TIME

path = "/Users/sudirekci/Documents/MATLAB/weak_measurements/data/";


%close all

no_atoms = 20

tau = 0.1;
w = 5.;
strong_meas = 0;
classical = 0;

no_samples = 2e5
noReps = 1;

g = sqrt(0.01);

T_arr = 10.^(log10(1.):0.05:log10(100.));
%T_arr = [1.20226443461741];

total_fisher_info = zeros(1, length(T_arr));
fi_quantum = zeros(1, length(T_arr));
mle_variance = zeros(length(T_arr), noReps);

good_probs = calculate_good_prob(T_arr, tau, no_atoms, g(1));

fi_weak_quantum = (8*g^2*T_arr.^3)/(3*tau)./...
    (1+ 0.66*g^2*T_arr./tau +2/3*g^4*T_arr.^2/tau^2);
fi_weak_classical = (8*g^2*T_arr.^3)/(3*tau);

fi_strong_quantum = 4*T_arr.^2./(1+g.^2.*T_arr/tau);

fi_strong_classical = 4*T_arr.^2+(8*g^2*T_arr.^3)/(3*tau);

fi_analytical_quantum = strong_meas*fi_strong_quantum + (1-strong_meas)*fi_weak_quantum;
fi_analytical_classical = strong_meas*fi_strong_classical + (1-strong_meas)*fi_weak_classical;

fi_analytical = classical*fi_analytical_classical + (1-classical)*fi_analytical_quantum;

mse_analytical = 1./(fi_analytical_classical*no_atoms).*good_probs+...
    (1-good_probs).*pi^2./(48*tau^2);


if isempty(gcp('nocreate'))

    parpool;

end

tic

for j = 1:noReps
    for i = 1:length(T_arr)

        g = sqrt(tau/T_arr(i))*(3/2)^(1/4);

        total_fisher_info(i) = avg_fisher_information_per_trajectory([g], ...
          no_samples, T_arr(i), tau, w, no_atoms, strong_meas, classical)/no_atoms;
    
        % mle_variance(i, j) = variance_of_mle(g, no_samples, T_arr(i), ...
        %       tau, w, no_atoms, strong_meas, classical);
    
        %disp(1/(mle_variance(i, j)*no_atoms*fi_analytical(i)));
        % disp(T_arr(i));
        
    end

    % disp(log(1./(mle_variance(:, j).*no_atoms)));

    %disp("*****************")
end

% mle_variance = mean(mle_variance, 2);

% disp(mle_variance);
% disp(log(1./(mle_variance.*no_atoms)));
% disp(total_fisher_info);

toc

% p = gcp;
% delete(p);

s = struct("tau_fisher", tau, "g_fisher" ,g, "T_arr_fisher", T_arr, ...
                  "total_fisher_info", total_fisher_info);
save('/Users/sudirekci/Documents/MATLAB/weak_measurements/data/fisher_info_no_strong_opt_g.mat', ...
    '-struct', 's');


figure
hold on
plot(g^2*T_arr/tau, total_fisher_info, Color=[255,218,8]/256, LineWidth=3)
plot(g^2*T_arr/tau, fi_weak_quantum, Color=[10,10,255]/256, LineWidth=3)
plot(g^2*T_arr/tau, 4*T_arr.^2, Color=[2,218,8]/256, LineWidth=3)
%plot(g^2*T_arr_mle/tau, log(1./(mle_variance.*no_atoms)), LineWidth=3, Color=[245,42,177]/256)
%plot(g^2*T_arr/tau, log(1./(mse_analytical.*no_atoms)), LineWidth=3, Color=[50, 168, 82]/256)
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
ax = gca;
ax.FontSize = 14; 
ax.FontName = "times";
xlabel('g^2 T/tau',FontName="times")
ylabel('ln(I/N)',FontName="times")
%legend("FI Quantum", "FI Classical", "1/MSE","1/MSE Analytic Fit", FontSize=12, FontName="times", Location='northwest')
title(strcat("N=",num2str(no_atoms), " g=", num2str(g(1))))






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
