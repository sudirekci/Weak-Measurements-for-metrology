
disp("------------------------------------------------");

T = 10;  % total time
tau = 0.1;  %weak meas time
w = 1;

g = 0.02;

no_atoms = 5;
no_trajectory_samples = 10;
strong_meas = 1;

[outcomes, prob_trajectory] = sample_trajectory(T, tau, g, rand(1, no_trajectory_samples), ...
    no_atoms, no_trajectory_samples, strong_meas);

[prob_w_given_tr] = simulate_given_trajectory(outcomes, ...
    w*ones(1, no_trajectory_samples), T, tau, g, no_atoms, ...
    no_trajectory_samples, strong_meas);


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


%% MLE ESTIMATION

%close all

no_atoms = 100;
T = 20;
w = 1;
tau = 0.1;

no_trajectory_samples = 200;

g_array = 10.^(-3.:0.1:-0.1);

strong_meas = 1;

tic


classical = 0;

[total_fisher_info] = avg_fisher_information_per_trajectory(g_array, ...
    no_trajectory_samples, T, tau, w, no_atoms, strong_meas, classical);

classical = 0;

[mle_variance] = variance_of_mle(g_array, no_trajectory_samples, T, ...
    tau, w, no_atoms, strong_meas, classical);


toc

figure
hold on
plot(g_array, 1./(mle_variance.*T^2.*no_atoms),LineWidth=3,...
    Color=[245,42,177]/256)
plot(g_array, smooth(total_fisher_info)./(no_atoms*T^2), ...
    Color=[255,218,8]/256, LineWidth=3)
set(gca, 'XScale', 'log')
ax = gca;
ax.FontSize = 22; 
ax.FontName = "times";
xlabel('Measurement Strength g',FontSize=24,FontName="times")
ylabel('1/Var(MLE)NT^2',FontSize=24,FontName="times")
legend("1/MLE","Fisher",FontSize=20,FontName="times")
%xlim([min(sin(2*g_array)), max(sin(2*g_array))])
%xticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6 ])

disp("MLE variance finished");

%% FI vs MLE wrt TIME


%close all

no_atoms = 50

tau = 0.1;
w = 5.;
strong_meas = 0;
classical = 0;

no_samples = 2e5;

g = sqrt(0.01);

T_arr = 10.^(log10(1.):0.05:log10(100.));
%T_arr = [69.1831   83.1764  100.0000];

total_fisher_info = zeros(1, length(T_arr));
fi_quantum = zeros(1, length(T_arr));
mle_variance = zeros(1, length(T_arr));

good_probs = calculate_good_prob(T_arr, tau, no_atoms, g(1));

fi_weak_quantum = (8*g^2*T_arr.^3)/(3*tau)./...
    (1+ 3/4*g^2*T_arr./tau +2/3*g^4*T_arr.^2/tau^2);
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

for i = 1:length(T_arr)

    total_fisher_info(i) = avg_fisher_information_per_trajectory([g], ...
       no_samples, T_arr(i), tau, w, no_atoms, strong_meas, classical);

    % mle_variance(i) = variance_of_mle(g, no_samples, T_arr(i), ...
    %       tau, w, no_atoms, strong_meas, classical);

    % disp(1/(mle_variance(i)*no_atoms*fi_analytical(i)));
    disp(T_arr(i));

end

toc

% p = gcp;
% delete(p);

total_fisher_info = total_fisher_info/no_atoms;
s = struct("tau_fisher", tau, "g_fisher" ,g, "T_arr_fisher", T_arr, "total_fisher_info", total_fisher_info);
save('/Users/sudirekci/Documents/MATLAB/weak_measurements/data/fisher_info_no_strong.mat', '-struct', 's');

% s = struct("tau_mle", tau, "g_mle" ,g, "T_arr_mle", T_arr, "no_atoms", no_atoms, "mle_variance", mle_variance);
% save('C:\Users\Endreslab\Documents\MATLAB\Su\weak_meas_data\mle_variance_no_strong_20_atoms.mat', '-struct', 's');


figure
hold on
plot(g^2*T_arr/tau, log(total_fisher_info/no_atoms), Color=[255,218,8]/256, LineWidth=3)
plot(g^2*T_arr/tau, log(fi_analytical), Color=[10,10,255]/256, LineWidth=3)
plot(g^2*T_arr/tau, log(1./(mle_variance.*no_atoms)), LineWidth=3, Color=[245,42,177]/256)
plot(g^2*T_arr/tau, log(1./(mse_analytical.*no_atoms)), LineWidth=3, Color=[50, 168, 82]/256)
set(gca, 'XScale', 'log')
ax = gca;
ax.FontSize = 14; 
ax.FontName = "times";
xlabel('g^2 T/tau',FontName="times")
ylabel('ln(I/N)',FontName="times")
legend("FI fit", "1/MSE","1/MSE Analytic Fit", FontSize=12, FontName="times", Location='northwest')
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
    good_probs = zeros(1, length(T_arr));

    for i = 1:length(T_arr)

        T = T_arr(i);
        sigma = sqrt(tau/(8*no_atoms*T));
        xBar = g/2;

        dx = 24*sigma/(L-1);

        x_array = xBar+(-12*sigma:dx:12*sigma);

        good_probs(i) = (dx*sqrt(1/(2*pi*sigma^2))*sum(exp(-(x_array-...
            xBar).^2/(2*sigma^2)).*(normcdf(x_array, 0, sigma)).^(T/(2*tau)-1)));

    end

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



function good_probs = calculate_good_prob_chi(T_arr, tau, no_atoms, g)

    L = 250;
    good_probs = zeros(1, length(T_arr));

    for i = 1:length(T_arr)

        T = T_arr(i);

        xBar = g/2;
        sigma_r = sqrt(tau/(8*no_atoms*T));

        dx = (6*sigma_r+xBar)/(L-1);

        x_array = 0:dx:6*sigma_r+xBar;

        good_probs(i) = dx*sum((1-exp(-x_array.^2/(2*sigma_r^2))).^(T/(2*tau)-1).*...
           x_array/sigma_r^2.*exp(-(x_array.^2+xBar^2)./(2*sigma_r^2)).*...
           besseli(0,x_array*xBar/sigma_r^2));

        %good_probs(i) = 1 - T/(4*tau)*exp(-xBar^2/(4*sigma_r^2));

    end

end
