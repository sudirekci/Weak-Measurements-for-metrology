
%close all

noSamples = 1e3;

N = 64;

omega_corner = pi; % from -pi/2 to pi/2
delta_omega = omega_corner;

T_arr = 10.^(-1:0.05:log10(80))/delta_omega;
%T_arr = [4.];
%delta_phi_arr = [1.2];

classical_bmse = zeros([1, length(T_arr)]);

rosenband = true;

tic

for l=1:length(T_arr)

    T = T_arr(l);
    %no_ensembles = floor(max(log2(T./(0.5)), -0.5)+1) + 1;
    no_ensembles =1;
    speeds = 2.^((no_ensembles-1):-1:0);

    if rosenband

        no_atoms = floor(N/no_ensembles)*ones([1, no_ensembles]);
        no_atoms(end) = no_atoms(end) + N - sum(no_atoms);
    else
        no_atoms = [44, 20];
    end

    classical_bmse(l) = cascade_sim_semianalytic(T, speeds, no_atoms, ...
        noSamples, omega_corner);

end

toc


figure
hold on
plot(delta_omega*T_arr, 1./(2*T_arr*delta_omega), LineWidth=3)
plot(delta_omega*T_arr, sqrt(classical_bmse)*sqrt(N)/delta_omega, LineWidth=3)
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
ax = gca;
ax.FontSize = 20; 
ax.FontName = "times";
xlabel("\delta\omega T")
lab = ylabel('$$\sqrt N \Delta\omega/\delta\omega$$');
set(lab,'Interpreter','latex') 
title("# Samples = " + num2str(noSamples))



function [bmse] = cascade_sim_freq(T, speeds, no_atoms, noSamples, omega_corner)

    w_samples = ((rand([noSamples, 1])-0.5)*omega_corner*T);

    %optimal_estimators = single(zeros([1, ceil(prod((no_atoms/2+1).^2))+5]));
    %trajectoryCount = single(zeros([1, ceil(prod((no_atoms/2+1).^2))+5]));

    r = (zeros([noSamples, 1]));

    for l = 1:length(speeds)

        num_atoms = no_atoms(l);
        num_1 = (floor(num_atoms/2));
        num_2 = (num_atoms - num_1);

        if num_atoms == 0
            continue
        end

        % r = (num_1+1)*r + binornd(num_1, (1+sin(2*w_samples./speeds(l)))/2);
        % r = (num_2+1)*r + binornd(num_2, (1+cos(2*w_samples./speeds(l)))/2);

        r = (num_1+1)*(num_2+1)*r + (num_2+1)*(binornd(num_1, ...
            (1+sin(2*w_samples./speeds(l)))/2)) + (binornd(num_2, ...
            (1+cos(2*w_samples./speeds(l)))/2));

    end

    r = r + 1;

    [G, ~] = findgroups(r);

    bmse = sum(accumarray(G, w_samples, [], @(x) sum((x - ...
        mean(x)).^2)))/(noSamples*T^2);

end



function [bmse] = cascade_sim_semianalytic(T, speeds, no_atoms, noSamples, ...
    omega_corner)

    phi_samples = ((rand([1, noSamples])-0.5)*omega_corner*T);

    freq_res = 350;
    phi_array = (-omega_corner/2:omega_corner/freq_res:omega_corner/2)*T;

    prob_funcs_plus = zeros(length(speeds)*2, freq_res+1);
    prob_funcs_minus = zeros(length(speeds)*2, freq_res+1);
    
    for l=1:length(speeds)

        prob_funcs_plus(2*l-1, :) = log((1+sin(2*phi_array./speeds(l)))/2);
        prob_funcs_plus(2*l, :) = log((1+cos(2*phi_array./speeds(l)))/2);

        prob_funcs_minus(2*l-1, :) = log((1-sin(2*phi_array./speeds(l)))/2);
        prob_funcs_minus(2*l, :) = log((1-cos(2*phi_array./speeds(l)))/2);

    end

    r = zeros([length(speeds)*2, 1, noSamples]);
    num_arr = zeros([length(speeds)*2, 1, 1]);

    for l = 1:length(speeds)

        num_atoms = no_atoms(l);
        num_1 = (floor(num_atoms/2));
        num_2 = (num_atoms - num_1);

        num_arr(2*l-1) = num_1;
        num_arr(2*l) = num_2;

        if num_atoms == 0
            continue
        end

        r(2*l-1, 1, :) = binornd(num_1, (1+sin(2*phi_samples./speeds(l)))/2);
        r(2*l, 1, :) = binornd(num_2, (1+cos(2*phi_samples./speeds(l)))/2);

    end

    prob_funcs_plus = r.*prob_funcs_plus;
    prob_funcs_minus = (num_arr-r).*prob_funcs_minus;

    prob_funcs_plus(isnan(prob_funcs_plus)) = 0;
    prob_funcs_minus(isnan(prob_funcs_minus)) = 0;

    likelihood = squeeze(sum(prob_funcs_plus + prob_funcs_minus,1)); 

    likelihood = exp(likelihood-max(likelihood)); % size: length(phi_array) x noSamples
    likelihood = likelihood./sum(likelihood, 1);

    bmse = sum((phi_samples - phi_array*likelihood).^2)/(noSamples*T^2);

end

