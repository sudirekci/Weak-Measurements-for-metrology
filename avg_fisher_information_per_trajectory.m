

% This function computes the average Fisher information per trajectory 
% (not the average of Fisher information over all trajectories)given some
% g_array. Set dw << 1.

% strong_meas: 1 -> perform strong projective measurement at the end

function [total_fisher_info] = avg_fisher_information_per_trajectory(g_array, ...
    no_trajectory_samples, T, tau, w, no_atoms, strong_meas, classical)

    dw = 10^(-8.);

    total_fisher_info = zeros(1, length(g_array));

    if classical

        sample_fun = @sample_classical_trajectory;
        simulate_fun = @simulate_classical_trajectory;

    else

        sample_fun = @sample_trajectory;
        simulate_fun = @simulate_given_trajectory_2;

    end

    parfor g_i=1:length(g_array)

        g = g_array(g_i);
        
        [outcomes, prob_trajectory] = sample_fun(T, tau, g, w, ...
            no_atoms, no_trajectory_samples, strong_meas);

        [prob_dw] = simulate_fun(outcomes, w+dw, ...
            T, tau, g, no_atoms, no_trajectory_samples, strong_meas);

        total_fisher_info(g_i) = sum((exp(prob_dw-prob_trajectory)-1).^2, ...
            "all")/dw^2/no_trajectory_samples;
    
    end


end