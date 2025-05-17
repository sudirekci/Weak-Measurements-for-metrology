

function [mle_variance] = variance_of_mle(g_array, no_samples, ...
    T, tau, w, no_atoms, strong_meas, classical)

    mle_variance = zeros(1, length(g_array));

    %w_array = 10.^(log10(pi/T):0.025:log10(pi/tau));

    if classical

        sample_fun = @sample_classical_trajectory;
        max_likelihood_fcn = @maximum_likelihood_est_classical;

    else

        sample_fun = @sample_trajectory;
        max_likelihood_fcn = @maximum_likelihood_est_fast;

    end
    
    for g_i=1:length(g_array)

        g = g_array(g_i);
    
        % sample random trajectories
        [outcomes, ~] = sample_fun(T, tau, g, w, no_atoms, no_samples, ...
            strong_meas);

        [w_est] = max_likelihood_fcn(outcomes, T, tau, g, ...
            no_atoms, no_samples, strong_meas);

        mle_variance(g_i) = sum((w - w_est).^2/no_samples, "all");
    
    end

end