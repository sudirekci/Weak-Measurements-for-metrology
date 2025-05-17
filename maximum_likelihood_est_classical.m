
% This function performs MLE over w_array given some outcomes array.

% Outcomes: (T/tau, no_samples)
% prob_w shape: (1, no_samples)

function [w_estimate] = maximum_likelihood_est_classical(outcomes, T, tau, ...
    g, no_atoms, no_samples, strong_meas)


    w_array = pi/T:0.15*pi/T:pi/(2*tau);

    tic
    prob_w = zeros(length(w_array), no_samples);

    ii = 1;

    for w_init = w_array

        [prob_w_given_tr] = simulate_classical_trajectory(outcomes, ...
            w_init, T, tau, g, no_atoms, no_samples, strong_meas);

        prob_w(ii, :) = sum(prob_w_given_tr, 1);
        ii = ii+1;

    end

    % figure
    % plot(w_array, mean(prob_w, 2))


    [~, ind] = max(prob_w, [], 1);
    w_0 = w_array(ind);

    w_estimate = zeros([1 , no_samples]);

    parfor i = 1:no_samples

        fun = @(x)-1*simulate_classical_trajectory(outcomes(:, i), x, ...
            T, tau, g, no_atoms, 1, strong_meas);

        x0 = w_0(i);

        w_estimate(i) = fminsearch(fun, x0);

    end
 
end