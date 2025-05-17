
% This function performs MLE over w_array given some outcomes array.

% Outcomes: (T/tau, no_atoms, no_samples)
% prob_w shape: (1, no_samples)

function [w_estimate] = maximum_likelihood_est_fast(outcomes, T, tau, ...
    g, no_atoms, no_samples, strong_meas)

    %w_array = 10.^(log10(pi/T):0.04:log10(pi/(2*tau)));
    w_array = pi/T:(pi/(2*tau)-pi/T)/137:pi/(2*tau);
    %disp(length(w_array));

    % tic

    [~, ind] = max(squeeze(sum(simulate_given_trajectory_2(outcomes, ...
             w_array, T, tau, g, no_atoms, no_samples, strong_meas), 1)), [], 2);

    w_0 = w_array(ind);

    % toc

    w_estimate = zeros([1 , no_samples]);

    options = optimset('TolX',1e-6);

    % tic

    parfor i = 1:no_samples

        fun = @(x)-1*sum(simulate_given_trajectory_2(outcomes(:, :, i), x, ...
            T, tau, g, no_atoms, 1, strong_meas), "all");

        x0 = w_0(i);

        w_estimate(i) = fminsearch(fun, x0, options);

    end

    % toc


end
