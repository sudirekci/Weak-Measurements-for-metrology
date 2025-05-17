
% This function samples one trajectory, returns the probability and
% outcomes. NO BACKACTION

% outcomes shape: (T/tau, no_samples)
% prob_trajectory shape: (no_atoms, no_samples)

function [outcomes, prob_trajectory] = sample_classical_trajectory(T, tau, ...
    g, w, no_atoms, no_samples, strong_meas)

    N = floor(T/tau);

    pplus = 1/2*(1+sin(2*(cat(1, g*ones(1,N-1)', strong_meas*pi/4 + ...
        (1-strong_meas)*g))).*cos(2*w*tau.*(1:N)'));

    outcomes = squeeze(sum(ceil(pplus-rand(N, no_atoms, no_samples)), 2));

    prob_trajectory = squeeze(sum(log(pplus).*outcomes + log(1-pplus).*...
        (no_atoms-outcomes), 1));

end