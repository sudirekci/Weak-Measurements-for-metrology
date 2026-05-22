
% This function simulates a given trajectory, given some outcomes and
% omega. NO BACKACTION

% prob_trajectory shape: (no_atoms, no_samples)
% prob_w_given_tr shape: (no_atoms, no_samples)
% outcomes shape: (T/tau, no_atoms, no_samples)

function [prob_w_given_tr] = simulate_classical_trajectory(outcomes, ...
    w_init, T, tau, g, no_atoms, no_samples, strong_meas, p_e)

    N = floor(T/tau);

    pplus = 1/2*(1+(1-2*p_e)*sin(2*(cat(1, g*ones(1,N-1)', strong_meas*pi/4 + ...
        (1-strong_meas)*g))).*cos(2*w_init*tau.*(1:N)'));

    prob_w_given_tr = squeeze(sum(log(pplus).*outcomes + log(1-pplus).*...
        (no_atoms-outcomes), 1));

end
