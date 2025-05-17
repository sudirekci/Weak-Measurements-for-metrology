
% This function samples one trajectory, returns the probability and
% outcomes.

% outcomes shape: (T/tau, no_atoms, no_samples)
% prob_trajectory shape: (no_atoms, no_samples)

function [outcomes, prob_trajectory] = sample_trajectory(T, tau, g, w, ...
    no_atoms, no_samples, strong_meas)

    N = floor(T/tau);

    outcomes = zeros(N, no_atoms, no_samples, "int8");
    prob_trajectory = ones(1, no_atoms, no_samples)*log(1/2)*N;

    %theta = zeros(1, no_atoms, no_samples);

    w_exp = exp(-1i*2*tau*w);
    theta_prime = w_exp;

    for j=1:N

        if strong_meas && j == N

            g = pi/4;

        end

        pplus = 1/2*(1+sin(2*g).*real(theta_prime));

        kplusind = 2*logical(ceil(pplus-rand(1, no_atoms, no_samples)))-1;

        prob_trajectory = prob_trajectory + log(kplusind.*(2*pplus-1)+1);

        denom = (1+kplusind.*sin(2*g).*real(theta_prime));

        theta_prime = ((real(theta_prime) + ...
            1i*cos(2*g)*imag(theta_prime)+kplusind.*sin(2*g))./denom).*w_exp;
       
        outcomes(j, :, :) = kplusind;

    end

    %prob_trajectory = reshape(prob_trajectory, [no_atoms, no_samples]);
    prob_trajectory = squeeze(prob_trajectory);

end
