
% This function simulates a given trajectory, given some outcomes and
% omega.

% prob_w_given_tr shape: (no_atoms, no_samples, length(w_init))
% outcomes shape: (T/tau, no_atoms, no_samples)

function [prob_w_given_tr] = simulate_given_trajectory_2(outcomes, ...
    w_init, T, tau, g, no_atoms, no_samples, strong_meas)


    N = floor(T/tau);

    %theta = ones(1, no_atoms, no_samples, length(w_init));
    prob_w_given_tr = ones(1, no_atoms, no_samples, length(w_init))*log(1/2)*N;

    cos_2g = cos(2*g);
    sin_2g = sin(2*g);

    w_init = exp(-1i*2*tau*reshape(w_init, [1,1,1,length(w_init)]));
    theta_prime = w_init;


    for j=1:N

        if strong_meas && j == N

            g = pi/4;

            cos_2g = 0;
            sin_2g = 1;

        end

        kplusind = sin_2g*double(outcomes(j, :, :));

        %theta_prime = theta.*w_init;

        p = (1+kplusind.*real(theta_prime));

        % the factor of 1/2 is already accounted for in the initialization
        prob_w_given_tr = prob_w_given_tr + log(p);

        %theta = (cos_2g^2*theta_prime + sin_2g^2*conj(theta_prime)+kplusind)./p;
        theta_prime = ((real(theta_prime) + 1i*cos_2g*imag(theta_prime)+kplusind)./p).*w_init;

    end

    %disp(abs(theta_prime));
    %prob_w_given_tr = reshape(prob_w_given_tr, [no_atoms, no_samples, length(w_init)]);

    prob_w_given_tr = squeeze(prob_w_given_tr);

end

