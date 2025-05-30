
close all

noSamples = 1e3;


N = 64;

delta_phi_arr = 10.^(log10(0.1):0.2:log10(20));

rosenband_bmse = zeros([1, length(delta_phi_arr)]);

tic

for l=1:length(delta_phi_arr)

    delta_phi = delta_phi_arr(l);
    %no_ensembles = ceil(delta_phi./0.7)
    no_ensembles = 2;
    speeds = 2.^((no_ensembles-1):-1:0);
    %no_atoms = ceil(N*2^(no_ensembles-1)/(2^(no_ensembles)-1))*ones([1, no_ensembles]);
    no_atoms

    no_atoms = floor(N/no_ensembles)*ones([1, no_ensembles]);
    no_atoms(end) = no_atoms(end) + N - sum(no_atoms);

    no_atoms = [10, 54];

    rosenband_bmse(l) = rosenband_sim(delta_phi, speeds, no_atoms, noSamples);

end

toc

figure
hold on
plot(delta_phi_arr, sqrt(rosenband_bmse)./delta_phi_arr*sqrt(N), LineWidth=3)
plot(delta_phi_arr, 1./(delta_phi_arr), LineWidth=3)
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
ax = gca;
ax.FontSize = 20; 
ax.FontName = "times";
xlabel("\delta\phi")
ylabel("N\Delta\phi/\delta\phi")





function [error_variance] = rosenband_sim(delta_phi, speeds, no_atoms, noSamples)

    samples = normrnd(0, delta_phi, [noSamples, 1]);

    for l = 1:length(speeds)

        %num_atoms = N*speeds(l);
        num_atoms = no_atoms(l);

        if num_atoms == 0
            continue
        end

        r1 = binornd(round(num_atoms), (1+sin(samples./speeds(l)))/2)./round(num_atoms);
        %r2 = binornd(round(num_atoms/2), (1+sin(samples./speeds(l)))/2)./round(num_atoms/2);

        beta = asin(2*r1-1);

        %beta = angle(r1-1/2+1j*(r2-1/2));

        samples = samples - beta*speeds(l);

    end

    error_variance = sum(samples.^2)/noSamples;

end