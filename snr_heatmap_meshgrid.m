%%% Calculation for std of ORR 
% calculates the mean and std of ORR with varying lambda_s (total number of
% photon counts)
% generates heat map showing SNR in dB scale

clear 
clc
close all
% ls = 1:30; % lambda1 + lambda2

% Init max values for lambda1 (FAD) and lambda2 (NADH)
l1m = 30;
l2m = 30; 

% rows are l2 (vertical), columns are l1 (horizontal)
m_var = zeros(l2m, l1m);
m_mu  = zeros(l2m, l1m);

% for ii = 1:length(ls)
for l2 = 1:l2m
    for l1 = 1:l1m 
        ls = l1 + l2;

        tol = 1e-16; % MATLAB computational precision
        Kmax = 1e10; % arbitrary large value, break once value equals zero
        
        tot = 0;
        nzero = 0;
        
        %%% calculate summation using poisson 
        for k = 1:Kmax
            
            temp = 1/k * poisspdf( k, ls );
            if (temp < tol && nzero), break, end
            if (temp > tol), nzero = 1; end
        
            tot = tot + temp;
    
        end
    
        % multiply summation by the scalar quantities to get rest of value
        
        % NOTE since we use poisspdf() above, that term includes the exp(-ls) term
        % that would've been in the following line
        res = ( 1 / (1-exp(-ls)) ) * ( l1*l2/ls^2 ) * tot;
    
        % add var and mu to the arrays
        % values of l1, l2 are also indices!
        m_var(l2, l1) = res;
        m_mu(l2, l1)  = l1 / (l1 + l2);
    end
end

% calculate std from variance
m_std = sqrt( m_var );

% calculate SNR
m_snr = m_mu ./ m_std;

% convert matrix to dB
m_mu  = 10*log10(m_mu);
m_var = 10*log10(m_var);
m_snr = 10*log10(m_snr);

%%% PLOT MEAN (SIGNAL) %%%
figure, imagesc(flipud(m_mu)) % excludes edge cases where l1, l2 = 0
yt = get(gca, 'YTick');
set(gca, 'YTickLabel', fliplr([0 yt ]) -5)

xlabel("FAD Photon Counts")
ylabel("NAD(P)H Photon Counts")

ax = gca;
ax.FontSize = 18;
ax.FontWeight = 'bold';

c = colorbar;
c.Label.String = "Mean (dB)";
c.FontSize = 18;
c.FontWeight = 'bold';
% clim([0 15])


%%% PLOT STANDARD DEVIATION (NOISE) %%%
figure, imagesc(flipud(m_var)) % excludes edge cases where l1, l2 = 0
hold on, plot(0.5:30.5, fliplr(0.5:30.5), 'Color', 'r', 'LineWidth', 4);
yt = get(gca, 'YTick');
set(gca, 'YTickLabel', fliplr([0 yt ]) -5)

xlabel("FAD Photon Counts")
ylabel("NAD(P)H Photon Counts")

ax = gca;
ax.FontSize = 18;
ax.FontWeight = 'bold';

c = colorbar;
c.Label.String = "Variance (dB)";
c.FontSize = 18;
c.FontWeight = 'bold';
% clim([0 15])


%%% PLOT SNR %%%
figure, imagesc(flipud(m_snr)) % excludes edge cases where l1, l2 = 0
yt = get(gca, 'YTick');
set(gca, 'YTickLabel', fliplr([0 yt ]) -5)

xlabel("FAD Photon Counts")
ylabel("NAD(P)H Photon Counts")

ax = gca;
ax.FontSize = 18;
ax.FontWeight = 'bold';

c = colorbar;
c.Label.String = "SNR (dB)";
c.FontSize = 18;
c.FontWeight = 'bold';
clim([0 15])

%%% Plot surface contour

fad_lim = l1m;
nadh_lim = l2m;
% create meshgrids
[FAD, NADH] = meshgrid( 1:l1m, 1:l2m );

figure, surf( FAD, NADH, m_snr )
xlabel("FAD Photon Counts")
ylabel("NAD(P)H Photon Counts")
zlabel("SNR (dB)")
