%%% Calculation for std of ORR 
% calculates the mean and std of ORR with varying lambda_s (total number of
% photon counts)
% generates heat map showing SNR in dB scale

clear 
clc
close all
% ls = 1:30; % lambda1 + lambda2

l1 = 5;
l2 = 5;

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

% add std and mu to the arrays
% values of l1, l2 are also indices!
std = sqrt(res);
mu = l1 / (l1 + l2);

% calculate SNR
snr = mu / std;

% convert matrix to dB
snr = 10*log10(snr);

fprintf("Calculated SNR: %.3f dB\n", snr);