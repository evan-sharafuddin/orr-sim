%% ORR Imaging Model
% Evan Sharafuddin
% Sorrells Lab
% WashU

% TODO need to resolve issues with the rounding
% TODO need to think of image quality metric

clear
clc
close all

%%% Set Model Params %%%
pc = 5; % photon count

%%% Load Image %%%
MAXPIXEL = 255; 
img = imread("images/phantom.png", "png");
img_norm = double( img ) / MAXPIXEL;

PRINT_FIGS = 1;
VERBOSE_PRINT = 0;

%%% Calculate ORR %%%
ORR = orr_model( img_norm, pc, PRINT_FIGS, VERBOSE_PRINT );

%%

%%% Quantify Image Quality %%%
%{
SOURCES:
https://books.google.com/books?id=VZvqqaQ5DvoC&pg=PA280#v=onepage&q&f=false
https://en.wikipedia.org/wiki/Signal-to-noise_ratio#Definition
https://scientificimaging.com/knowledge-base/signal-and-noise/
https://journals.asm.org/doi/10.1128/aem.02536-07

DEFINITIONS: 
either mean/std or mean^2/std^2
20log10 or 10log10?
%}

DB_PREFIX = 10; % TODO need to determine if its 10 or 20

% (1) Per pixel SNR using multiple image calculations
num_iter = 10;

ORRs = ORR;

for i = 1:num_iter-1
    ORRs = cat( 3, ORRs, orr_model(img_norm, pc) ); % concatenate along third dimension
end

mu_ORRs = mean(ORRs, 3);
sigma_ORRs = std(ORRs, 0, 3);

% NOTE this calculation will get you inf's (x/0) and Nan's (0/0)
SNR_pixelwise = mu_ORRs ./ sigma_ORRs;

% anything Nan cooresponds to a signal value of 0, meaning we can set SNR
% to 0
SNR_pixelwise( isnan(SNR_pixelwise) ) = 0;

% for visualization, set infinite SNR to maximum noninfinite value
% NOTE this is not necissary if using imagesc()
% SNR_pixelwise( isinf(SNR_pixelwise) ) = max(max(SNR_pixelwise(~isinf(SNR_pixelwise))));

% convert to dB
dB_SNR_pixelwise = DB_PREFIX * log10(SNR_pixelwise);

% plot colormap
figure, imagesc(dB_SNR_pixelwise)
title( sprintf("(%d PC, %d iter) Pixelwise SNR of Simulated ORR image", pc, num_iter) )
c = colorbar;
c.Label.String = "Pixelwise SNR values (dB)";
c.Label.FontWeight = "bold";
caxis([-10 10])
% xlim([150 450]), ylim([150 450])
set(gcf, 'Position',  [100, 100, 300, 500]*2.5)
saveas(gcf, sprintf('SNR_px_%d.png', pc), 'png')

% (2) Per pixel error based on reference image
% NOTE currently having issues with edge cases driving up the error values
error = abs( (ORR - img_norm) ./ img_norm );

% both Nan and inf are due to division by zero, which only occurs when
% pixel values are around zero
error(isnan(error) | isinf(error)) = 0;

% figure, imagesc(error);
% title( sprintf("(%d PC, %d iter) Pixelwise error values for simulated ORR image", pc, num_iter) )
% c = colorbar;
% % NOTE get spikes in values but these are only when ref image has a value
% % of nearly zero (i.e., at a boarder) and we don't really care
% caxis([0 1]); % set to 2 as max value on the colorbar
% c.Label.String = "Pixelwise error values from ref image";
% c.Label.FontWeight = "bold";
% xlim([150 450]), ylim([150 450])
% set(gcf, 'Position',  [100, 100, 300, 500]*2.5)
% saveas(gcf, sprintf('error_px_%d.png', pc), 'png')

figure, imagesc(error);
title( sprintf("(%d PC, %d iter) Pixelwise error values for simulated ORR image", pc, num_iter), ...
    'FontSize', 30, 'FontWeight', 'bold');
c = colorbar;
% NOTE get spikes in values but these are only when ref image has a value
% of nearly zero (i.e., at a boarder) and we don't really care
caxis([0 1]); % set to 2 as max value on the colorbar
c.Label.String = "Pixelwise Error";
c.Label.FontSize = 35;
c.Label.FontWeight = "bold";
% xlim([150 450]), ylim([150 450])
set(gca, 'FontSize', 20); % Sets font size for axes numberings and labels
set(gcf, 'Position', [100, 100, 300, 500]*2.5)
ax = gca;
ax.PositionConstraint = "outerposition";
saveas(gcf, sprintf('error_px_%d.png', pc), 'png')



% (3) Mean squared error based on reference image
% TODO not sure if this is a completely valid metric
rmse_ORR = sqrt( sum(sum( (ORR - img_norm).^2 ))/numel(ORR) );
fprintf("RMSE error metric for %d photon counts: %.5f\n", pc, rmse_ORR)



