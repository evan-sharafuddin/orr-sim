
function ORR = orr_model( ...
    img_norm, ... Normalized "ideal image" that will be used to drive the ORR simulator
    pc,  ... Sum of photon counts in NADH and FAD channels (assumes same for each pixel)
    DISP_FIGURES, ... Flag to display figures
    VERBOSE ... Verbose output
)

ERR = -1;

if nargin < 4
    VERBOSE = 0;
end

if nargin < 3
    DISP_FIGURES = 0; % default behavior is to not print figures
end

if nargin < 2
    pc = 50; % default photon count value 
end 

if nargin < 1
    fprintf("Insufficent amount of arguments, aborting\n")
    ORR = ERR;
    return
end

%%% Calculate FAD and NAD(P)H Intensities %%%
% Assume that FAD+NAD(P)H = photon_count

FAD = img_norm * pc;
NADH = pc - FAD;

% round values
FAD  = round( FAD );
NADH = round( NADH );

%%% Calculate Shot Noise %%%
% ASSUMPTION: treat phantom image as expected value for photon count per 
% pixel
FAD_sn  = poissrnd(FAD);
NADH_sn = poissrnd(NADH);


%%% Calculate Gaussian Noise %%%

%{

ASSUMPTION: we can model effects external to our system as AGWN
* Assume detector noise is primary source of external noise
    * 100 cps is a typical figure for this
* typical image size is 390x390 pixels; typically acquisition time is 120
s; typical dwell time is 10^-3 s
--> model noise as a binary image of random dark counts (each pixel is
either 0 (no dark counts) or 1 (single dark count)). Assume that multiple
dark counts at a given pixel is unlikely
--> sum(dark_counts) = (img_dim)^2 * (100 cps) * (10^-3 s)


REFERENCES:
https://doi.org/10.1117/3.725073
https://en.wikipedia.org/wiki/Additive_white_Gaussian_noise

%}

%%% param

%%
img_width = size(img_norm, 1);
img_height = size(img_norm, 2);
t_dwell = 1e-3; % dwell time per pixel [s]
detector_cps = 100; % typical detector noise [dark counts/s]
num_dark_counts = img_width * img_height * t_dwell * detector_cps;

% change to generate appropriate dark count binary image
% NOTE: want there to be a num_dark_counts/img_width/img_height percent chance
% that a given pixel is 1 
%   we are generating values following a guassian distribution, and then will
%   assign to zero or one
%   can use area under the pdf to design the gaussian from which we draw
%   values from

mu = 0; sigma = 1; % choose normal distribution

tol = 50; % how close to the num_dark_counts value our binary image has to be

chance = num_dark_counts/img_width/img_height;
prob_is_greater_x = @(x) 1 - cdf('Normal', x, mu, sigma);
fun = @(x) prob_is_greater_x(x) - chance; 
cutoff = fzero( fun, 0 );

%%% create binary image of dark counts
% loop until tolerance requirement is reached
iter1 = 0;
iter2 = 0;


while true
    iter1 = iter1 + 1;
    vals = normrnd( mu, sigma, [(img_width*img_height), 1]);
    vals_bin = vals >= cutoff;
    if abs( sum(vals_bin) - num_dark_counts ) <= tol
        vals_bin_FAD = vals_bin;
        break;
    end
end

while true
    iter2 = iter2 + 1;
    vals = normrnd( mu, sigma, [(img_width*img_height), 1]);
    vals_bin = vals >= cutoff;
    if abs( sum(vals_bin) - num_dark_counts ) <= tol        
        vals_bin_NADH = vals_bin;
        break;
    end
end

if VERBOSE 
    fprintf("Took %d iterations to find dark count image for FAD\n", iter1)
    fprintf("Took %d iterations to find dark count image for NADH\n", iter2)
end

% create binary image
dark_counts = reshape(vals_bin, img_width, img_height);
% figure, imshow(dark_counts)

%%% Add Together Noise %%%
FAD_dc = reshape(vals_bin_FAD, img_width, img_height);
NADH_dc = reshape(vals_bin_NADH, img_width, img_height);

FAD_ct = FAD_sn + FAD_dc;
NADH_ct = NADH_sn + NADH_dc;

% zero out all negative values (not physically possible)
% FAD_ct(FAD_ct < 0) = 0;
% NADH_ct(NADH_ct < 0) = 0;


%%% Calculate ORR %%%
ORR = double( FAD_ct ./ ( FAD_ct + NADH_ct ) );
% replace NaN vals with 0
ORR( isnan(ORR) ) = 0;


%%% Display figures %%%

if DISP_FIGURES 
    
    %%% Display Figures Individually
    figure
    imshow(img_norm)
    title("Ground Truth Image", 'FontWeight', 'bold', ...
          'FontSize', 35)
    % xlim([150 450]), ylim([150 450])
    set(gcf, 'Position',  [100, 100, 300, 500]*2.5)
    ax = gca;
    ax.PositionConstraint = "outerposition";
    saveas(gcf, 'phantom_og.png', 'png')
    
    figure
    imshow(ORR)
    title(sprintf("Simulated ORR Image | %d PC", pc), 'FontWeight', 'bold', ...
          'FontSize', 35)
    % xlim([150 450]), ylim([150 450])
    set(gcf, 'Position',  [100, 100, 300, 500]*2.5)
    saveas(gcf, sprintf('phantom_%d.png', pc), 'png')


    % %%% Display intermediaries 
    % figure
    % ha = tight_subplot(1, 2, 0.05, 0.05, 0.05);
    % axes(ha(1)), imshow( FAD/max(max(FAD)) )
    % title("FAD (no noise)")
    % axes(ha(2)), imshow( NADH/max(max(NADH)) )
    % title("NADH (no noise)")
    % 
    % figure
    % ha = tight_subplot(1, 2, 0.05, 0.05, 0.05);
    % axes(ha(1)), imshow( FAD_sn/max(max(FAD_sn)) )
    % title("FAD (with Poisson shot noise)")
    % axes(ha(2)), imshow( NADH_sn/max(max(NADH_sn)) )
    % title("NADH (with Poisson shot noise)")
    % 
    % figure
    % ha = tight_subplot(1, 2, 0.05, 0.05, 0.05);
    % axes(ha(1)), imshow( (FAD + FAD_dc) / max(max( FAD + FAD_dc )) )
    % title("FAD (with dark count noise)")
    % axes(ha(2)), imshow( (NADH + NADH_dc) / max(max( NADH + NADH_dc )) )
    % title("NADH (with dark count noise)")
    % 
    % 
    % % Display original & final images
    % figure
    % ha = tight_subplot(1, 2, 0.05, 0.05, 0.05);
    % axes(ha(1)), imshow(img_norm)
    % title("Original Phantom Image (normalized)")
    % axes(ha(2)), imshow(ORR)
    % title("Simulated ORR Image")
end 

end 


%{ 

imshow() NUANCES
* if arg is a double(), then image will be displayed on scale of [0, 1]. As
far as positive values go, anything over 1 will be treated as a fully white
pixel
*   this is important to keep in mind whether we want to normalize wrt
photon_count, or wrt the maximum value in the matrix (which can be bigger
than photon_count due to noise)

* NOTE: the below only applies to FAD and NADH intermediary image
generations. The final ORR image is already normalized between 0 and 1 via
the formula used to calculate the values
* RUNNING ASSUMPTION: normalize image based on maximum value, as image
generation will be based on this instead of a given value of photon counts
(which might be hard to determine accurately in practice)
*   impact: this will reduce dynamic range, especially if there is a pixel
that has a really high value due to additive noise

%}