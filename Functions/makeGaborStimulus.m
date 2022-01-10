function gaborPatch = makeGaborStimulus(pixelsPerCM,sizeGrating,phaseGrating,spatialFrequency,annulus,distFromScreen,environment,whichScreen)
% spatialFrequency: cycles/degree



width = degrees2pixels(sizeGrating, distFromScreen, pixelsPerCM, whichScreen);

if annulus == 1
    inner_degree = 1;
else
    inner_degree =0;
end
size_inner_most_circle = degrees2pixels(inner_degree, distFromScreen, [], whichScreen);
start_linear_decay_in_degree = .5; 
start_linear_decay = degrees2pixels(start_linear_decay_in_degree, distFromScreen, [], whichScreen);

nCycles = sizeGrating*spatialFrequency; % number of cycles in a stimulus

% compute the pixels per grating period
pixelsPerGratingPeriod = width / nCycles;

spatialFrequency = 1 / pixelsPerGratingPeriod; % How many periods/cycles are there in a pixel?
radiansPerPixel = spatialFrequency * (2 * pi); % = (periods per pixel) * (2 pi radians per period)

%gray = 127;

% adjust contrast
% black = 0;
% absoluteDifferenceBetweenBlackAndGray = abs(gray - black);
%[background, Lmin_rgb, Lmax_rgb] = calibrate_lum(contrastGrating, environment, scanner);
%virtual_background = (Lmin_rgb+Lmax_rgb)/2;


halfWidthOfGrid = width / 2;
widthArray = (-halfWidthOfGrid) : halfWidthOfGrid;  % widthArray is used in creating the meshgrid.

% Creates a two-dimensional square grid.  For each element i = i(x0, y0) of
% the grid, x = x(x0, y0) corresponds to the x-coordinate of element "i"
% and y = y(x0, y0) corresponds to the y-coordinate of element "i"
[x y] = meshgrid(widthArray);

% Creates a sinusoidal grating, where the period of the sinusoid is 
% approximately equal to "pixelsPerGratingPeriod" pixels.
% Note that each entry of gratingMatrix varies between minus one and
% one; -1 <= gratingMatrix(x0, y0)  <= 1

% the grating is oriented vertically unless otherwise specified.

stimulusMatrix = sin(radiansPerPixel * x + phaseGrating);
% Because the difference between Lmax_rgb and background is not equal to
% the difference between Lmin_rgb and background (the differences are equal
% in luminance, but not in rgb), correct the range of the positive values
% of stimulusMatrix.

%rgb_range_up = Lmax_rgb - background;
rgb_range_up = background;
%rgb_range_down = background - Lmin_rgb;
rgb_range_down = background;

rgb_range_up_down_ratio = rgb_range_up/rgb_range_down;
stimulusMatrix(stimulusMatrix > 0) = stimulusMatrix(stimulusMatrix > 0) * rgb_range_up_down_ratio;

% Make a fading annulus, to use as a mask.
annulusMatrix = makeLinearMaskCircleAnn(width+1,width+1,size_inner_most_circle,start_linear_decay,width/2);
stimulusMatrix = stimulusMatrix .* annulusMatrix;

%Make the grating
gaborPatch = background + rgb_range_down * stimulusMatrix;
end

% c = contrastGrating
% background
% Lmax_rgb
% Lmin_rgb
% Lmax = max(max(gaborPatch))
% Lmin = min(min(gaborPatch))
% Lmean = mean(mean(gaborPatch))
