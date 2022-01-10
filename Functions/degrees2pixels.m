function [nPixels, nPixelsUnrounded] = degrees2pixels(degrees, distFromScreen_inCm, pixels_perCm)
% [nPixels, nPixelsUnrounded] = degrees2pixels(degrees, distFromScreen_inCm, pixels_perCm)
%
% Converts degrees to pixels, given distance from screen in centimeters
% and pixels per centimeter. Result is rounded by default. Use
% nPixelsUnrounded to get the unrounded result of the calculation.
%
% Default value for distance from screen is 60 cm.
% If not specified, pixels per cm is measured using Screen('Resolution', 
% screenNumber) and Screen('DisplaySize', screenNumber). Default for
% screenNumber is max(Screen('Screens')).

if ~exist('distFromScreen_inCm','var') || isempty(distFromScreen_inCm)
    distFromScreen_inCm = 60;
end

%pixels_perCm = 27.49; 

% convert degrees to centimeters
sizeOnScreen_inCm = 2 * distFromScreen_inCm * tan((degrees/2) * (pi/180));

% convert centimeters to pixels
nPixelsUnrounded = sizeOnScreen_inCm * pixels_perCm;
nPixels = round(nPixelsUnrounded);
end