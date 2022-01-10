function [background, Lmin_rgb, Lmax_rgb] = calibrate_lum2(contrast)
% some RGB & corresponding luminance levels from the dummy scanner screen,
% with brightness (+30%) and contrast (+0%) toned down and the cardboard
% screen in place to prevent high levels of background lighting.
Lmin_rgb = NaN;
Lmax_rgb = NaN;
rgb = [0        30      60      90      120     150     180     190     210     220     230     240     250     255];

lum(1,:) = [0.18 3.23 11.815 18.7 31.5 50.7 71.3 77.05 90.55 96.9 103 116 135.5 133];
lum = mean(lum,1);

% interpolate over the whole rgb range
lum = interp1(rgb,lum,0:255,'spline');

%middle of the luminance range:
medium_lum = (min(lum) + max(lum))/2;
background = 256;

%what colour should the background be? medium luminance, or darker?
while lum(background) > medium_lum%/2
    background = background - 1;
end

if exist('contrast', 'var')
    Lmin = medium_lum - (medium_lum-min(lum))*contrast;
    Lmax = medium_lum + (max(lum)-medium_lum)*contrast;
    Lmin_rgb = 1;
    while lum(Lmin_rgb) < Lmin
        Lmin_rgb = Lmin_rgb + 1;
    end
    Lmax_rgb = 255;
    while lum(Lmax_rgb) > Lmax
        Lmax_rgb = Lmax_rgb - 1;
    end
    
end
