function draw_fixationdot(window,background, width,height,dotsize,color,colorInside,distFromScreen,pixelsPerCM)
if nargin <=4
    dotsize = [0.2 0.05];
end
if nargin <=5
    color = 0;
end
if nargin<=6
    colorInside = color;
end

x = width/2;
y = height/2;

% Screen('DrawDots', cfg.win, [cfg.width/2, cfg.height/2], dotSize, color);


d1 = dotsize(1); % diameter of outer circle (degrees)
d2 = dotsize(2); % diameter of inner circle (degrees)
[~,ppd] = degrees2pixels(1,distFromScreen,pixelsPerCM);

Screen('FillOval', window, color, [x-d1/2 * ppd, y-d1/2 * ppd, x+d1/2 * ppd, y+d1/2 * ppd], d1 * ppd);
Screen('DrawLine', window, background, x-d1/2 * ppd, y, x+d1/2 * ppd, y, min(d2 * ppd,5));
Screen('DrawLine', window, background, x, y-d1/2 * ppd, x, y+d1/2 * ppd, min(d2 * ppd,5));
Screen('FillOval', window, colorInside, [x-d2/2 * ppd, y-d2/2 * ppd, x+d2/2 * ppd, y+d2/2 * ppd], d2 * ppd);
end