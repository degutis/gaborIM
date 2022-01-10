function [width,vertDist,horzDist] = visualAngleCalc(visualAngle, eccentricity, stimulusAngle, distanceScreen, ignoreVertOffset)

%values for testing
% visualAngle = 2.34;
% eccentricity = 6.4;
% stimulusAngle = 29;
% distanceScreen = 40;

%% Explanation of calculation
% d = distance from eye to fixation dot
% E = eccentricity
% V = visual angle
% Ecc_hyp = hypotenuse distance from fixation cross to the center of the square

% We know that the distance from the eye to the center of the stimulus can
% be determined by the eccentricity in this way:
% cos(E) = d/distnaceEyetoCenterSquare
% tan(E) = Ecc_hyp/d 
% Note, that the other end of this triangle is the diagonal and NOT the
% horizontal or vertical distance to the stimulus square.

%Moving onto the visual angle: 
% The visual angle calculation should take into account that the stimulus
% is both off-center horizontally and off-center vertically, thus, just like in the
% eccentricity case, the diagonal distance from center to stimulus should
% the basis of the calculation. 

% Let's define a variable a --> 1/2 the hypotenuse of the square stimulus,
% meaning (2a)^2 = width_stimulus^2 + width_stimulus^2 (since width=height)
% Let's set up a few equations with the things that we know:
% as mentioned: tan(E) = Ecc_hyp/d so d*tan(E)=Ecc_hyp

% Let's define two more variables Theta1 and Theta2 in such a way:
% d*tan(Theta1) = Ecc_hyp-a (so this would measure the distance from eyes
% to the closer edge of the stimulus)
% AND d*tan(Theta2) = Ecc_hyp+a (distance from eyes to the farther edge of
% stimulus) 
% We also know that: V = Theta2-Theta1 based on the way we defined them

% So let's set up some equations basing them on the eccentricity measurement:
% tan(Theta1) = (Ecc_hyp-a)/d
% tan(Theta2) = (Ecc_hyp+a)/d
% since V=Theta2-Theta1 --> Theta2 = V+Theta1
% tan(V+Theta1) = (Ecc_hyp+a)/d

% The rest is just algebra plugging in tan(Theta1) = (Ecc_hyp-a)/d into the
% other equation and solving for a. 
% And knowing that tan(A+B) = (tan(A)+tan(B))/(1-tan(A)tan(B))

%% Calculation
%Eccentricity based calculation
Ecc_hyp = tan(degtorad(eccentricity))*distanceScreen;
vertDist = Ecc_hyp*sin(degtorad(stimulusAngle));
horzDist = Ecc_hyp*cos(degtorad(stimulusAngle));

if ignoreVertOffset == 1 
    % current assumption ignores the vertical offset of the stimulus
    % instead of using diagonal distance to the image, we will use the
    % horizontal distance (thus ignoring the vertical component). 
    % As such, the width of the stimulus becomes 2*a
    Ecc_hyp = horzDist;
    multipleOfA = 2; 
else
    % In this case the width of the stimulus is sqrt(2)*a following the
    % rules of a 45/45/90 degree triangle. 2*a in this case is the
    % hypotenuse of the triangle, and we need the side. 
    multipleOfA = sqrt(2);
end

%Visual Angle based calculation

f1 = (distanceScreen^2)*((tan(degtorad(visualAngle))^2)+1);
f2 = (Ecc_hyp^2)*(tan(degtorad(visualAngle))^2);
a = (sqrt(f1+f2)-distanceScreen)/tan(degtorad(visualAngle));

width = multipleOfA*a; 


end
