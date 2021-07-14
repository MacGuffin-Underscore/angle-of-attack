%% Branch 1: Rhett A Smith | 10.29.2020
% Collaborators: Navin Solomon
% Changes:
%   1) Changed the base code into a function that can be called for a range of values for both the 
%      thickness ratio and angle of attack.
%   2) Changed the name of variables that do not work within functions 
%       \\M and beta are both functions)
%   3) Changed sections to include freestream [inf] as section 1 \\personal preference)
%   4) Changed everything to work in deg, got sick of working in rad
%   5) Added seperation of questions and displays the solutions
%   6) Added a loop to calculate for varying tc and alpha values
%
% Known Problems:
%   1) govern() is still broken, working with govern2
%   2) Program is slow, taking over 4 seconds to run

%% ASE 4343 Compressible Aerodynamics
% angleOfAttack.m
% November 9, 2020

% Initialize
clc; clearvars; format shortg; tic
global gam
gam = 1.4;

%% 1) Test Case: alphadeg = 5, tc = 0.088
[Cl_test, Cd_test] = ClCd(5, 0.088);
disp(['1) For the test case of alpha = 5 and tc = 0.088'])
disp(['   Cl = ',num2str(Cl_test)])
disp(['   Cd = ',num2str(Cd_test)])
disp(' ')

%% 2-3) alpha of 0-10deg and tc of [0.044 , 0.088 , 0.176]
% Initialize alpha and thickness ratio
alpha = [0:1:10];
tc = [0.044 , 0.088 , 0.176];

% Loop
for i = 1:length(tc)
    for j = 1:length(alpha)
        [Cl(i,j) , Cd(i,j)] = ClCd(alpha(j) , tc(i));
        C_ld(i,j) = Cl(i,j)/Cd(i,j);
    end
end

%% Plots
% Cl vs alpha
figure(1); hold on; title('Cl vs alpha')
plot(alpha(:),Cl(1,:))
plot(alpha(:),Cl(2,:))
plot(alpha(:),Cl(3,:))
legend(['tc = ',num2str(tc(1))],['tc = ',num2str(tc(2))],['tc = ',num2str(tc(3))],'Location','northwest')
xlabel('alpha(deg)') ; ylabel('Cl')
hold off

% Cd vs alpha
figure(2); hold on; title('Cd vs alpha')
plot(alpha(:),Cd(1,:))
plot(alpha(:),Cd(2,:))
plot(alpha(:),Cd(3,:))
legend(['tc = ',num2str(tc(1))],['tc = ',num2str(tc(2))],['tc = ',num2str(tc(3))],'Location','northwest')
xlabel('alpha(deg)') ; ylabel('Cd')
hold off

% Cl/Cd vs alpha
figure(3); hold on; title('Cl/Cd vs alpha')
plot(alpha(:),C_ld(1,:))
plot(alpha(:),C_ld(2,:))
plot(alpha(:),C_ld(3,:))
legend(['tc = ',num2str(tc(1))],['tc = ',num2str(tc(2))],['tc = ',num2str(tc(3))],'Location','northeast')
xlabel('alpha(deg)') ; ylabel('Cl/Cd')
hold off

disp('4) I give my partner a grade of: ENTER NUMBER HERE')
disp(' ')
disp(['This program took ',num2str(toc*1000),'ms to run'])


%% Functions
function [Cl , Cd] = ClCd(alpha , tc)
%% Function Information
% ------------------------------------
% Finds the Cl and Cd of a diamond airfoil
% ------------------------------------
% Inputs:
%   alphadeg: Angle of attack in degrees
%   tc: Thickness ratio of a diamond airfoil
%
% Outputs:
%   Cl: Coefficent of lift
%   Cd: Coefficent of drag
% -------------------------------------

%% Known Values
% airfoil properties \\unused?)
len = 1;                                % length of section
xc = 0.5;                               % thickness location

% angles
eta = round(atand(tc));                 % rounded foil half angle in deg

% air properties
M_inf = 2;                              % Free-stream mach number
gam = 1.4;

%% Pressure Symbolics Setup

syms p1 p2 p3 p4 p5

%% In-Line Functions
% Pressure Ratio from Mach Solver

Pratio = @(Ma,Mb) ((1 + ((gam-1)/2) * Ma^2)/(1 + ((gam-1)/2) * Mb^2))^(gam/(gam-1));


% Normal to Mach Solver

NormalM = @(Mn) sqrt((1 + ((gam-1)/2) * Mn^2)/(gam*Mn^2 - (gam-1)/2));


%% Section 2 \\Open to freestream)
theta_2 = eta - alpha;
if theta_2 >=0
    [M2 , P2_inf] = obliqueShock(theta_2 , M_inf);
else
    theta_2 = abs(theta_2);
    [M2 , P2_inf] = expansionFan(theta_2 , M_inf);
end

p2 = P2_inf * p1;

%% Section 3 \\Connected to section 2)
theta_3 = 2 * eta;

[~ , P32] = expansionFan(theta_3 , M2);

p3 = P32 * P2_inf * p1;

%% Section 4 \\Open to freestream)
theta_4 = eta + alpha;
if theta_4 >=0
    [M4 , P4_inf] = obliqueShock(theta_4 , M_inf);
else
    theta_4 = abs(theta_4);
    [M4 , P4_inf] = expansionFan(theta_4 , M_inf);
end

p4 = P4_inf*p1;

%% Section 5 \\Connected to section 4)
theta_5 = 2 * eta;

[~ , P54] = expansionFan(theta_5 , M4);

p5 = P54 * P4_inf * p1;

%% Lift and Drag Coefficient Calculation

syms l c

l = (0.5*c)/(cosd(eta));

L = (p4 - p3)*(l*cosd(eta+alpha)) + (p5 - p2)*(l*cosd(eta-alpha));

D = (p4 - p3)*(l*sind(eta+alpha)) + (p2 - p5)*(l*sind(eta-alpha));

Cl = double(L / ((gam/2)*p1*M_inf^2*c));

Cd = double(D / ((gam/2)*p1*M_inf^2*c));
end

function [M2 , P21] = obliqueShock(th , M1)
%% Function Information
% ------------------------------------
% Finds the PRESSURE RATIO and MACH NUMBER Across a Oblique Shock
% ------------------------------------
% Inputs:
%   th: Theta in degrees
%   M1: Mach number of section 1
%
% Outputs:
%   P21: Pressure ratio section 2/1
%   M2: Mach number of section 2
% -------------------------------------
global gam
if th == 0
    M2 = M1;
    P21 = 1;
else
    be = thetaBetaM(th , M1);

    Mn1 = M1*sind(be);                                          % Normal mach number at section 1
    Mn2 = sqrt((1+((gam-1)/2)*Mn1^2)/(gam*Mn1^2-((gam-1)/2)));  % Normal mach number at section 2

% OUTPUT CALCULATIONS

    M2 = Mn2/sind(be-th);
    P21 = 1+((2*gam)/(gam+1))*(Mn1^2-1);
end
end

function [M2 , P21] = expansionFan(th , M1)
%% Function Information
% ------------------------------------
% Finds the PRESSURE RATIO and MACH NUMBER Across a Expansion Fan
% ------------------------------------
% Inputs:
%   th: Theta in degrees
%   M1: Mach number of section 1
%
% Outputs:
%   P21: Pressure ratio section 2/1
%   M2: Mach number of section 2
% -------------------------------------
global gam
if th == 0
    M2 = M1;
    P21 = 1;
else
    PM = @(M) sqrt((gam+1)/(gam-1))*atand(sqrt(((gam-1)/(gam+1))*(M^2-1)))-atand(sqrt(M^2-1)); 
%                                                                       \\Prandtl Meyer function)
    nu1 = PM(M1);               % Nu at mach of section 1
    nu2 = th + nu1;             % Nu at mach of section 2
    M2_fun = @(M) PM(M)-nu2;    % Mach of section 2

% OUTPUT CALCULATIONS
    opts = optimset('Diagnostics','off', 'Display','off');
    M2 = fsolve(M2_fun,M1,opts);     % Solves the P-M equation for M2 \\with an initial guess of M1)
    P21 = ((1+((gam-1)/2)*M1^2)/(1+((gam-1)/2)*M2^2))^(gam/(gam-1));
end
end

function be = thetaBetaM(th , M)
%% Function Information
% ------------------------------------
% Finds the Value of BETA for a given THETA and MACH NUMBER
% ------------------------------------
% Inputs:
%   th: Theta in degrees
%   M: Mach number
%
% Outputs:
%   be: Beta in degrees
% -------------------------------------
global gam

% govern1 = @(be) atand(2*cotd(be) * (M^2 * sind(be)^2 - 1)/(M^2 * (gam + cosd(2*be) + 2))) - th;
govern2 = @(be) tand(be) * ( ((gam+1)*M^2) / (2 * (M^2 * sind(be)^2 -1)) -1) - cotd(th);
%   \\Bad news, govern 1 is broke)

% OUTPUT CALCULATIONS
opts = optimset('Diagnostics','off', 'Display','off');
be = fsolve(govern2,45,opts);

end
