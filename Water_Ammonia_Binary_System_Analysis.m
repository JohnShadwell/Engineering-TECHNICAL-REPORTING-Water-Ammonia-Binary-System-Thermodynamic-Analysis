%% Step 1: Calculate the value of a12 for the binary mixture

clear;
clc;
close all;

%Constants
R = 0.08314; %uses liters and bar

TcA = 405.6; %phys prop data pulled from Table A.1.2 for ammonia
PcA = 112.77;
wA  = 0.250;

TcW = 647.3; %same data for water
PcW = 220.48;
wW  = 0.344;

T = linspace(250, 500, 1000); %creates 1000 temperature points to represent a12(T)

mu12 = -0.34918; %empirical data from a 2012 study (included in written report) relating to the binary mixing parameter k12
nu12 = 2.4467e-4;

kappaA = 0.37464 + 1.54226*wA - 0.26992*wA^2; %PR EOS kappa vals to plug into function for alpha
kappaW = 0.37464 + 1.54226*wW - 0.26992*wW^2;


aA  = zeros(1,length(T)); %preallocating arrays for temperature dependent vals
aW  = zeros(1,length(T));
a12 = zeros(1,length(T));
a12_withoutk = zeros(1,length(T));
k12 = zeros(1,length(T));


for i = 1:length(T) %iterate through temperatures

    k12(i) = mu12 + nu12 * T(i); %binary mixing parameter k12 as an empirically derived function of T (same 2012 study)

    alphaA = (1 + kappaA*(1 - sqrt(T(i)/TcA)))^2; %alpha values from PR EOS
    alphaW = (1 + kappaW*(1 - sqrt(T(i)/TcW)))^2;

    aA(i) = 0.45724 * (R^2 * TcA^2 / PcA) * alphaA; %a values for pure species from PR EOS
    aW(i) = 0.45724 * (R^2 * TcW^2 / PcW) * alphaW;

    a12(i) = (1 - k12(i)) * sqrt(aA(i) * aW(i)); %mixing parameter from quadratic mixing 
    a12_withoutk(i) = sqrt(aA(i) * aW(i));

end

% Plots a12 as a function of T
figure;
plot(T, a12);
hold on;
plot(T, a12_withoutk);
xlabel('Temperature (K)');
ylabel('a12 (L^2*bar/mol^2)');
title('Mixing Parameter a12 vs Temperature');
legend('Including k12', 'Without k12');




%% Step 2a: Fugacity coefficient as a function of pressure

T = 350; %constant conditions
yA = 0.5;
yW = 0.5;

Pvar = linspace(1, 100, 1000); %1000 pressures from 1 to 100 bar

phiA = zeros(size(Pvar));   % NH3 fugacity coefficient
phiW = zeros(size(Pvar));   % H2O fugacity coefficient

for i = 1:length(Pvar) %iterate through pressures
    P = Pvar(i); %given pressure
    [phiA(i), phiW(i)] = phiSolve(T, P, yA, yW, R, TcA, PcA, wA, TcW, PcW, wW, mu12, nu12); %solve for fugacity cofficients at given pressure
end

figure; %Plots fugacity coefficients at each pressure
plot(Pvar, phiA);
hold on;
plot(Pvar, phiW);
xlabel('Pressure (bar)');
ylabel('Fugacity Coefficient');
title('Fugacity Coefficients vs Pressure (T = 350 K, yA = 0.5)');
legend('NH3', 'H2O');


%% Step 2b
P  = 20; %constant conditions
yA = 0.5;
yW = 0.5;

Tvar = linspace(250, 450, 1000); %1000 temperatures from 280 to 450 K

phiA = zeros(size(Tvar));   % NH3 fugacity coefficient vs T
phiW = zeros(size(Tvar));   % H2O fugacity coefficient vs T

for i = 1:length(Tvar) %iterates through temperatures
    T = Tvar(i); %given temperature
    [phiA(i), phiW(i)] = phiSolve(T, P, yA, yW, R, TcA, PcA, wA, TcW, PcW, wW, mu12, nu12); %calculates fugacity coefficients at given temperature
end

figure; %plots fugacity coefficients at each temperature
plot(Tvar, phiA);
hold on;
plot(Tvar, phiW);
xlabel('Temperature (K)');
ylabel('Fugacity Coefficient');
title(sprintf('Fugacity Coefficients vs Temperature (P = 20 bar, yA = 0.5)'));
legend('NH3', 'H2O');


%% Step 2c: Fugacity coefficients as a function of composition

T = 350; %constant conditions
P = 20;

y1var = linspace(0, 1, 1000);   %1000 points from 0 mol NH3/mol to 1 mol NH3/mol
phiA = zeros(size(y1var));
phiW = zeros(size(y1var));

for i = 1:length(y1var) %iterate through compositions
    yA = y1var(i); %given NH3 composition
    yW = 1 - yA; %H2O composition
    [phiA(i), phiW(i)] = phiSolve(T, P, yA, yW, R, TcA, PcA, wA, TcW, PcW, wW, mu12, nu12); %computes fugacity coefficients based on composition
end

figure; %plots fugacity coefficients at each composition
plot(y1var, phiA); 
hold on;
plot(y1var, phiW);
xlabel('yNH3');
ylabel('Fugacity Coefficient');
title(sprintf('Fugacity Coefficients vs Composition (T = 350 K, P = 20 bar)'));
legend('NH3','H2O');


function [phiA, phiW] = phiSolve(T, P, yA, yW, R, TcA, PcA, wA, TcW, PcW, wW, mu12, nu12) %computes fugacity coeffeficients for a binary mixture with the given parameters

    kappaA = 0.37464 + 1.54226*wA - 0.26992*wA^2; %pure species kappa vals from PR EOS
    kappaW = 0.37464 + 1.54226*wW - 0.26992*wW^2;

    alphaA = (1 + kappaA*(1 - sqrt(T/TcA)))^2; %pure species alpha vals from PR EOS
    alphaW = (1 + kappaW*(1 - sqrt(T/TcW)))^2;

    aA = 0.45724 * (R^2 * TcA^2 / PcA) * alphaA; %pure species a parameter from PR EOS
    aW = 0.45724 * (R^2 * TcW^2 / PcW) * alphaW;

    bA = 0.07780 * (R * TcA / PcA); %pure species b parameter from PR EOS
    bW = 0.07780 * (R * TcW / PcW);

    k12 = mu12 + nu12 * T; %binary mixing parameter k12 from 2012 study

    a12 = (1 - k12) * sqrt(aA * aW); %mixing parameter a12 from quadratic mixing

    am = yA^2 * aA + 2*yA*yW*a12 + yW^2 * aW; %a parameter from quadratic mixing for PR EOS

    bm = yA*bA + yW*bW; %b parameter from quadratic mixing for PR EOS

    A = am * P / (R^2 * T^2); %A value for cubic form PR EOS (needed to find Z for fugacity coefficient formula)
    B = bm * P / (R * T); %B value for cubic form PR EOS
    c3 = 1; %z^3 coefficient for cubic PR EOS
    c2 = -(1 - B); %z^2 coefficient for cubic PR EOS
    c1 = (A - 3*B^2 - 2*B); %z coefficient for cubic PR EOS
    c0 = -(A*B - B^2 - B^3); %constant for cubic PR EOS

    Zroots = roots([c3 c2 c1 c0]); %finds all roots of PR EOS with cubic form = 0
    Zreal  = Zroots(abs(imag(Zroots)) < 1e-12); %checks if the imaginary part of each root = 0 to determine whether the root is real or imaginary, and then returns only the real roots
    Z = max(real(Zreal)); %takes the largest real root, which is the Z for vapor phase

    %applies the equation for fugacity coeffecient for PR EOS
    phiA = exp((bA/bm)*(Z - 1) - log(Z - B) + (am/(2*sqrt(2)*bm*R*T)) * ((bA/bm) - 2*(yA*aA + yW*a12)/am) * log((Z + (1+sqrt(2))*B) / (Z + (1-sqrt(2))*B)));
    phiW = exp((bW/bm)*(Z - 1) - log(Z - B) + (am/(2*sqrt(2)*bm*R*T)) * ((bW/bm) - 2*(yA*a12 + yW*aW)/am) * log((Z + (1+sqrt(2))*B) / (Z + (1-sqrt(2))*B)));
end


%% Step 3: Using Generalized Correlations from Table

T = 500; %temperature (constant)
Pideal = 5; %ideal pressure
Pnonideal = 50; %nonideal pressure

%calls interpolater function with given table values (found manually with
%Tr and Tc values in Appendix C)
phiAideal = tableSolver(T/TcA, Pideal/PcA, wA, [-0.002 -0.004, -0.002 -0.003], [0.000 0.001 0.001 0.001], [1.2 1.3], [0.025 0.05]); 
phiAnonideal = tableSolver(T/TcA, Pnonideal/PcA, wA, [-0.021 -0.042, -0.016 -0.033], [0.004 0.009 0.007 0.014], [1.2 1.3], [0.25 0.5]);
phiWideal = tableSolver(T/TcW, Pideal/PcW, wW, [-0.003 -0.009, -0.003 -0.007], [-0.003 -0.007, -0.002 -0.005], [0.75 0.8], [0.01 0.025]);
phiWnonideal = tableSolver(T/TcW, Pnonideal/PcW, wW, [-0.035 -0.237, -0.029 -0.076], [-0.030 -0.677 -0.020 -0.056], [0.75 0.8], [0.1 0.25]);

%prints results
fprintf('Ideal Coefficients:\nWater = %.5f\nAmmonia = %.5f\n\nNonideal Coefficients:\nWater = %.5f\nAmmonia = %.5f\n\n',phiWideal,phiAideal,phiWnonideal,phiAnonideal);

function phi = tableSolver(Tr, Pr, w, phi0vals, phi1vals, Trvals, Prvals) %calculates fugacity coefficient with given conditions and table values using generalized correlations

Prdist = (Pr - Prvals(1)) / (Prvals(2) - Prvals(1)); %used for double interpolation
Trdist = (Tr - Trvals(1)) / (Trvals(2) - Trvals(1));

tempInter1 = phi0vals(1) + (phi0vals(2)-phi0vals(1)) * Prdist; %interpolate with pressures for phi0
tempInter2 = phi0vals(3) + (phi0vals(4)-phi0vals(3)) * Prdist;
logphi0 = tempInter1 + (tempInter2-tempInter1) * Trdist; %interpolate with temperature for phi0

tempInter1 = phi1vals(1) + (phi1vals(2)-phi1vals(1)) * Prdist; %interpolate with pressures for phi1
tempInter2 = phi1vals(3) + (phi1vals(4)-phi1vals(3)) * Prdist;
logphi1 = tempInter1 + (tempInter2-tempInter1) * Trdist; %interpolate with temperature for phi1

phi = 10^(logphi0 + w * logphi1); %calculates phi
end


%% Step 4: Accuracy of Lewis Fugacity Rule
yA = 0.5; %redefine composition
yW = 0.5;

[phiAActualIdeal, phiWActualIdeal] = phiSolve(T, Pideal, yA, yW, R, TcA, PcA, wA, TcW, PcW, wW, mu12, nu12); %finds phi values from PR EOS from step 2
[phiAActualNonideal, phiWActualNonideal] = phiSolve(T, Pnonideal, yA, yW, R, TcA, PcA, wA, TcW, PcW, wW, mu12, nu12);

fAActualIdeal = yA * phiAActualIdeal * Pideal; %finds fugacity values from PR EOS from step 2
fWActualIdeal = yW * phiWActualIdeal * Pideal;
fAActualNonideal = yA * phiAActualNonideal * Pnonideal;
fWActualNonideal = yW * phiWActualNonideal * Pnonideal;

fALewisIdeal = yA * phiAideal * Pideal; %finds fugacity values from generalized correlations from step 3 and using lewis rule
fWLewisIdeal = yW * phiWideal * Pideal;
fALewisNonideal = yA * phiAnonideal * Pideal;
fWLewisNonideal = yW * phiWnonideal * Pideal;

errorAIdeal = abs(((fALewisIdeal - fAActualIdeal)/fAActualIdeal) * 100); %calculates percent error values for each case
errorWIdeal = abs(((fWLewisIdeal - fWActualIdeal)/fWActualIdeal) * 100);
errorANonideal = abs(((fALewisNonideal - fAActualNonideal)/fAActualNonideal) * 100);
errorWNonideal = abs(((fWLewisNonideal - fWActualNonideal)/fWActualNonideal) * 100);

%prints results
fprintf("\nIdeal Values (T = 500 K, P = 5 bar)\n");
fprintf("Ammonia:\n");
fprintf("PR mixture fugacity = %.5f bar\n", fAActualIdeal);
fprintf("Lewis-rule fugacity = %.5f bar\n", fALewisIdeal);
fprintf("%% Error = %.5f %%\n\n", errorAIdeal);

fprintf("Water:\n");
fprintf("PR mixture fugacity = %.5f bar\n", fWActualIdeal);
fprintf("Lewis-rule fugacity = %.5f bar\n", fWLewisIdeal);
fprintf("%% Error = %.5f %%\n", errorWIdeal);

fprintf("\nNon-Ideal Values (T = 500 K, P = 50 bar)\n");
fprintf("Ammonia:\n");
fprintf("PR mixture fugacity = %.5f bar\n", fAActualNonideal);
fprintf("Lewis-rule fugacity = %.5f bar\n", fALewisNonideal);
fprintf("%% Error = %.5f %%\n\n", errorANonideal);

fprintf("Water:\n");
fprintf("PR mixture fugacity = %.5f bar\n", fWActualNonideal);
fprintf("Lewis-rule fugacity = %.5f bar\n", fWLewisNonideal);
fprintf("%% Error = %.5f %%\n\n", errorWNonideal);