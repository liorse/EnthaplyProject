% Input : T - Oven Temperature [K]
%       : PP - Partial pressure at Oven Temperature [Pascal] Input for
%       : Df - Final Diameter [nm]
%       solver
function ZeroResult = DebugginingEnthalpyFromExperiment(T, PP, Di, Df)
%clear all; clc
% input: [Temperature (k) particle final sizes (nm)]
% output: Partial pressure(pascal) in oven temperature for a specific particle made of
% single type of molecule


% Calculation according to paper J. Phys. Chem A. Vol 114. No. 13. 2010
% Inputs for experiment are
% T  - temperature of the thermal denuder
%T               = [40.22+273.15];           % Temperature in which the experiment was perfromed. will become an array
% Di - Initial particle diameter
Di = Di * 1e-9; % diameter in m
% Df - Final particle diameter
%Df = 95.3e-9; % diameter in m

% Output
% Saturation vapor pressure
% dH - Delta H = Latent heat of sublimation
% partial pressure in Pascal

% Physical parameters for containing molecular crystal
% adipic methylamine salt (2:1)
MolMass         = 0.2083;       % Molecular weight [kg/mol]
Density         = 1154.5;        % Molecular crystal density [kg/m^3]
Gamma_i           = 0.1;        % ------------- [J/m^2]
Sigmaii         = 7.22e-10;       % ------------- [m]
Epsilonii_kb    = 646.3; %708.21;        % Epsilon/Kb    [K]
Alpha           = 1;             % accomodation coef. is set to 1 in this calculation

% Physical parameters for containing molecular crystal
% Addipic Acid for an example
% MolMass         = 0.14614;       % Molecular weight [kg/mol]
% Density         = 1.36e3;        % Molecular crystal density [kg/m^3]
% Gamma_i           = 0.06;        % ------------- [J/m^2]
% Sigmaii         = 6.3e-10;       % ------------- [m]
% Epsilonii_kb    = 818.27; %708.21;        % Epsilon/Kb    [K]
% Alpha           = 1;             % accomodation coef. is set to 1 in this calculation

AirMass         = 0.02897;       % Air Mass [kg/mol]
MeltingPoint    = 426;           % Melting point [K]
SigmaAir        = 3.617e-10;     % Interparticle distance in which the potential is zero
R               = 8.3144621;     % Gas constant [J/(mol*K)]
EpsilonAir_kb   = 97.0;          % Potential well depth
TotalPressure   = 101325;        % 1 atm pressure in pascal
Kb              = 1.38065e-23;   % Boltzmann constant J/K
Diameter_intvar = 20;            % Temperary assignment for an diameter integration variable
Nd              = 218e-12;       % Nitrogen molecular diameter in m
L               = 0.6;           % Tube length (m)
r               = (25.4e-3)/8;   % Tube Radius (m)
Q               = 300*(1e-6)/60; % Flow rate   (m^3/sec)
T_ambient       = 298.15;        % Room Temperature (K)



Df=Df*1e-9; % T is in Kelvin, Df is in meters

Ufront = 2*Q/(pi*r^2)*(T/T_ambient); % molecule velocity inside thermal denuder
dt = L/Ufront; % molecule residence time inside the thermal denuder in seconds.

% Calculation according to paper J. Phys. Chem A. Vol 114. No. 13. 2010
Sigma_i_air = (Sigmaii+SigmaAir)/2;  % Interparticle distance in which the potential is zero
Epsilon_i_air = sqrt(Epsilonii_kb*EpsilonAir_kb); % Potential well depth
T_star = T/Epsilon_i_air;   % kb*T/epsilon

% This formula for Omega_i_air is coming from page 866, Transport phenomena second edition by Bird
Omega_i_air = 1.06036./(T_star.^0.1561)+0.193./(exp(0.47635*T_star))+1.03587./(exp(1.52996*T_star))+1.76474./exp(3.89411*T_star);

% The following is taken from:
% D. J. Rader , P. H. McMurry & S. Smith (1987): Evaporation Rates of
% Monodisperse Organic Aerosols in the 0.02- to 0.2-?m-Diameter Range, Aerosol Science and
% Technology, 6:3, 247-260
D_i_air = (sqrt(6.022141e23) * 2*sqrt(2*pi)/(16*pi*Omega_i_air))*(sqrt((Kb*T).^3*(1/MolMass+1/AirMass)))/(TotalPressure*Sigma_i_air^2);

% Avi's versoin for D_i_air
%D_i_air = 3*sqrt(2*pi)*sqrt(Kb^3*T^3*(1/MolMass+1/AirMass)*6.022141e23)/(16*pi*Omega_i_air*TotalPressure*Sigma_i_air^2)

% Molecule A mean velocity
C_A = sqrt(8*Kb*T/(pi*(MolMass/6.022141e23)));

p_0 = (Density*R*T./(4*D_i_air*dt*MolMass));

% Integrand
MFP = @(D_AB, C_A) 3*D_AB/C_A;
Knudsen1 = @(Lamdai, D_p) 2*Lamdai./D_p;
f_number1 = @(Kn_i,Alpha) (1 + Kn_i)./(1+0.3773.*Kn_i+1.33.*Kn_i.*(1+Kn_i)/Alpha);

fun = @(D) (D./f_number1(Knudsen1(MFP(D_i_air,C_A), D), Alpha)).*exp(-4*Gamma_i*MolMass./(D.*Density*R*T));

% Compelete experssion for partial pressure -
% the - PP is used to find the diameter in an outside script
ZeroResult = -p_0*integral(fun, Di, Df) -  PP;
%ZeroResult = abs(p_0*integral(fun, Di, Df) +  PP);


end