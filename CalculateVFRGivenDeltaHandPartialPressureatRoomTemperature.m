% Given Oven Temperature, Partial pressure at room temperature, DeltaH
% Calculate using the Clausius-Clapiron the corresponding partial pressures
% at the different given Oven temperature, and then calculate what is the
% final diameter and hence VFR (Dfinal^3/Dinit^3) as a function of
% temperature.
%clear all;clc
a=[];
d=1;
% constants
R       = 8.3144621;    % Gas constant [J/(mol*K)]

% inputs
PPRT    =1e-4;      % Partial pressure at room temperature [Pascal]


T       = linspace(298,500,200);
DeltaH  = linspace(100e3,80e3,200);    
%DeltaH = 100e3 ;% Enthalpy [J/mol]
T_amb   = 298;          % Room temperature [K]
Di      = 50;        % Initial diameter in nm

% Calculate corresponding partial pressure in oven temperatures according
% the Clausius-Clapiron equation (double checked correctly !)
p = exp(-(DeltaH./(R*T))+(DeltaH/(R*T_amb)+log(PPRT)));


Diameterfinal=[];
options = optimset('Display', 'off');
for i=1:length(T)
    f = @(Df)DebugginingEnthalpyFromExperiment(T(i),p(i), Di, Df);
    DiameterTemp = fzero(f,80,options);
    Diameterfinal = [Diameterfinal DiameterTemp];
    VFRtemp = DiameterTemp.^3/(Di^3);
    if (VFRtemp < 0.02)
        T = T(1:length(Diameterfinal));
        break
    end
    %Diameterfinal = [Diameterfinal fminbnd(f,0,Di,options)];
end


figure(2)
hold all
plot(T, Diameterfinal.^3/(Di^3),'-')
xlabel('Oven Temperature [K]');
ylabel('VFR');





