% Given Oven Temperature, Partial pressure at room temperature, DeltaH
% Calculate using the Clausius-Clapiron the corresponding partial pressures
% at the different given Oven temperature, and then calculate what is the
% final diameter and hence VFR (Dfinal^3/Dinit^3) as a function of
% temperature.
% clear all;clc
function VFRResult = VFR(a, b, T)
%load data.mat
% constants
R       = 8.3144621;    % Gas constant [J/(mol*K)]

% inputs
PPRT    = 2.77e-5;      % Partial pressure at room temperature [Pascal]
%T       = data(:,1);    % Oven temperature [K]
%T       = linspace(298,400,400);
%DeltaH  = linspace(100e3,50e3,10);        % Enthalpy [J/mol]
DeltaH  = a*T+b;        % Enthalpy [J/mol]
T_amb   = 298;          % Room temperature [K]
Di      = 101.8;        % Initial diameter in nm

% Calculate corresponding partial pressure in oven temperatures according
% the Clausius-Clapiron equation (double checked correctly !)
p = exp(-(DeltaH./(R*T))+(DeltaH/(R*T_amb)+log(PPRT)));

tic
Diameterfinal=[];
options = optimset('Display', 'off');
for i=1:length(T)
    f = @(Df)DebugginingEnthalpyFromExperiment(T(i),p(i), Di, Df);
    Diameterfinal = [Diameterfinal fzero(f,80,options)];
    %Diameterfinal = [Diameterfinal fminbnd(f,0,Di,options)];
end
toc

VFRResult = Diameterfinal.^3/(Di^3);
%figure(2)
%hold off
%plot(T, Diameterfinal.^3/(Di^3),'*')
%xlabel('Oven Temperature [K]');
%ylabel('VFR');