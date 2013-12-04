%% Given VFR and temperature, find the best fit paramters that connect enthalpy and temperature
% for a given matarial

% Given Oven Temperature, Partial pressure at room temperature, DeltaH
% Calculate using the Clausius-Clapiron the corresponding partial pressures
% at the different given Oven temperature, and then calculate what is the
% final diameter and hence VFR (Dfinal^3/Dinit^3) as a function of
% temperature.
clear all;clc

load data.mat

% Plot Measured Data of VFR as a function of T
hold off
plot(MeasuredData(:,1)+273.15,MeasuredData(:,2),'*')

% Non linear fit data
% let a and b be:
% Don't fit on the part that is close to zero, Set the threshold for zero
Threshold = 0.06;
T = MeasuredData(MeasuredData(:,2)>Threshold,1)+273.15;
VFRm = MeasuredData(MeasuredData(:,2)>Threshold,2)';


%Coef=[-4.9e2, 2.8e5];
Chi2 = @(Coef) sum((VFR(Coef(1),Coef(2),T)-VFRm).^2);

coef = fminsearch(Chi2, [-4.9e2 2.8e5])

%% Plot Result of best fit

subplot(2,1,1)
hold off
%VFRfit = VFR(-69.618110996500008,1.398873976775031e5,T);
VFRfit = VFR(coef(1),coef(2),T);
plot(T, VFRfit) 
hold all
plot(MeasuredData(:,1)+273.15,MeasuredData(:,2),'*')
xlabel('Temp (K)')
ylabel('VFR')
subplot(2,1,2)
plot(T, coef(1)*T+coef(2))
xlabel('Temp (K)')
ylabel('Enthalpy (J/mol)')

