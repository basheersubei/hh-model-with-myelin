%% main.m
% where parameters are configured
% performs myelinated HH simulation (runs doSim)
% ~~~~~~~~~ BEGIN: ~~~~~~~~~ 

clear all;
close all;

global factor Io myelin_curve_freq myelin_thickness timeDesired pulseStart pulseWidth N Vrest GK_max GNa_max GL ENa EK EL n0 m0 h0 axon_radius axon_length Cm Cm_myelin Ri Rout Rm_myelin membrane_thickness;


GK_max   = 36e-3;   %S/cm^2  (maximum specific potassium conductance)
GNa_max  = 120e-3;  %S/cm^2  (maximum specific sodium conductance)
GL       = 0.3e-3;  %S/cm^2  (specific leakage conductance)
ENa = 55e-3;   %volts
EK  = -75e-3;  %volts
EL  = -33e-3;  %volts

tMax = 0.1;
pulseWidth = 0.04;
pulseStart = 0.04;
number_of_time_intervals= 1e3;
timeDesired = linspace(0,tMax,number_of_time_intervals); %desired time vector for ode integrator



%% modify parameters here based on problem
Io = 4e-7;
%Io = 0;
%Io = 0.02;

N = 30;
Vrest = -60e-3; %in volts

n0 = 0.3534;  %initial n-gate probability of being in open state
m0 = 0.06915; %initial m-gate probability of being in open state
h0 = 0.5142;  %initial h-gate probability of being in open state

myelin_curve_freq = 3;
axon_radius = 1e-2; %in cm
%axon_radius = 1e-3;
%axon_radius = 1;
%total_myelin_volume = 6.6e-5;

factor = 1.0;

%myelin_thickness = 1e-6;
axon_length = 1; %in cm
%axon_length = 100;
membrane_thickness = 5e-7; % in cm
Cm = 1e-6; % specific capacitance of membrane in F/cm^2
Cm_myelin = 5e-9 * factor;
Ri = 5; % specific intercellular resistance in ohm*cm
Rout = 1e5; % specific resistance exiting the tube, also called Ro in ohm*cm
Rm_myelin = 4e5 / factor; % specific myelin resistance in ohm*cm^2

%myelin_thicknesses =  4 * (2 * pi * axon_radius) .* factors; %in cm;
myelin_thickness = axon_radius * factor;
%% perform simulation with required parameters
doSim();

%sound notification that simulation is done
load chirp 
sound(y,Fs)
