%% getProbCoeff.m
% returns alpha and beta parameters for gate models given a voltage
% parameter V_m voltage
% returns an array of the alpha and beta coefficients
% V_m is to be taken in Volts (returns in 1/s)
function [coeffOut] = getProbCoeff(V_m)

global Vrest;
%convert V_m and Vrest to milliVolts
V = (V_m - Vrest ) .*1000;

% remember, each of these variables is a column vector of size V_m
% these variables are in units 1/ms
alpha_n = (0.1 - 0.01.*V) ./ (exp( 1 - 0.1 .*V) - 1 );
beta_n = 0.125 .* exp(-V ./ 80);
alpha_m = (2.5 -0.1 .*V ) ./ (exp( 2.5 - 0.1 .*V) - 1);
beta_m = 4 .* exp(-V ./ 18);
alpha_h = 0.07 .* exp(-V ./ 20);
beta_h = 1 ./ (exp(3 - 0.1 .*V) + 1);

coeffOut = [alpha_n beta_n alpha_m beta_m alpha_h beta_h] .* 1000; % returns coeff values in units 1/s
end