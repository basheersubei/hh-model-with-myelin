%% MyelinDyn.m
% dyn function to be thrown in ode23s
% returns the characteristic equation of a HH myelinated axon
% uses HH for unmyelinated parts and passive cable for myelinated parts
% parameters: time t (used to calculate drive or vector b), and current
% states [Vm; n; m; h] <-- these are all column vectors of length N.
% ~~~~~~~~~ BEGIN: ~~~~~~~~~ 
function s_dot = MyelinDyn(t, s)

%% import global variables
 global N c_m g_i g_out g_m Io pulseWidth pulseStart GK_max GNa_max GL ENa EK EL myelination;
 
%% extract state variables
Vm = s(1:N);
n = s(N+1:2*N);
m = s(2*N+1:3*N);
h = s(3*N+1:4*N);

%% now construct G matrix (in preparation for VmDot diff eq)

%first find G of all the channels
G_Na = GNa_max .* (m .^3) .* h;
G_K = GK_max .* n .^4;

% GL is static and already given
G_total = G_Na + G_K + GL;

%check if the node has myelin
%get correct value for membrane conductance
G_total(myelination>0) = g_m(myelination>0); %where there is myelination, set to g_m manually and without using gate conductances.



%build vector for gi diagonal
gi_diagonal = zeros(1, N-1); %notice length of this vector is N-1
for i = 1:length(gi_diagonal)
   gi_diagonal(i) = g_i;  
end

%build vector for middle diagonal (diagonal index 0)
middle_diagonal = zeros(1,N); %note that middle_diagonal is one element bigger than gi_diagonal
for i = 1:length(middle_diagonal)
    if(i==1)
       middle_diagonal(i) = -G_total(i) - g_i;
    elseif(i==N)
        middle_diagonal(i) = -G_total(i) -g_i - g_out;
    else
        middle_diagonal(i) = -G_total(i) - (2*g_i);
    end
end


% Construct G matrix from 3 diagonal matrices.
G = diag(gi_diagonal, -1) + diag(middle_diagonal, 0) + diag(gi_diagonal, 1);

%% now build b vector (drive)
drive = ENa .* G_Na + EK .* G_K + EL * GL;

%make sure to zero out drive where there is myelin (since there are no
%gates there)
drive(myelination>0) = 0;

% add input current pulse at the right times
if(t> pulseStart && t< (pulseWidth+pulseStart)) 
    drive(1) = drive(1) + Io/2;
end

%% first, solve V_m diff equation using given V_m, G's, and E's.
VmDot = (1./c_m) .* ((G * Vm) + drive);

%% then, solve the diff equations for n, m, and h using alphas and betas.

coefficients = getProbCoeff(Vm); %returns all coefficients as columns of length N each

alpha_n = coefficients(:, 1);
beta_n = coefficients(:, 2);
alpha_m = coefficients(:, 3);
beta_m = coefficients(:, 4);
alpha_h = coefficients(:, 5);
beta_h = coefficients(:, 6);

n_dot = -(alpha_n + beta_n) .* n + alpha_n;
m_dot = -(alpha_m + beta_m) .* m + alpha_m;
h_dot = -(alpha_h + beta_h) .* h + alpha_h;

%% repackage statesDot from characteristic equations
s_dot = [VmDot; n_dot; m_dot; h_dot];

end