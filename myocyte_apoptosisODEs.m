function [dydt, algvars] = myocyte_apoptosisODEs(t,y,params)
% BME 6315: myocte apoptosis
% Vignesh Valaboju

% Assign names for parameter values and state variables
% [k1f,k1r,k2f,k2r,kcat,Km,k4,k5,Vratio,L] = params{:};
[k1f, k1r, k2f, k2r, k3, km, k4, k5, k6, k7, k8f, k8r, k9, k10, k11, k12] = params{:};

% R = y(1);
% LR = y(2);
% E = y(3);
% LRE = y(4);
% S = y(5);
% P = y(6);
% Pnuc = y(7);
ESR = y(1); 
FADD = y(2); 
ESR_FADD = y(3); 
procas8 = y(4);
procas8_ESR_FADD = y(5); 
E = y(6);
cas8 = y(7);
cas9 = y(8); 
cytc = y(9);
bax_bak = y(10); 
bax = y(11);
bak = y(12); 
bcl = y(13);
ISR = y(14); 
NICD = y(15); 
JaggedNotch = y(16); 
cas367 = y(17);

% Reaction rates
% react1 = k1f*L*R - k1r*LR;      % [uM/s] L+R<->LR
% react2 = k2f*LR*E - k2r*LRE;    % [uM/s] LR+E<->LRE
% react3 = kcat*LRE*S/(Km+S);     % [uM/s] LRE catalyzes S->P
% react4 = k4*P;                  % [uM/s] P->S
% react5 = k5*(P-Pnuc);           % [uM/s] Pnuc <-> P
J1 = k1f*ESR*FADD - k1r*ESR_FADD;
J2 = k2f*procas8*ESR_FADD - k2r*procas8_ESR_FADD;
J3 = k3*procas8_ESR_FADD*E/(km + E);
J4 = k4*cas8;
J5 = k5*cas9;
J6 = k6*cytc;
J7 = k7*bax_bak;
J8 = k8f*bax*bak* - k8r*bax_bak;
J9 = k9*bcl; 
J10 = k10*ISR;
J11 = k11*NICD;
J12 = k12*JaggedNotch;

% Differential equations;
% dR = -react1;             % [uM/s] free receptor
% dLR = react1-react2;            % [uM/s] ligand-receptor complex
% dE = -react2;                   % [uM/s] free enzyme
% dLRE = react2;                  % [uM/s] ligand-receptor-enzyme complex
% dS = react4-react3;             % [uM/s] substrate
% dP = react3-react4-react5;      % [uM/s] product
% dPnuc = react5/Vratio;          % [uM/s] nuclear product
% dydt = [dR;dLR;dE;dLRE;dS;dP;dPnuc]; % Reassemble differential equations
% algvars = [react1,react2,react3,react4,react5]; % optional for seeing fluxes
dESR = -J1;
dFADD = -J1; 
dprocas8 = -J2;
dESR_FADD = J1-J2;
dprocas8_ESR_FADD = J2; 
dE = -J3;
dcas8 = J3-J4; 
dcas9 = -J5+J6;
dcytc = -J6+J7; 
dbax_bak = J8-J7; 
dbax = -J8; 
dbak = J9-J8; 
dbcl = -J9+J10; 
dISR = -J10; 
dNICD = -J11 + J12; 
dJaggedNotch = -J12;
dcas367 = -J4-J5;

dydt = [dESR; dFADD; dprocas8; dESR_FADD; dprocas8_ESR_FADD; dE; dcas8; dcas9; dcytc; dbax_bak; dbax; dbak; dbcl; dISR; dNICD; dJaggedNotch; dcas367];
algvars = [J1, J2, J3, J4, J5, J6, J7, J8, J9, J10, J11, J12];











