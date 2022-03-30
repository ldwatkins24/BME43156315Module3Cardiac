% function myocyte_ApoptosisCompute
% solve the ODEs for myocyte apoptosis 
% Vignesh Valaboju
clear all; close all;


%% define parameters
% k1f = 1;    % [uM^-1 s^-1] react1 forward rate constant
% k1r = 1;    % [s^-1] react1 reverse rate constant
% k2f = 1;    % [uM^-1 s^-1] react2 forward rate constant
% k2r = 1;    % [s^-1] react2 reverse rate constant
% kcat = 1;   % [s^-1] catalytic rate constant for enzyme
% Km = 1;     % [uM] Michaelis constant for enzyme
% k4 = 1;     % [s^-1] react4 rate constant
% k5 = 1;
% Vratio = 0.1;
% 
% L = 10;   % [uM] concentration of ligand
% Rtot = 1;   % [uM] total concentration of receptor
% Etot = 10;   % [uM] total concentration of enzyme
% Stot = 1;   % [uM] total concentration of substrate
% params = {k1f,k1r,k2f,k2r,kcat,Km,k4,k5,Vratio,L};

% rate constants
k1f = 1; k1r = 1;
k2f = 1; k2r = 1; 
k3 = 1; km = 1; 
k4 = 1; 
k5 = 1; 
k6 = 1; 
k7 = 1;
k8f = 1; k8r = 1; 
k9 = 1; 
k10 = 1; 
k11 = 1; 
k12 = 1;
params = {k1f, k1r, k2f, k2r, k3, km, k4, k5, k6, k7, k8f, k8r, k9, k10, k11, k12};

%initial concentrations 
ESR = 10; %make ESR very big so it never runs out  
FADD = 10; %make FADD very big so it never runs out  
ESR_FADD = 1; 
procas8 = 10; %make procaspase 8 very big so it never runs out
procas8_ESR_FADD = 1; 
E = 10; %make the enzyme very big so it never runs out
cas8 = 1;
cas9 = 1; 
cytc = 1;
bax_bak = 1; 
bax = 10; %make bax very big so it never runs out
bak = 1; 
bcl = 1;
ISR = 1; 
NICD = 1; 
JaggedNotch = 1; 
cas367 = 0; %starting off with none of our measured output 

%% Run single simulation
% 'R','LR','E','LRE','S','P','Pnuc'
% y0 = [Rtot; 0; Etot; 0; Stot; 0; 0];   % or you could load saved data like: y0 = load('yfinal.dat');
y0 = [ESR; FADD; procas8; ESR_FADD; procas8_ESR_FADD; E; cas8; cas9; cytc; bax_bak; bax; bak; bcl; ISR; NICD; JaggedNotch; cas367];
tspan = [0 10];
options = [];%odeset('MaxStep',5e-3);
[t,y] = ode23(@myocyte_apoptosisODEs,tspan,y0,options,params);

yfinal = y(end,:)';
% save -ascii 'yfinal.dat' yfinal;    % save final values to a file

% plot timecourse
subplot(1,3,1);
plot(t,y);
xlabel('Time (sec)'); ylabel('y(t) (\muM)'); %legend('R','LR','E','LRE','S','P','Pnuc');

% optional: re-evaluate ODE function after solving ODEs to calculate algebraic variables 
for tstep=1:length(t),
    [~,algvars(tstep,:)]=myocyte_apoptosisODEs(t(tstep),y(tstep,:),params);
end
subplot(1,3,2);
plot(t,algvars(:,1),t,algvars(:,2),t,algvars(:,3),t,algvars(:,4),t,algvars(:,5));
xlabel('Time (sec)'); ylabel('fluxes(t) (\muM/s)');
% legend('react1','react2','react3','react4','react5');

%% Run dose response over a range of total ligand concentrations
% paramRange = 10.^[-2:.1:2];
% for i=1:length(paramRange)
%     L = paramRange(i);
%     params = {k1f,k1r,k2f,k2r,kcat,Km,k4,k5,Vratio,L};
%     y0 = [Rtot; 0; Etot; 0; Stot; 0; 0];
%     tspan = [0 10];
%     options = [];
%     [t,y] = ode15s(@toyModel2ODEfunc,tspan,y0,options,params);
%     Pnucfinal(i) = y(end,end);
% end
% subplot(1,3,3);
% semilogx(paramRange,Pnucfinal);
% xlabel('Ligand (\muM)'); ylabel('Steady state nuclear Product (\muM)');