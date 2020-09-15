% EnergeticsModelScript
 function out = EnergeticsModelScript(TAN, CRtot, TEP, Ox_capacity, x_ATPase)

% INPUTS
%   TAN: total adenine nucleotide pool (M per liter cell)
%   CRtot: total creatine pool (M per liter cell)
%   TEP: total exchangeable phosphate pool (M per liter cell)
%   Ox_capacity: oxidative capacity (relative to control)
%   x_ATPase: ATP hydrolysis rate: M / s / (liter cytosol)
%
% OUTPUTS
%   MgATP_cytoplasm: Mg-ATP concentration (mM)
%   MgADP_cytoplasm: Mg-ADP concentration (mM)
%   fPi_cytoplasm: free unchelated Pi concentration (mM)
%   J_ATPase_tissue: ATP hydrolysis rate in units of (mol/s/(l cell))
%   MVO2_tissue: oxygen consumption rate in units of (umol o2/min/g tissue)
%   dGrATPase: ATP hydrolysis potential in kJ / mol
%   PCrATP: CrP/ATP ratio (unitless)
%   ATP_cyto: cyto [ATP]
%   ADP_cyto: cyto [ADP]
%   Pi_cyto: cyto [Pi]

% TAN = input(1); % (M per liter cell)
% CRtot = input(2); % (M per liter cell)
% TEP = input(3); % (M per liter cell)
% Ox_capacity = input(4); % oxidative capacity relative to control sham
% x_ATPase = input(5); % ATP hydrolysis rate: M / s / (liter cytosol)

% Arbitrary initial Cr and PCr values

Cr = CRtot*0.7; % M per liter cell
PCr = CRtot*0.3; % M per liter cell

%% Set the desired temperature
T = 37.5;

%% pH buffering parameters
BX(1) = 0.02;  % (M) matrix
BX(2) = 1;     % (M) cytoplasm
BX(3) = 1;     % (M) intermembrane space
%% Proton cytoplasm binding constant(s)
K_BX(1) = 1e-7; % (M)
K_BX(2) = 10^(-7.2); % (M)
K_BX(3) = 10^(-7.2); % (M)

%% Tissue and cell composition constants
% All constants based on [https://doi.org/10.1152/ajpheart.00478.2003]
% with the mitochonrial content scaled by measured oxidative capacity

Rc_t = 0.73; % cell volume/tissue volume 
rho = 1.05; % tissue mass density 
Vcyto = 0.6801 + 0.2882 *(1 - Ox_capacity); % cytosolic volume fraction (mL cyto / mL cell) 
Wcyto = 0.8425; % cytosol water fraction (mL water / mL cytosol)
Vmito = 0.2882 * Ox_capacity; % mito volume fraction (mL mito / mL cell)
Wmito_x = 0.6514; % mL matrix water / mito volume
Wmito_i = 0.0724; % mL IMS water / mito volume
var = [T Ox_capacity];

%% defining mitochondrial total pool concentrations

Ntot = 2.97e-3; % NADH
Qtot = 20e-3; % coenzyme Q
Atot = 10e-3; % ADP
Ctot = 2e-4; % cytochrome C

%% Assigning initial conditions to satisfy measured metabolite pools

% matrix: fixed
ATP_x = 1e-3;%2e-3;
ADP_x = Atot - ATP_x;
Pi_x = 1e-3;%1.82e-3;

% cytoplasm: determined by metabolite pools
ATP_c = (TAN-Vmito*Wmito_x*Atot)/(Vcyto*Wcyto+Vmito*Wmito_i);% (ATP_cell-Vmito*Wmito_x*ATP_x)/(Vcyto*Wcyto+Vmito*Wmito_i)
AMP_c = 1e-9;
ADP_c = 1e-9;
PCr_c = PCr/(Vcyto*Wcyto);
Cr_c = Cr/(Vcyto*Wcyto);
CRt = (PCr_c+Cr_c);

% im:
ATP_i = ATP_c;
ADP_i = ADP_c;
AMP_i = AMP_c;

Pi_c = (TEP-Vcyto*Wcyto*(ATP_c*2+ADP_c+PCr_c)-Vmito*Wmito_i*(ATP_i*2+ADP_i)-Vmito*Wmito_x*(ATP_x+Pi_x))/(Vcyto*Wcyto+Vmito*Wmito_i);
Pi_i = Pi_c;


%% Assign initial conditions for in vitro simulations

x0 = [];
% MATRIX
x0(1) = Ntot - 1e-12; % NAD, matrix
x0(2) = 1e-12; % NADH, matrix
x0(3) = ADP_x; % ADP, matrix
x0(4) = ATP_x; % ATP, matrix
x0(5) = Pi_x; % Pi, matrix
x0(6) = Qtot - 1e-12; % coQ
x0(7) = 1e-12; % coQH2
% CYTOPLASM
x0(8) = ATP_c; % ATP, cytoplasm 
x0(9) = ADP_c; % ADP, cytoplasm
x0(10) = Pi_c; % Pi, cytoplasm 
x0(11) = PCr_c; % Pi, cytoplasm 
x0(12) = Cr_c; % Pi, cytoplasm 
x0(13) = AMP_c; % ADP, cytoplasm
% IM SPACE
x0(14) = Ctot - 1e-12; % cytoC (ox), im 
x0(15) = 1e-12; % cytoC (red), im 
x0(16) = Pi_i; % Pi, cytoplasm 
x0(17) = ADP_i; % ADP, cytoplasm
x0(18) = ATP_i; % ATP, cytoplasm 
x0(19) = AMP_i; % ADP, cytoplasm
% CATIONS

x0(20) = 10^-(7.4); % H, matrix
x0(21) = 1.2e-3; % Mg, matrix

x0(22) = 0.060; % K, matrix
x0(23) = 10^(-7.15); % H, cytoplasm 
x0(24) = 1e-3; % Mg, cytoplasm 
x0(25) = 130e-3; % K, cytoplasm 
x0(26) = 10^(-7.15); % H, im
x0(27) = 1e-3; % Mg, im 
x0(28) = 130e-3; % K, im 
% MEMBRANE POTENTIAL
x0(29) =  0 ; % membrane potential, IM to matrix


%% Assign parameter values (from Bazil et al. 2016)
params = [];
params(1) = 1.1273e-1 ;    % dehydrogenase
params(2) = 0.25;    % SDH
params(3) = x_ATPase; % x_ATPase THE ATP HYDROLYSIS RATE
% params(4) = ; % parameter not used
params(5) =  812 ;    % F1F0
params(6) = 1.66953e-4 ;    % ETC1Activity
params(7) =  0.00158489319246111 ;    % ETC3Activity
params(8) = 8.91e-2;    % ETC4Activity
params(9) =  0.3e3 ;    % proton leak
params(10) =  3.34e7 ;    % PiH
params(11) = 0.14125 ;    % ANT
params(12) = 4760000; % k1_KH
params(13) =  1000 ;    % AMPPERM
params(14) = 5000 ;    % ADPPERM
params(15) = 1000 ;    % ATPPERM
params(16) = 1000 ;    % PIPERM
params(17) =  1e9 ;    % HPERM
params(18) =  100 ;    % KPERM
params(19) =  100 ;    % MPERM
params(20) = 6.7568e-06 ;    % C, im to matrix


%% Make sure all the variables except the membrane potential cannot become negative
varlist = 1:28;
options = odeset('NonNegative', varlist);

[t,x] = ode15s(@dXdT_energetics, [0 10000], x0, options, var, BX, K_BX, params,1 ); % to steady state
[y,J] = dXdT_energetics(t(end),x(end,:)',var,BX,K_BX,params,0);

%% Results...

ATP_cyto = x(end,8);
ADP_cyto = x(end,9);
AMP_cyto = x(end,13);
Pi_cyto= x(end,10);
ATP_im = x(end,18);
ADP_im = x(end,17);
AMP_im = x(end,19);
Pi_im= x(end,16);
ATP_matrix = x(end,4);
ADP_matrix = x(end,3);
Pi_matrix= x(end,5);
PCr_cyto = x(end,11);
Cr_cyto = x(end,12);

% % double check pool concentrations 
% % (useful in debugging comment out for use in multiscle model)
% totCr = (PCr_cyto + Cr_cyto)*Vcyto*Wcyto;
% TAN = (ATP_cyto + ADP_cyto + AMP_cyto)*Vcyto*Wcyto +...
%     (ATP_im + ADP_im + AMP_im)*Vmito*Wmito_i + (ATP_matrix + ADP_matrix )*Vmito*Wmito_x;
% TEP = (2*ATP_cyto + ADP_cyto + Pi_cyto + PCr_cyto)*Vcyto*Wcyto +...
%     (2*ATP_im + ADP_im + Pi_im)*Vmito*Wmito_i + (ATP_matrix + Pi_matrix )*Vmito*Wmito_x;

% PCr/ATP
PCr_cell = Vcyto*Wcyto*PCr_cyto;
Cr_cell = Vcyto*Wcyto*Cr_cyto;
ATP_cell = Vcyto*Wcyto*ATP_cyto+Vmito*Wmito_x*ATP_matrix+Vmito*Wmito_i*ATP_im;
PCrATP = PCr_cell/ATP_cell;
% Fluxes
J(3)=J(3);
lcell = 1e-3*Rc_t/rho; % l of cell per g tissue; Vinnakota etal Am J Physiol Heart Circ Physiol. 2004 May;286(5):H1742-9
J_ATPase_cell = J(3)*Vcyto; % J_ATPASE_cytoplasm % J_ATPase: mol/s/(l cell)

J_ATPase_tissue = J_ATPase_cell*lcell; % mol/s/(g tissue)

J_C4 = J(9); % J_ETC4_im_to_matrix
J_ANT = J(12);
J_HLEAK = J(10);
MVO2_tissue = (J_C4/2)*Vmito*Rc_t/1000/rho*60*1e6; % umol o2/min/g tissue

% 
MgATP_cytoplasm = y(1); % mol/(l cyto water)
MgADP_cytoplasm = y(2);
fPi_cytoplasm = y(3);
dGrATPase = y(4);

out = [MgATP_cytoplasm*1000 MgADP_cytoplasm*1000 fPi_cytoplasm*1000 J_ATPase_tissue MVO2_tissue dGrATPase PCrATP ATP_cyto ADP_cyto Pi_cyto];

% end