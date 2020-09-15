% function DXDT = dXdT_cardiovascular_mechanics(t,x,...)
%
%
% Output parameters:
%   DXDT  time derivatives of the model
%
% Mandatory input parameters:
%   t     time
%   x     state variables at time t%   
%
% State Variables:
% []
%
% Parameters:
% []


function dXdT = dXdT_cardiovascular_mechanics(t,x,adjvar,input)

%% Parameters
CO_target = input(1);
%% Parameters
stim_period = input(2); % Added for Calcium handling (9/5 BM)
% Heart Model Parameters
Vw_LV     = input(3); % LV wall volume, mL 
Vw_SEP    = input(4); % Septal wall volume, mL 
Vw_RV     = input(5); % RV wall volume, mL 
R_TAC = input(6); % resistance of TAC , mmHg*sec/mL

Lsref     = 1.9; % Resting SL, micron
eta       = 1; % visconsity, mmHg s /micron
L_rest_pas = 0.9;

% Sarcomere geometry parameters
L_thick = 1.67; % Length of thick filament, um
L_hbare = 0.10; % Length of bare region of thick filament, um
L_thin  = 1.20; % Length of thin filament, um
deltaR  = 0.010; % um
% 

% Metabolite levels (use the values from Tewari's paper) (9/5 BM)
MgATP_LV  = input(7); % cytosolic Mg-ATP concentration, mM
MgADP_LV  = input(8); % cytosolic Mg-ADP concentration, mM
Pi_LV     = input(9); % cytosolic Pi concentration, mM
MgATP_SEP = input(10); % cytosolic Mg-ATP concentration, mM
MgADP_SEP = input(11); % cytosolic Mg-ATP concentration, mM
Pi_SEP    = input(12); % cytosolic Mg-ATP concentration, mM
MgATP_RV  = input(13); % cytosolic Mg-ATP concentration, mM
MgADP_RV  = input(14); % cytosolic Mg-ATP concentration, mM
Pi_RV     = input(15); % cytosolic Mg-ATP concentration, mM
a    = input(16);
b  = input(17);
c  = input(18);
Ca0     = input(19);
Amref_LV = input(20);
Amref_SEP = input(21);
Amref_RV = input(22);
flag_swap_metabolite = input(23);
K_Pi = 4.00; K_T = 0.4897; K_D = 0.194;% Used the values from Tewari etal JMCC (9/5 BM)

%% para 4

ka      = 559.5568; % myosin-actin attach rate constant, 1/sec
kd      =  304.6708;% myosin-actin detach rate constant, 1/sec
k1      = 112.3727; % transition A1 to A2 rate constant, 1/sec
km1     = 21.296; % transition A2 to A1 rate constant, 1/sec
k2      = 811.72; % transition A2 to A3 rate constant, 1/sec
km2     = 43.25;% transition A3 to A2 rate constant, 1/sec
k3      = 144.5586; % transition A3 to P rate constant, 1/sec
alpha1  = 10.0; % Stretch sensing parameter for k1 and k?1, 1/um
alpha2  = 9.1; % Stretch sensing parameter for k2 and k?2, 1/um
alpha3  =  0.1*59.3; % Stretch sensing parameter for k3, 1/um
s3      = 9.9e-3;  % Strain at which k3 is minimum, um

K_coop =      9.6846; % Campbell et al, Biophysical Journal 2018,
k_on =   101.1850; % Campbell et al, Biophysical Journal 2018
k_off = 723.8520; % manually tuned parameter!

% transitions between super relaxed state and non relaxed state
ksr = adjvar(4);
kforce = adjvar(5)/7.5;  %dived by kPa to mmHg conversion rate
kmsr =    50.032; % made-up number

kstiff1 = 1.4 * 1.7351e+03 * 7.5; % kPa/um (9/5 BM)
kstiff2 = 1.4 * 45545 * 7.5; % kPa/um (9/5 BM)

k_passive = 25; % mN / mm^2 / micron % for mean SHAM rat and TAC rat 1

% correcting rate constants for metabolite levels in LV, SEP, and RV
kd_LV  = kd*(Pi_LV/K_Pi)/(1.0 + Pi_LV/K_Pi);
k1_LV  = k1/(1.0 + Pi_LV/K_Pi);
km2_LV = km2*(MgADP_LV/K_D)/(1.0 + MgADP_LV/K_D + MgATP_LV/K_T);
k3_LV  = k3*(MgATP_LV/K_T)/(1.0 + MgATP_LV/K_T + MgADP_LV/K_D);

kd_SEP  = kd*(Pi_SEP/K_Pi)/(1.0 + Pi_SEP/K_Pi);
k1_SEP  = k1/(1.0 + Pi_SEP/K_Pi);
km2_SEP = km2*(MgADP_SEP/K_D)/(1.0 + MgADP_SEP/K_D + MgATP_SEP/K_T);
k3_SEP  = k3*(MgATP_SEP/K_T)/(1.0 + MgATP_SEP/K_T + MgADP_SEP/K_D);

kd_RV  = kd*(Pi_RV/K_Pi)/(1.0 + Pi_RV/K_Pi);
k1_RV  = k1/(1.0 + Pi_RV/K_Pi);
km2_RV = km2*(MgADP_RV/K_D)/(1.0 + MgADP_RV/K_D + MgATP_RV/K_T);
k3_RV  = k3*(MgATP_RV/K_T)/(1.0 + MgATP_RV/K_T + MgADP_RV/K_D);

% Lumped circulatory parameters
C_Ao = 0.0022045;  % Proximal aortic compliance, mL/mmHg
C_SA = 0.0077157; % Systemic arterial compliance, mL/mmHg

C_SV = 2.5; % Systemic venous compliance, mL/mmHg  DAB 10/7/2018
C_PV = 0.25; % Pulmonary venous compliance, mL/mmHg
C_PA = 0.013778; % Pulmonary arterial compliance, mL/mmHg
R_Ao   = 2.5; % resistance of aorta , mmHg*sec/mL
R_SA   = adjvar(7);% mmHg*sec/mL; % Systemic vasculature resistance, mmHg*sec/mL
if flag_swap_metabolite == 1
    R_PA   = 12/CO_target*60; % Pulmonary vasculature resistance, mmHg*sec/mL % Match the old code(9/5 BM) DAB change 9/15
else 
    R_PA   = adjvar(7)*12/88; % Pulmonary vasculature resistance, mmHg*sec/mL % Match the old code(9/5 BM) DAB change 9/15
end
R_SV   = 0.25; 
R_PV   = 0.25; 
R_vlv  = 0.05; %  valve resistance, mmHg*sec/mL
R_AV   = R_vlv + R_TAC; % resistance across aortic valve
R_tAo  = 0.5;
R_tSA  = 4;

Kse    = 50000; % series element elastance, mmHg/micron (Changed to match the 

%% Variables

xm_LV  = x(1); % LV heart geometry variable, cm
xm_SEP = x(2); % septum heart geometry variable, cm
xm_RV  = x(3); % RV heart geometry variable, cm
ym     = x(4); % Heart geometry variable, cm
SL_LV  = x(5); % sarcomere length, LV, micron
SL_SEP = x(6); % sarcomere length, septum, micron
SL_RV  = x(7); % sarcomere length, RV, micron
V_LV   = x(8); % volume LV, mL
V_RV   = x(9); % volume RV, mL

P1_0_LV = x(10); % 0th moment state A1, LV
P1_1_LV = x(11); % 1st moment state A1, LV
P1_2_LV = x(12); % 2nd moment state A1, LV
P2_0_LV = x(13); % 0th moment state A2, LV
P2_1_LV = x(14); % 1st moment state A2, LV
P2_2_LV = x(15); % 2nd moment state A2, LV
P3_0_LV = x(16); % 0th moment state A3, LV
P3_1_LV = x(17); % 1st moment state A3, LV
P3_2_LV = x(18); % 2nd moment state A3, LV
N_LV    = x(19); % nonpermissive fraction LV
U_NR_LV = x(20); % U_NR represents the Non relaxed state

P1_0_SEP = x(21); % 0th moment state A1, SEP
P1_1_SEP = x(22); % 1st moment state A1, SEP
P1_2_SEP = x(23); % 2nd moment state A1, SEP
P2_0_SEP = x(24); % 0th moment state A2, SEP
P2_1_SEP = x(25); % 1st moment state A2, SEP
P2_2_SEP = x(26); % 2nd moment state A2, SEP
P3_0_SEP = x(27); % 0th moment state A3, SEP
P3_1_SEP = x(28); % 1st moment state A3, SEP
P3_2_SEP = x(29); % 2nd moment state A3, SEP
N_SEP    = x(30); % nonpermissive fraction SEP
U_NR_SEP = x(31); % U_NR represents the Non relaxed state

P1_0_RV = x(32); % 0th moment state A1, RV
P1_1_RV = x(33); % 1st moment state A1, RV
P1_2_RV = x(34); % 2nd moment state A1, RV
P2_0_RV = x(35); % 0th moment state A2, RV
P2_1_RV = x(36); % 1st moment state A2, RV
P2_2_RV = x(37); % 2nd moment state A2, RV
P3_0_RV = x(38); % 0th moment state A3, RV
P3_1_RV = x(39); % 1st moment state A3, RV
P3_2_RV = x(40); % 2nd moment state A3, RV
N_RV    = x(41); % nonpermissive fraction RV
U_NR_RV = x(42); % U_NR represents the Non relaxed state

V_SV   = x(43); % volume of systemic veins
V_PV   = x(44); % volume of pulmonary veins
V_SA   = x(45); % volume of systemic arterys
V_PA   = x(46); % volume of pulmonary arterys
V_Ao   = x(47); % volume of proximal aorta
% % 

%% Ca Calculations
phi = mod(t+0.0001,stim_period)/stim_period;
Ca_i = (a/phi)*exp(-b*(log(phi)+c)^2) + Ca0;
%---------------------------------------------------------------------

%% Heart and Sarcomere Model

Vm_LV  = (pi/6)*xm_LV*(xm_LV^2 + 3*ym^2);
Vm_SEP = (pi/6)*xm_SEP*(xm_SEP^2 + 3*ym^2);
Vm_RV  = (pi/6)*xm_RV*(xm_RV^2 + 3*ym^2);
Am_LV  = pi*(xm_LV^2 + ym^2);
Am_SEP = pi*(xm_SEP^2 + ym^2);
Am_RV  = pi*(xm_RV^2 + ym^2);
Cm_LV  = 2*xm_LV/(xm_LV^2 + ym^2);
Cm_SEP = 2*xm_SEP/(xm_SEP^2 + ym^2);
Cm_RV  = 2*xm_RV/(xm_RV^2 + ym^2);
z_LV   = 3*Cm_LV*Vw_LV/(2*Am_LV);
z_SEP  = 3*Cm_SEP*Vw_SEP/(2*Am_SEP);
z_RV   = 3*Cm_RV*Vw_RV/(2*Am_RV);

epsf_LV = 0.5*log(Am_LV/Amref_LV) - 0.083333*z_LV^2 - 0.019*z_LV^4;
epsf_SEP = 0.5*log(Am_SEP/Amref_SEP) - 0.083333*z_SEP^2 - 0.019*z_SEP^4;
epsf_RV = 0.5*log(Am_RV/Amref_RV) - 0.083333*z_RV^2 - 0.019*z_RV^4;
SLo_LV = Lsref*exp(epsf_LV); 
SLo_SEP = Lsref*exp(epsf_SEP); 
SLo_RV = Lsref*exp(epsf_RV);

% Collagen force
SLcollagen = 2.25; % threshold for collagen activation,(um)
PConcollagen = 0.01; % contriubtion of collagen (unitless)
PExpcollagen = 70; % expresion of collagen (1/um)
sigma_collagen_LV  = PConcollagen*(exp(PExpcollagen*(SLo_LV - SLcollagen)) - 1).*(SLo_LV > SLcollagen);
sigma_collagen_SEP = PConcollagen*(exp(PExpcollagen*(SLo_SEP - SLcollagen)) - 1).*(SLo_SEP > SLcollagen);
sigma_collagen_RV  = PConcollagen*(exp(PExpcollagen*(SLo_RV - SLcollagen)) - 1).*(SLo_RV > SLcollagen);

sigmapas_LV  = k_passive*(SLo_LV/2-L_rest_pas)  + sigma_collagen_LV ;
sigmapas_SEP = k_passive*(SLo_SEP/2-L_rest_pas) + sigma_collagen_SEP;
sigmapas_RV  = k_passive*(SLo_RV/2-L_rest_pas)   + sigma_collagen_RV;

% Sarcomere geometry
sovr_ze = min(L_thick*0.5, SL_LV*0.5);
sovr_cle = max(SL_LV*0.5 - (SL_LV-L_thin),L_hbare*0.5);
L_sovr = sovr_ze - sovr_cle; % Length of single overlap region
N_overlap_LV = L_sovr*2/(L_thick - L_hbare);
% sov_thin_LV = L_sovr/L_thin;

sovr_ze = min(L_thick*0.5, SL_SEP*0.5);
sovr_cle = max(SL_SEP*0.5 - (SL_SEP-L_thin),L_hbare*0.5);
L_sovr = sovr_ze - sovr_cle; % Length of single overlap region
N_overlap_SEP = L_sovr*2/(L_thick - L_hbare);
% sov_thin_SEP = L_sovr/L_thin;

sovr_ze = min(L_thick*0.5, SL_RV*0.5);
sovr_cle = max(SL_RV*0.5 - (SL_RV-L_thin),L_hbare*0.5);
L_sovr = sovr_ze - sovr_cle; % Length of single overlap region
N_overlap_RV = L_sovr*2/(L_thick - L_hbare);
% sov_thin_RV = L_sovr/L_thin;

% Active forces
sigmaact_LV  = N_overlap_LV*(kstiff2*deltaR*(P3_0_LV) + kstiff1*(P2_1_LV + P3_1_LV)); % mmHg * normalised force
sigmaact_SEP = N_overlap_SEP*(kstiff2*deltaR*(P3_0_SEP) + kstiff1*(P2_1_SEP+P3_1_SEP)); % mmHg * normalised force
sigmaact_RV  = N_overlap_RV*(kstiff2*deltaR*(P3_0_RV) + kstiff1*(P2_1_RV+P3_1_RV)); % mmHg * normalised force

% Total forces
sigmaf_LV = -Kse*(SL_LV - SLo_LV);
sigmaf_SEP = -Kse*(SL_SEP - SLo_SEP);
sigmaf_RV = -Kse*(SL_RV - SLo_RV);

Tm_LV  = (Vw_LV*sigmaf_LV/(2*Am_LV))*(1 + (z_LV^2)/3 + (z_LV^4)/5);
Tm_SEP = (Vw_SEP*sigmaf_SEP/(2*Am_SEP))*(1 + (z_SEP^2)/3 + (z_SEP^4)/5);
Tm_RV  = (Vw_RV*sigmaf_RV/(2*Am_RV))*(1 + (z_RV^2)/3 + (z_RV^4)/5);

Tx_LV  = Tm_LV*2*xm_LV*ym/(xm_LV^2 + ym^2);
Tx_SEP = Tm_SEP*2*xm_SEP*ym/(xm_SEP^2 + ym^2);
Tx_RV  = Tm_RV*2*xm_RV*ym/(xm_RV^2 + ym^2);

Ty_LV  = Tm_LV*(-xm_LV^2 + ym^2)/(xm_LV^2 + ym^2);
Ty_SEP = Tm_SEP*(-xm_SEP^2 + ym^2)/(xm_SEP^2 + ym^2);
Ty_RV  = Tm_RV*(-xm_RV^2 + ym^2)/(xm_RV^2 + ym^2);

% ventricular pressure
ptrans1 = 2*Tx_LV/ym;
ptrans3 = 2*Tx_RV/ym;
P_LV = -ptrans1;
P_RV = ptrans3;

%% Lumped circulatory model
P_SV = V_SV/C_SV;
P_PV = V_PV/C_PV;
P_PA = V_PA/C_PA;

% Ao valves closed equations
QOUT_LV = 0;
P_Ao = (C_SA*R_Ao*R_SA*V_Ao + C_SA*R_Ao*R_tSA*V_Ao + C_SA*R_SA*R_tSA*V_Ao + C_Ao*R_SA*R_tAo*V_SA + C_Ao*C_SA*P_SV*R_tSA*R_tAo)/(C_Ao*C_SA*(R_Ao*R_SA + R_Ao*R_tSA + R_SA*R_tSA + R_SA*R_tAo + R_tSA*R_tAo));
P_SA = (C_Ao*R_Ao*R_SA*V_SA + C_SA*R_SA*R_tSA*V_Ao + C_Ao*R_SA*R_tAo*V_SA + C_Ao*C_SA*P_SV*R_Ao*R_tSA + C_Ao*C_SA*P_SV*R_tSA*R_tAo)/(C_Ao*C_SA*(R_Ao*R_SA + R_Ao*R_tSA + R_SA*R_tSA + R_SA*R_tAo + R_tSA*R_tAo));
Q_Ao = -(C_Ao*R_SA*V_SA - C_SA*R_SA*V_Ao - C_SA*R_tSA*V_Ao + C_Ao*C_SA*P_SV*R_tSA)/(C_Ao*C_SA*(R_Ao*R_SA + R_Ao*R_tSA + R_SA*R_tSA + R_SA*R_tAo + R_tSA*R_tAo));
if (P_Ao < P_LV)*(V_LV>0) 
  % Ao valve open equations 
  P_SA    = (C_Ao*R_Ao*R_SA*R_tAo*V_SA + C_Ao*R_Ao*R_SA*R_AV*V_SA + C_SA*R_SA*R_tSA*R_AV*V_Ao + C_Ao*R_SA*R_tAo*R_AV*V_SA + C_Ao*C_SA*P_SV*R_Ao*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_Ao*R_tSA*R_AV + C_Ao*C_SA*P_LV*R_SA*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_tSA*R_tAo*R_AV)/(C_Ao*C_SA*(R_Ao*R_SA*R_tAo + R_Ao*R_SA*R_AV + R_Ao*R_tSA*R_tAo + R_Ao*R_tSA*R_AV + R_SA*R_tSA*R_tAo + R_SA*R_tSA*R_AV + R_SA*R_tAo*R_AV + R_tSA*R_tAo*R_AV));
  QOUT_LV = -(C_SA*R_Ao*R_SA*V_Ao + C_SA*R_Ao*R_tSA*V_Ao + C_SA*R_SA*R_tSA*V_Ao + C_Ao*R_SA*R_tAo*V_SA - C_Ao*C_SA*P_LV*R_Ao*R_SA - C_Ao*C_SA*P_LV*R_Ao*R_tSA - C_Ao*C_SA*P_LV*R_SA*R_tSA - C_Ao*C_SA*P_LV*R_SA*R_tAo - C_Ao*C_SA*P_LV*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_tSA*R_tAo)/(C_Ao*C_SA*(R_Ao*R_SA*R_tAo + R_Ao*R_SA*R_AV + R_Ao*R_tSA*R_tAo + R_Ao*R_tSA*R_AV + R_SA*R_tSA*R_tAo + R_SA*R_tSA*R_AV + R_SA*R_tAo*R_AV + R_tSA*R_tAo*R_AV));
  Q_Ao    = -(C_Ao*R_SA*R_tAo*V_SA + C_Ao*R_SA*R_AV*V_SA - C_SA*R_SA*R_AV*V_Ao - C_SA*R_tSA*R_AV*V_Ao - C_Ao*C_SA*P_LV*R_SA*R_tAo - C_Ao*C_SA*P_LV*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_tSA*R_AV)/(C_Ao*C_SA*(R_Ao*R_SA*R_tAo + R_Ao*R_SA*R_AV + R_Ao*R_tSA*R_tAo + R_Ao*R_tSA*R_AV + R_SA*R_tSA*R_tAo + R_SA*R_tSA*R_AV + R_SA*R_tAo*R_AV + R_tSA*R_tAo*R_AV));
end
QIN_LV  = max((P_PV - P_LV)/R_PV,0); 
QIN_RV  = max((P_SV - P_RV)/R_SV,0); 
QOUT_RV = max((P_RV - P_PA)/R_vlv,0)*(V_RV>0); 

%--------------------------------------------
% TriSeg
dXdT(1) = (-V_LV - 0.5*Vw_LV - 0.5*Vw_SEP + Vm_SEP - Vm_LV)/V_LV; %
dXdT(2) = (Tx_LV + Tx_SEP + Tx_RV); % 
dXdT(3) = (V_RV + 0.5*Vw_RV + 0.5*Vw_SEP + Vm_SEP - Vm_RV)/V_LV; % 
dXdT(4) = (Ty_LV + Ty_SEP + Ty_RV);  %


%% Myofiber Mechanics: SL_LV

% Calculations for stretch-senstive rates    
f_alpha1o = (P1_0_LV - alpha1*P1_1_LV + 0.5*(alpha1*alpha1)*P1_2_LV);
f_alpha1i = (P1_1_LV - alpha1*P1_2_LV);

f_alpha0o = (P2_0_LV + alpha1*P2_1_LV + 0.5*alpha1*alpha1*P2_2_LV);
f_alpha0i = (P2_1_LV + alpha1*P2_2_LV);

f_alpha2o = (P2_0_LV - alpha2*P2_1_LV + 0.5*(alpha2*alpha2)*P2_2_LV);
f_alpha2i = (P2_1_LV - alpha2*P2_2_LV);

f_alpha3o = (P3_0_LV + alpha3*(s3*s3*P3_0_LV + 2.0*s3*P3_1_LV + P3_2_LV));
f_alpha3i = (P3_1_LV + alpha3*(s3*s3*P3_1_LV + 2.0*s3*P3_2_LV));

dSL_LV = (sigmaf_LV - sigmapas_LV - sigmaact_LV)/eta;% 

P0_LV = 1.0 - N_LV - P1_0_LV - P2_0_LV - P3_0_LV; % DAB 10/8/2019
dXdT(10) = ka*P0_LV*U_NR_LV*N_overlap_LV - kd_LV*P1_0_LV - k1_LV*f_alpha1o + km1*f_alpha0o; % DAB 10/8/2019
dXdT(11) = dSL_LV*P1_0_LV - kd_LV*P1_1_LV - k1_LV*f_alpha1i + km1*f_alpha0i;
dXdT(12) = 2*dSL_LV*P1_1_LV - kd_LV*P1_2_LV - k1_LV*P1_2_LV + km1*P2_2_LV;

dXdT(13) = -km1*f_alpha0o - k2*f_alpha2o + km2_LV*P3_0_LV + k1_LV*f_alpha1o;
dXdT(14) = dSL_LV*P2_0_LV - k1_LV*f_alpha0i - k2*f_alpha2i + km2_LV*P3_1_LV + k1_LV*f_alpha1i;
dXdT(15) = 2*dSL_LV*P2_1_LV - km1*P2_2_LV - k2*P2_2_LV + km2_LV*P3_2_LV + k1_LV*P1_2_LV;

dXdT(16) = +k2*f_alpha2o - km2_LV*P3_0_LV - k3_LV*f_alpha3o;
dXdT(17) = dSL_LV*P3_0_LV + k2*f_alpha2i - km2_LV*P3_1_LV - k3_LV*f_alpha3i;
dXdT(18) = 2*dSL_LV*P3_1_LV + k2*P2_2_LV - km2_LV*P3_2_LV - k3_LV*P3_2_LV;


U_SR_LV = 1 - U_NR_LV;
Jon = k_on*Ca_i*N_LV*(1 + K_coop*(1 - N_LV));
Joff = k_off*P0_LV*(1 + K_coop*N_LV);
dXdT(19) = - Jon + Joff; % dN_LV / dt
dXdT(20) = ksr * (1 + kforce * sigmaact_LV) * U_SR_LV - kmsr*U_NR_LV ; % DAB 10/8

%% Myofiber Mechanics: SL_SEP

% Calculations for stretch-senstive rates    
f_alpha1o = (P1_0_SEP - alpha1*P1_1_SEP + 0.5*(alpha1*alpha1)*P1_2_SEP);
f_alpha1i = (P1_1_SEP - alpha1*P1_2_SEP);

f_alpha0o = (P2_0_SEP + alpha1*P2_1_SEP + 0.5*alpha1*alpha1*P2_2_SEP);
f_alpha0i = (P2_1_SEP + alpha1*P2_2_SEP);

f_alpha2o = (P2_0_SEP - alpha2*P2_1_SEP + 0.5*(alpha2*alpha2)*P2_2_SEP);
f_alpha2i = (P2_1_SEP - alpha2*P2_2_SEP);

f_alpha3o = (P3_0_SEP + alpha3*(s3*s3*P3_0_SEP + 2.0*s3*P3_1_SEP + P3_2_SEP));
f_alpha3i = (P3_1_SEP + alpha3*(s3*s3*P3_1_SEP + 2.0*s3*P3_2_SEP));

% % % XB ODEs
 dSL_SEP = (sigmaf_SEP - sigmapas_SEP - sigmaact_SEP)/eta;% 
 
P0_SEP = 1 - N_SEP - P1_0_SEP - P2_0_SEP - P3_0_SEP;  % DAB 10/8/2019
dXdT(21) = ka*P0_SEP*U_NR_SEP*N_overlap_SEP - kd_SEP*P1_0_SEP - k1_SEP*f_alpha1o + km1*f_alpha0o; % DAB 10/8/2019
dXdT(22) = dSL_SEP*P1_0_SEP - kd_SEP*P1_1_SEP - k1_SEP*f_alpha1i + km1*f_alpha0i;
dXdT(23) = 2*dSL_SEP*P1_1_SEP - kd_SEP*P1_2_SEP - k1_SEP*P1_2_SEP + km1*P2_2_SEP;

dXdT(24) = -km1*f_alpha0o - k2*f_alpha2o + km2_SEP*P3_0_SEP + k1_SEP*f_alpha1o;
dXdT(25) = dSL_SEP*P2_0_SEP - k1_SEP*f_alpha0i - k2*f_alpha2i + km2_SEP*P3_1_SEP + k1_SEP*f_alpha1i;
dXdT(26) = 2*dSL_SEP*P2_1_SEP - km1*P2_2_SEP       - k2*P2_2_SEP + km2_SEP*P3_2_SEP + k1_SEP*P1_2_SEP;

dXdT(27) = +k2*f_alpha2o - km2_SEP*P3_0_SEP - k3_SEP*f_alpha3o;
dXdT(28) = dSL_SEP*P3_0_SEP + k2*f_alpha2i - km2_SEP*P3_1_SEP - k3_SEP*f_alpha3i;
dXdT(29) = 2*dSL_SEP*P3_1_SEP + k2*P2_2_SEP       - km2_SEP*P3_2_SEP - k3_SEP*P3_2_SEP;

U_SR_SEP = 1 - U_NR_SEP;
Jon = k_on*Ca_i*N_SEP*(1 + K_coop*(1 - N_SEP));
Joff = k_off*P0_SEP*(1 + K_coop*N_SEP);
dXdT(30) = - Jon + Joff; 
dXdT(31) = ksr*(1 + kforce * sigmaact_SEP) * U_SR_SEP - kmsr*U_NR_SEP ; % DAB 10/8

%% Myofiber Mechanics: SL_RV

% Calculations for stretch-senstive rates    
f_alpha1o = (P1_0_RV - alpha1*P1_1_RV + 0.5*(alpha1*alpha1)*P1_2_RV);
f_alpha1i = (P1_1_RV - alpha1*P1_2_RV);

f_alpha0o = (P2_0_RV + alpha1*P2_1_RV + 0.5*alpha1*alpha1*P2_2_RV);
f_alpha0i = (P2_1_RV + alpha1*P2_2_RV);

f_alpha2o = (P2_0_RV - alpha2*P2_1_RV + 0.5*(alpha2*alpha2)*P2_2_RV);
f_alpha2i = (P2_1_RV - alpha2*P2_2_RV);

f_alpha3o = (P3_0_RV + alpha3*(s3*s3*P3_0_RV + 2.0*s3*P3_1_RV + P3_2_RV));
f_alpha3i = (P3_1_RV + alpha3*(s3*s3*P3_1_RV + 2.0*s3*P3_2_RV));

% % XB ODEs
 dSL_RV = (sigmaf_RV - sigmapas_RV - sigmaact_RV)/eta;% 

P0_RV = 1.0 - N_RV - P1_0_RV - P2_0_RV - P3_0_RV;  % DAB 10/8/2019
dXdT(32) = ka*P0_RV*U_NR_RV*N_overlap_RV - kd_RV*P1_0_RV - k1_RV*f_alpha1o + km1*f_alpha0o;  % DAB 10/8/2019
dXdT(33) = dSL_RV*P1_0_RV - kd_RV*P1_1_RV - k1_RV*f_alpha1i + km1*f_alpha0i;
dXdT(34) = 2*dSL_RV*P1_1_RV - kd_RV*P1_2_RV - k1_RV*P1_2_RV + km1*P2_2_RV;

dXdT(35) = -km1*f_alpha0o - k2*f_alpha2o + km2_RV*P3_0_RV + k1_RV*f_alpha1o;
dXdT(36) = dSL_RV*P2_0_RV - k1_RV*f_alpha0i - k2*f_alpha2i + km2_RV*P3_1_RV + k1_RV*f_alpha1i;
dXdT(37) = 2*dSL_RV*P2_1_RV - km1*P2_2_RV       - k2*P2_2_RV + km2_RV*P3_2_RV + k1_RV*P1_2_RV;

dXdT(38) = +k2*f_alpha2o - km2_RV*P3_0_RV - k3_RV*f_alpha3o;
dXdT(39) = dSL_RV*P3_0_RV + k2*f_alpha2i - km2_RV*P3_1_RV - k3_RV*f_alpha3i;
dXdT(40) = 2*dSL_RV*P3_1_RV + k2*P2_2_RV       - km2_RV*P3_2_RV - k3_RV*P3_2_RV;


U_SR_RV = 1.0 - U_NR_RV;
Jon = k_on*Ca_i*N_RV*(1 + K_coop*(1 - N_RV));
Joff = k_off*P0_RV*(1 + K_coop*N_RV);
dXdT(41) = - Jon + Joff; 
dXdT(42) = ksr*(1 + kforce * sigmaact_RV) * U_SR_RV - kmsr * U_NR_RV ; % DAB 10/8

% Myofiber Mechanics
dXdT(5) = dSL_LV;%(Kse*(SLo_LV - SL_LV) - sigmapas_LV - sigmaact_LV)/eta;
dXdT(6) = dSL_SEP;%(Kse*(SLo_SEP - SL_SEP) - sigmapas_SEP - sigmaact2)/eta;
dXdT(7) = dSL_RV;%(Kse*(SLo_RV - SL_RV) - sigmapas_RV - sigmaact3)/eta;

% Lumped circulation model variables
dXdT(8) = QIN_LV - QOUT_LV; % V_LV
dXdT(9) = QIN_RV - QOUT_RV; % V_RV
dXdT(43) = (P_SA - P_SV)/R_SA - QIN_RV;  % V_SV
dXdT(44) = (P_PA - P_PV)/R_PA - QIN_LV;  % V_PV
dXdT(45) = Q_Ao - (P_SA - P_SV)/R_SA; % V_SA 
dXdT(46) = QOUT_RV - (P_PA - P_PV)/R_PA; % V_PA 
dXdT(47) = QOUT_LV - Q_Ao; % V_Ao

dXdT = dXdT(:);


