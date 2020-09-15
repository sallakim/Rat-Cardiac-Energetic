% Run rat mechano-energetic model for all rats
% DriverGeneticAlgorithm_MechanoEnegergetic
clear; 

%% specifying what simulation should be run, all rats or just a single rat.
rat_number = 10; % assgining the rat number

% flag_swap_metabolite is set to 1 for the simulations with the replaced
% metabolite. 
flag_swap_metabolite = 0;

% Reading the adujstable scale variables for each individual rats

adjvar_all_rest = xlsread('Adjustable_parameters_table_rest.xlsx','B2:J21'); % rat 9 is mean sham rat
adjvar_all_swap = xlsread('Adjustable_parameters_table_swap.xlsx','B2:J21'); % rat 9 is mean sham rat

adjvar_all = adjvar_all_rest;
if flag_swap_metabolite ==1
 adjvar_all = adjvar_all_swap;  
end
rate_of_XB_turnover_mean_sham = 5.0147;

if rat_number<=9
    shamRat = 1;
    delta_p = 0;
else
    shamRat = 0;
end
%% adjustable variables 
% adjvar = [Reference area LV , Reference area Septal,  Reference area RV, ksr& kforce, Blood volume, R_SA, R_TAC, ATP_tune_Coeff]

adjvar = adjvar_all(rat_number,:);

%% Read the experimental data for SHAM and TAC rats from the excel file 
data = xlsread('data1.xlsx','A3:W23');
BW  = data(rat_number , 1); % g
LVW = data(rat_number , 2); % mg
RVW = data(rat_number , 3); % mg
LW = data(rat_number , 5); % mg
HR  = data(rat_number , 6); % beats/min

edLV_target = data(rat_number , 13); % uL
esLV_target = data(rat_number , 14); % uL

SV_LV_target = edLV_target - esLV_target;
CO_target = SV_LV_target / 1000 * HR;
EF_LV_target = SV_LV_target / edLV_target * 100;
TAN = data(rat_number,16)/1000; % mole/L cell
CRtot = data(rat_number,18)/1000; % mole/L cell
Ox_capacity = data(rat_number,21)/data(9,21); 
Ox_capacity_sham = 1; 
Ox_capacity_TAC = data(20,21)/data(9,21);
% Average sham
TAN_sham = data(9,16)/1000; % mole/L cell
CRtot_sham = data(9,18)/1000; % mole/L cell
TAN_TAC =  data(20,16)/1000; % mole/L cell
CRtot_TAC = data(20,18)/1000; % mole/L cell

Ao = 10.260e-3; % (M per liter cell)
Co = 43.007e-3; % (M per liter cell)
Po = 35.446e-3; % (M per liter cell)
% Default settings for meant SHAM)

TEP = Po - (0.283e-3)*(Ao-TAN)/(0.082e-3); % (M per liter cell)
TEP_sham = Po - (0.283e-3)*(Ao-TAN_sham)/(0.082e-3); % (M per liter cell)
TEP_TAC = Po - (0.283e-3)*(Ao-TAN_TAC)/(0.082e-3); % (M per liter cell)

if flag_swap_metabolite ==1
    if shamRat == 0
    TAN = TAN_sham;
    CRtot = CRtot_sham;
    TEP = TEP_sham;
    Ox_capacity = Ox_capacity_sham;
    else 
    TAN = TAN_TAC;
    CRtot = CRtot_TAC;
    TEP = TEP_TAC;
    Ox_capacity = Ox_capacity_TAC;
    end
end

if shamRat == 0
    preV = data(rat_number , 11); % mm/s
    postV = data(rat_number , 12); % mm/s
    postV = postV/1000; %m/s
    preV = preV/1000; %m/s
    rho_blood = 1060; % kg/m^3
    delta_p = 0.5*(postV^2-preV^2)*rho_blood; % Pa
    delta_p = 0.0075*delta_p; % mmHg
    if delta_p > 60 % replacing the delta_p with Max delta_P dor the animals without measured post TAC velocity
        delta_p = 31.48;
    end
end

CO_target = 95; % ml/min
% MAP_target = 93.33; %mmHg taregt mean arterial pressure based on MAP = DBP +[1/3(SBP - DBP)];

R_TAC = adjvar(8);

tune_ATPase_LV =  adjvar(9)* (1/ 0.6801) *1.0e-3; % ATP hydrolysis rate: M / s / (liter cytosol)

Amref_LV  = adjvar(1); % LV midwall reference surface area, cm^2
Amref_SEP = adjvar(2); % SEP midwall reference surface area, cm^2
Amref_RV  = adjvar(3); % RV midwall reference surface area, cm^2

% Assign initial condtion for LV and RV

V_LV  = edLV_target/1000 + 0.2;% intial value for V_LV and V_RV assumed to be equal to edLV_target
V_RV  = edLV_target/1000 + 0.2;%


Vw_LV = (LVW*2/3)/1000/1.05;
Vw_SEP =(LVW/3)/1000/1.05;
Vw_RV = RVW/1000/1.05;

%% Run the energetics model to get the metabolite concentrations

    energtics_output  = EnergeticsModelScript(TAN, CRtot, TEP, Ox_capacity, tune_ATPase_LV);
   
    MgATP_cyto = energtics_output(1);
    MgADP_cyto = energtics_output(2);
    
    fPi_cytoplasm = energtics_output(3); 
    MVO2_tissue = energtics_output(5);
    dGrATPase = energtics_output(6); % kJ / mol

    PCrATP =  energtics_output(7); %(unitless)
    ATP_cyto = energtics_output(8); %cyto [ATP]
    ADP_cyto = energtics_output(9); %cyto [ADP]
    Pi_cyto = energtics_output(10)*1000;%cyto [Pi]
   

%% Run cardiovascular mechanics model
Lsref = 1.9;
k3      = 144.5586; % transition A3 to P rate constant, 1/sec
K_T = 0.4897; 
K_D = 0.194;% Used the values from Tewari etal JMCC (9/5 BM)
alpha3  =      0.1*59.3; % Stretch sensing parameter for k3, 1/um
s3      = 9.9e-3;  % Strain at which k3 is minimum, um

%% Extract Ca coeficient based on the Ca data for diferent simulation 
para_fitted_Ca = [2	3	4	5	6	7	8	9	10;
0.0838	0.1306	0.1802	0.2557	0.3099	0.3613	0.408	0.4539	0.4602;
0.7513	0.8756	1.0274	1.4988	1.6107	1.6741	1.7902	2.1398	1.9832;
2.237	2.0486	1.948	1.852	1.6932	1.6773	1.5988	1.4952	1.4524;
0.1865	0.1815	0.1709	0.1693	0.161	0.1661	0.1425	0.1458	0.1222];
freq_all = para_fitted_Ca(1,:);
A_HR_pchip = pchip(freq_all,para_fitted_Ca(2,:));
A_HR = ppval(A_HR_pchip,HR/60);
B_HR_pchip = pchip(freq_all,para_fitted_Ca(3,:));
B_HR = ppval(B_HR_pchip,HR/60);
C_HR_pchip = pchip(freq_all,para_fitted_Ca(4,:));
C_HR = ppval(C_HR_pchip,HR/60);
Ca0_HR_pchip = pchip(freq_all,para_fitted_Ca(5,:));
Ca0_HR = ppval(Ca0_HR_pchip,HR/60);

stim_period = 1/(HR/60);
M = speye(47);
M(1,1) = 0;
M(2,2) = 0;
M(3,3) = 0;
M(4,4) = 0; 
input = [CO_target stim_period Vw_LV Vw_SEP Vw_RV R_TAC MgATP_cyto MgADP_cyto Pi_cyto MgATP_cyto MgADP_cyto Pi_cyto MgATP_cyto MgADP_cyto Pi_cyto A_HR B_HR C_HR Ca0_HR Amref_LV Amref_SEP Amref_RV flag_swap_metabolite];
%  options = odeset('Mass',M,'RelTol',1e-5,'AbsTol',1e-5,'MaxStep',stim_period/50,'OutputFcn',@odeprint);
options = odeset('Mass',M,'MassSingular','yes','RelTol',1e-7,'AbsTol',1e-7,'MaxStep',stim_period/200);

xm_LV   = -0.60;
xm_SEP  = 0.40;
xm_RV   = 1.0;
ym    = 0.50;

SL_LV  = 2.2;
SL_SEP = 2.2;
SL_RV  = 2.2;

V_SA = adjvar(6)* 3.0;
V_SV = adjvar(6)* 4.80;
V_PA = adjvar(6)* 0.5;
V_PV = adjvar(6)* 1.0; 
V_Ao = adjvar(6)* 1.0;

P1_0_LV = 0; % 0th moment state A1, LV
P1_1_LV = 0; % 1st moment state A1, LV
P1_2_LV = 0; % 2nd moment state A1, LV
P2_0_LV = 0; % 0th moment state A2, LV
P2_1_LV = 0; % 1st moment state A2, LV
P2_2_LV = 0; % 2nd moment state A2, LV
P3_0_LV = 0; % 0th moment state A3, LV
P3_1_LV = 0; % 1st moment state A3, LV
P3_2_LV = 0; % 2nd moment state A3, LV
N_LV = 1;
U_NR_LV = 0;
P1_0_RV = 0; % 0th moment state A1, LV
P1_1_RV = 0; % 1st moment state A1, LV
P1_2_RV = 0; % 2nd moment state A1, LV
P2_0_RV = 0; % 0th moment state A2, LV
P2_1_RV = 0; % 1st moment state A2, LV
P2_2_RV = 0; % 2nd moment state A2, LV
P3_0_RV = 0; % 0th moment state A3, LV
P3_1_RV = 0; % 1st moment state A3, LV
P3_2_RV = 0; % 2nd moment state A3, LV
N_RV = 1;
U_NR_RV = 0;
P1_0_SEP = 0; % 0th moment state A1, LV
P1_1_SEP = 0; % 1st moment state A1, LV
P1_2_SEP = 0; % 2nd moment state A1, LV
P2_0_SEP = 0; % 0th moment state A2, LV
P2_1_SEP= 0; % 1st moment state A2, LV
P2_2_SEP = 0; % 2nd moment state A2, LV
P3_0_SEP = 0; % 0th moment state A3, LV
P3_1_SEP = 0; % 1st moment state A3, LV
P3_2_SEP = 0; % 2nd moment state A3, LV
N_SEP = 1;
U_NR_SEP = 0;


init = [xm_LV ,xm_SEP ,xm_RV ,ym , SL_LV, SL_SEP, SL_RV, V_LV, V_RV, ...
       P1_0_LV, P1_1_LV, P1_2_LV ,P2_0_LV, P2_1_LV, P2_2_LV, P3_0_LV, P3_1_LV, P3_2_LV, N_LV, U_NR_LV,...
       P1_0_SEP,P1_1_SEP,P1_2_SEP,P2_0_SEP,P2_1_SEP,P2_2_SEP,P3_0_SEP,P3_1_SEP,P3_2_SEP,N_SEP,U_NR_SEP,...
       P1_0_RV, P1_1_RV, P1_2_RV, P2_0_RV, P2_1_RV, P2_2_RV, P3_0_RV, P3_1_RV, P3_2_RV, N_RV, U_NR_RV,...
       V_SV, V_PV ,V_SA ,V_PA, V_Ao]';

 
 opts = optimset('MaxFunEvals',10000,'MaxIter',1000);
%     TrisegEquations(init(1:4),Vw_LV,Vw_SEP,Vw_RV,SL_LV,SL_SEP,SL_RV,V_LV,V_RV,Amref_LV,Amref_SEP,Amref_RV);
   x_triseg = fsolve(@TrisegEquations,init(1:4),opts,Vw_LV,Vw_SEP,Vw_RV,SL_LV,SL_SEP,SL_RV,V_LV,V_RV,Amref_LV,Amref_SEP,Amref_RV);
   init(1:4) = x_triseg;

    [ts,ys] = ode15s(@dXdT_cardiovascular_mechanics,[0 120*stim_period],init,options,adjvar,input);
% V_LV   = ys(:,8); % volume LV, mL
% V_RV   = ys(:,9); % volume RV, mL
% V_SV   = ys(:,43); % volume of systemic veins
% V_PV   = ys(:,44); % volume of pulmonary veins
% V_SA   = ys(:,45); % volume of systemic arterys
% V_PA   = ys(:,46); % volume of pulmonary arterys
% V_Ao   = ys(:,47); % volume of proximal aorta
% V_T = V_LV + V_RV + V_SV + V_PV + V_SA + V_PA + V_Ao;
% 
% P_PV = V_PV/C_PV;
% P_PA = V_PA/C_PA;
% 
%  figure(1); plot(ts,V_LV,ts,ys(:,11)); title('ventricular volumes')
% % figure(2); plot(ts,(V_PV + V_PA)./V_T);
% figure(3); plot(ts,P_PA,ts,P_PV); title('pulmonary pressures');
init = ys(end,:);
[t,Y] = ode15s(@dXdT_cardiovascular_mechanics,[0 1*stim_period],init,options,adjvar,input);

% % Assignig the solution of the ODE's to the variables

xm_LV  = Y(:,1); % LV heart geometry variable, cm
% xm_SEP = Y(:,2); % septum heart geometry variable, cm
xm_RV  = Y(:,3); % RV heart geometry variable, cm
ym     = Y(:,4); % Heart geometry variable, cm
SL_LV  = Y(:,5); % sarcomere length, LV, micron
SL_SEP = Y(:,6); % sarcomere length, septum, micron
SL_RV  = Y(:,7); % sarcomere length, RV, micron
V_LV   = Y(:,8); % volume LV, mL
V_RV   = Y(:,9); % volume RV, mL
% 
% P1_0_LV = Y(:,10); % 0th moment state A1, LV
% P1_1_LV = Y(:,11); % 1st moment state A1, LV
% P1_2_LV = Y(:,12); % 2nd moment state A1, LV
% P2_0_LV = Y(:,13); % 0th moment state A2, LV
% P2_1_LV = Y(:,14); % 1st moment state A2, LV
% P2_2_LV = Y(:,15); % 2nd moment state A2, LV
P3_0_LV = Y(:,16); % 0th moment state A3, LV
P3_1_LV = Y(:,17); % 1st moment state A3, LV
P3_2_LV = Y(:,18); % 2nd moment state A3, LV
% N_LV    = Y(:,19); % non-permissive fraction LV
% U_NR_LV = Y(:,20); % U_NR represents the Non relaxed state

% P1_0_SEP = Y(:,21); % 0th moment state A1, SEP
% P1_1_SEP = Y(:,22); % 1st moment state A1, SEP
% P1_2_SEP = Y(:,23); % 2nd moment state A1, SEP
% P2_0_SEP = Y(:,24); % 0th moment state A2, SEP
% P2_1_SEP = Y(:,25); % 1st moment state A2, SEP
% P2_2_SEP = Y(:,26); % 2nd moment state A2, SEP
P3_0_SEP = Y(:,27); % 0th moment state A3, SEP
P3_1_SEP = Y(:,28); % 1st moment state A3, SEP
P3_2_SEP = Y(:,29); % 2nd moment state A3, SEP
% N_SEP    = Y(:,30); % nonpermissive fraction SEP
% U_NR_SEP = Y(:,31); % U_NR represents the Non relaxed state

% P1_0_RV = Y(:,32); % 0th moment state A1, RV
% P1_1_RV = Y(:,33); % 1st moment state A1, RV
% P1_2_RV = Y(:,34); % 2nd moment state A1, RV
% P2_0_RV = Y(:,35); % 0th moment state A2, RV
% P2_1_RV = Y(:,36); % 1st moment state A2, RV
% P2_2_RV = Y(:,37); % 2nd moment state A2, RV
% P3_0_RV = Y(:,38); % 0th moment state A3, RV
% P3_1_RV = Y(:,39); % 1st moment state A3, RV
% P3_2_RV = Y(:,40); % 2nd moment state A3, RV
% N_RV    = Y(:,41); % nonpermissive fraction RV
% U_NR_RV = Y(:,42); % U_NR represents the Non relaxed state

V_SV   = Y(:,43); % volume of systemic veins
V_PV   = Y(:,44); % volume of pulmonary veins
V_SA   = Y(:,45); % volume of systemic arterys
V_PA   = Y(:,46); % volume of pulmonary arterys
V_Ao   = Y(:,47); % volume of proximal aorta
% V_T = V_LV + V_RV + V_SV + V_PV + V_SA + V_PA + V_Ao;

% PlotTriSeg(xm_LV,xm_SEP,xm_RV,ym,t)
% Lumped circulatory parameters
C_Ao = 0.0022045;  % Proximal aortic compliance, mL/mmHg
C_SA = 0.0077157; % Systemic arterial compliance, mL/mmHg
C_SV = 2.5; % Systemic venous compliance, mL/mmHg  DAB 10/7/2018
C_PV = 0.25; % Pulmonary venous compliance, mL/mmHg
C_PA = 0.013778; % Pulmonary arterial compliance, mL/mmHg
R_Ao   = 2.5; % resistance of aorta , mmHg*sec/mL
R_SA   = adjvar(7);% mmHg*sec/mL; % Systemic vasculature resistance, mmHg*sec/mL
R_vlv  = 0.05; %  valve resistance, mmHg*sec/mL
R_AV   = R_vlv + R_TAC; % resistance across aortic valve
R_tAo  = 0.5;
R_tSA  = 4;
Kse    = 50000; % series element elastance, mmHg/micron (Changed to match the value in Tewari's code) (9/5 BM)

%  Pulmonary Pressures
P_PV = V_PV/C_PV;
P_SV = V_SV/C_SV;
P_PA = V_PA/C_PA;
P_SA = V_SA/C_SA;

Am_LV = pi*(xm_LV.^2 + ym.^2);
% Am_SEP = pi*(xm_SEP.^2 + ym.^2);
Am_RV = pi*(xm_RV.^2 + ym.^2);
Cm_LV = 2*xm_LV./(xm_LV.^2 + ym.^2);
% Cm_SEP = 2*xm_SEP./(xm_SEP.^2 + ym.^2);
Cm_RV = 2*xm_RV./(xm_RV.^2 + ym.^2);
z_LV = 3*Cm_LV.*Vw_LV./(2*Am_LV);
% z_SEP = 3*Cm_SEP.*Vw_SEP./(2*Am_SEP);
z_RV = 3*Cm_RV.*Vw_RV./(2*Am_RV);

epsf_LV = (1/2)*log(Am_LV./Amref_LV) - (1/12)*z_LV.^2 - 0.019*z_LV.^4;
% epsf_SEP = (1/2)*log(Am_SEP./Amref_SEP) - (1/12)*z_SEP.^2 - 0.019*z_SEP.^4;
epsf_RV = (1/2)*log(Am_RV./Amref_RV) - (1/12)*z_RV.^2 - 0.019*z_RV.^4;
SLo_LV = Lsref*exp(epsf_LV);
% SLo_SEP = Lsref*exp(epsf_SEP);
SLo_RV = Lsref*exp(epsf_RV);

% % Total forces
sigmaf_LV = -Kse*(SL_LV - SLo_LV);
% sigmaf_SEP = -Kse*(SL_SEP - SLo_SEP);
sigmaf_RV = -Kse*(SL_RV - SLo_RV);

% % equilibrium of forces at junction circle
Tm_LV = (Vw_LV.*sigmaf_LV./(2*Am_LV)).*(1 + z_LV.^2/3 + z_LV.^4/5);
% Tm_SEP = (Vw_SEP.*sigmaf_SEP./(2*Am_SEP)).*(1 + z_SEP.^2/3 + z_SEP.^4/5);
Tm_RV = (Vw_RV.*sigmaf_RV./(2*Am_RV)).*(1 + z_RV.^2/3 + z_RV.^4/5);
sinalpha_LV = 2*xm_LV.*ym./(xm_LV.^2 + ym.^2);
% sinalpha_SEP = 2*xm_SEP.*ym./(xm_SEP.^2 + ym.^2);
sinalpha_RV = 2*xm_RV.*ym./(xm_RV.^2 + ym.^2);
% cosalpha_LV = (-xm_LV.^2 + ym.^2)./(xm_LV.^2 + ym.^2);
% cosalpha_SEP = (-xm_SEP.^2 + ym.^2)./(xm_SEP.^2 + ym.^2);
% cosalpha_RV = (-xm_RV.^2 + ym.^2)./(xm_RV.^2 + ym.^2);
Tx_LV = Tm_LV.*sinalpha_LV;
% Tx_SEP = Tm_SEP.*sinalpha_SEP;
Tx_RV = Tm_RV.*sinalpha_RV;
% Ty_LV = Tm_LV.*cosalpha_LV;
% Ty_SEP = Tm_SEP.*cosalpha_SEP;
% Ty_RV= Tm_RV.*cosalpha_RV;

% % ventricular pressure
ptrans_LV = 2*Tx_LV./ym;
ptrans_RV = 2*Tx_RV./ym;
P_LV = -ptrans_LV;
P_RV = ptrans_RV;
% Ao valves closed equations
P_Ao_closed = (C_SA*R_Ao*R_SA*V_Ao + C_SA*R_Ao*R_tSA*V_Ao + C_SA*R_SA*R_tSA*V_Ao + C_Ao*R_SA*R_tAo*V_SA + C_Ao*C_SA*P_SV*R_tSA*R_tAo)/(C_Ao*C_SA*(R_Ao*R_SA + R_Ao*R_tSA + R_SA*R_tSA + R_SA*R_tAo + R_tSA*R_tAo));
% P_SA_closed = (C_Ao*R_Ao*R_SA*V_SA + C_SA*R_SA*R_tSA*V_Ao + C_Ao*R_SA*R_tAo*V_SA + C_Ao*C_SA*P_SV*R_Ao*R_tSA + C_Ao*C_SA*P_SV*R_tSA*R_tAo)/(C_Ao*C_SA*(R_Ao*R_SA + R_Ao*R_tSA + R_SA*R_tSA + R_SA*R_tAo + R_tSA*R_tAo));
% Ao valve open equations 
P_Ao_open = (C_SA*R_Ao*R_SA*R_AV*V_Ao + C_SA*R_Ao*R_tSA*R_AV*V_Ao + C_SA*R_SA*R_tSA*R_AV*V_Ao + C_Ao*R_SA*R_tAo*R_AV*V_SA + C_Ao*C_SA*P_LV*R_Ao*R_SA*R_tAo + C_Ao*C_SA*P_LV*R_Ao*R_tSA*R_tAo + C_Ao*C_SA*P_LV*R_SA*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_tSA*R_tAo*R_AV)/(C_Ao*C_SA*(R_Ao*R_SA*R_tAo + R_Ao*R_SA*R_AV + R_Ao*R_tSA*R_tAo + R_Ao*R_tSA*R_AV + R_SA*R_tSA*R_tAo + R_SA*R_tSA*R_AV + R_SA*R_tAo*R_AV + R_tSA*R_tAo*R_AV));
% P_SA_open = (C_Ao*R_Ao*R_SA*R_tAo*V_SA + C_Ao*R_Ao*R_SA*R_AV*V_SA + C_SA*R_SA*R_tSA*R_AV*V_Ao + C_Ao*R_SA*R_tAo*R_AV*V_SA + C_Ao*C_SA*P_SV*R_Ao*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_Ao*R_tSA*R_AV + C_Ao*C_SA*P_LV*R_SA*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_tSA*R_tAo*R_AV)/(C_Ao*C_SA*(R_Ao*R_SA*R_tAo + R_Ao*R_SA*R_AV + R_Ao*R_tSA*R_tAo + R_Ao*R_tSA*R_AV + R_SA*R_tSA*R_tAo + R_SA*R_tSA*R_AV + R_SA*R_tAo*R_AV + R_tSA*R_tAo*R_AV));
% QOUT_LV = -(C_SA*R_Ao*R_SA*V_Ao + C_SA*R_Ao*R_tSA*V_Ao + C_SA*R_SA*R_tSA*V_Ao + C_Ao*R_SA*R_tAo*V_SA - C_Ao*C_SA*P_LV*R_Ao*R_SA - C_Ao*C_SA*P_LV*R_Ao*R_tSA - C_Ao*C_SA*P_LV*R_SA*R_tSA - C_Ao*C_SA*P_LV*R_SA*R_tAo - C_Ao*C_SA*P_LV*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_tSA*R_tAo)/(C_Ao*C_SA*(R_Ao*R_SA*R_tAo + R_Ao*R_SA*R_AV + R_Ao*R_tSA*R_tAo + R_Ao*R_tSA*R_AV + R_SA*R_tSA*R_tAo + R_SA*R_tSA*R_AV + R_SA*R_tAo*R_AV + R_tSA*R_tAo*R_AV));

P_Ao = P_Ao_open.*(P_LV>P_Ao_open) + ...
       P_Ao_closed.*(P_LV<=P_Ao_open);
% QOUT_LV = QOUT_LV.*(P_LV>P_Ao_open);
% 
% % CO_sim =  (max(V_LV)-min(V_LV))*HR;
% % 
% % CO_mean = mean(QOUT_LV);

SV_LV_sim = max(1e3*V_LV) - min(1e3*V_LV);
EF_LV_sim = SV_LV_sim/max(1e3*V_LV) * 100;
    
%% calculate the error
edLV_sim =  max(1e3*V_LV);
esLV_sim =  min(1e3*V_LV);
edRV_sim =  max(1e3*V_RV);
esRV_sim =  min(1e3*V_RV);

g2_LV = (MgATP_cyto/K_T)/(1.0 + MgATP_cyto/K_T + MgADP_cyto/K_D);
k3_LV = k3*g2_LV;%
g2_SEP = (MgATP_cyto/K_T)/(1.0 + MgATP_cyto/K_T + MgADP_cyto/K_D);
k3_SEP = k3*g2_SEP;%
f_alpha3o_LV  = (P3_0_LV + alpha3*(s3*s3*P3_0_LV + 2.0*s3*P3_1_LV + P3_2_LV)); 
f_alpha3o_SEP = (P3_0_SEP + alpha3*(s3*s3*P3_0_SEP + 2.0*s3*P3_1_SEP + P3_2_SEP)); 

% detachment rates
ti = 0:0.00001:stim_period;
MAP = mean(interp1(t,P_Ao,ti));
% MeanSVP = mean(interp1(t,P_SV,ti))
% MeanPAP = mean(interp1(t,P_PA,ti))
% MeanPVP = mean(interp1(t,P_PV,ti))

r_LV  = interp1(t,k3_LV*f_alpha3o_LV,ti);
r_SEP = interp1(t,k3_SEP*f_alpha3o_SEP,ti);

% LV X-bridge turnover rate
Vw_LV_W = (2/3)*LVW/1000;
Vw_SEP_W= (1/3)*LVW/1000;
rate_of_XB_turnover_ave = (Vw_LV_W*mean(r_LV) + Vw_SEP_W*mean(r_SEP))/(Vw_LV_W + Vw_SEP_W);

% unit convert to oxygen consumption
ATP_ase_mechannics_Averge_LV_SEP = (1.3/rate_of_XB_turnover_mean_sham)*rate_of_XB_turnover_ave; %  1.31 Kstiff - ATP hydrolized (mmol/s/(L cell)) per X-bridge turnover rate in LV
%% Calculation of the error
SL_LV_MAX = max(SL_LV);
SL_SEP_MAX = max(SL_SEP);
SL_RV_MAX = max(SL_RV);

x_ATPase_err = (adjvar(9)- ATP_ase_mechannics_Averge_LV_SEP)/ATP_ase_mechannics_Averge_LV_SEP;

%% Plotting MechanoEnergetics Results
%% Plotting the pressure time course figure
ylimMax = max(P_LV)+ 20;
 
figure(1); clf;
hold on;
title(['rat number = ', num2str(rat_number)]);
if shamRat == 1
plot(t,P_LV,'k',t,P_Ao,'r-',t,P_SA,'g','linewidth',2);

else
plot(t,P_LV,'k',t,P_Ao,'r-',t,P_SA,'g','linewidth',2);
hold on
plot([0 1/(HR/60)],[max(P_Ao) max(P_Ao)],'k--',[0 1/(HR/60)],[max(P_Ao)+4*delta_p max(P_Ao)+4*delta_p],'k--','linewidth',2); hold off
       
end

legend('P_{LV}','P_{Ao}','P_{SA}');

set(gca,'Fontsize',20,'linewidth',2)
box off;
 legend('boxoff') 
set(gca,'Ylim',[0 ylimMax]);
ylabel('P (mmHg)','fontsize',20,'FontWeight','bold')

xlabel('t (sec)','fontsize',20,'FontWeight','bold')
set(gca,'Fontsize',20,'linewidth',2)

%% Plotting the PV loop figure
figure(2); clf;
 
plot(V_LV,P_LV,'k',[edLV_target edLV_target]./1000,[0 ylimMax],'k--',[esLV_target esLV_target]./1000,[0 ylimMax],'k--','linewidth',2); 
set(gca,'Xlim',[0 0.70]);
set(gca,'Ylim',[0 ylimMax]);

ylabel('P (mmHg)','fontsize',20,'FontWeight','bold')
xlabel('V (mL)','fontsize',20,'FontWeight','bold')
set(gca,'Fontsize',20,'linewidth',2)
box off;

txt_ATPase = ['x ATPase =', num2str(ATP_ase_mechannics_Averge_LV_SEP)];
disp(txt_ATPase)
txt_XB_rate = ['rate turn XB =', num2str(rate_of_XB_turnover_ave)];
disp(txt_XB_rate)
disp(MAP)
disp(EF_LV_sim) 
% % save data to excel file 
% CO_sim =  (max(V_LV)-min(V_LV))*HR;
% k1sr = adjvar(4);
% kforce = adjvar(5);  %unit converted by kPa to mmHg conversion rate
% 
% result_table(rat_number,:) = [rat_number, Pi_cyto, MgATP_cyto, MgADP_cyto, EF_LV_sim, MAP, P_LV(end), dGrATPase, CO_sim];
% table_S1_3 (rat_number,:) = [TAN CRtot TEP Ox_capacity adjvar(9)]
% table_S1_4 (rat_number,:) = [MgATP_cyto MgADP_cyto fPi_cytoplasm MVO2_tissue dGrATPase PCrATP ATP_cyto ADP_cyto Pi_cyto];
% table_S5_1 (rat_number,:) = [HR Vw_LV Vw_SEP Vw_RV rate_of_XB_turnover_ave edLV_sim esLV_sim];
% adjustable_parameter_table(rat_number,:) =[rat_number, Amref_LV, Amref_SEP, Amref_RV, k1sr,kforce, adjvar(6),R_SA,R_TAC ,adjvar(9)];
% 
