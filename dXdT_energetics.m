% function [f,J,] = dXdT(t,x,T,BX,K_BX,par)
%
%
%
% Output parameters:
%   f     time derivatives of the model
%   J     flux
%
% Mandatory input parameters:
%   t     time
%   x     state variables at time t
%   var(1) temperature in degreesCelcius
%   var(2) mitochondrial oxidative capacity (unitless, relative to control)
%   BX    buffer sizes
%   K_BX  proton buffer dissociation constants ( matrix cytoplasm im )
%   par   parameter vector for the free parameters
%   
%
% State Variables: (mitochondrial matrix metabolic, inter-membrane space,
% and cytosolic metabolite state variables)
% [NAD_matrix , NADH_matrix , ADP_matrix , ATP_matrix , Pi_matrix , coQ_matrix , coQH2_matrix , ATP_c , ADP_c , Pi_c , phosphocreatine_c , creatine_c , AMP_c , H2O2aq_matrix , SOaq_matrix , cytocox_im , cytocred_im , O2aq_matrix , Pi_im , ADP_im , ATP_im , AMP_im ]
%
% Parameters: mitochondrial energetics 
% [x_DH , x_SDH , x_ATPase , k1_F1F0 , ETC1Activity , ETC3_activity , ETC4_activity , x_HLE , x_PIH , x_ANT , k1_KH , x_AMPPERM , x_ADPPERM , x_ATPPERM , x_PIPERM , x_HPERM , x_KPERM , x_MPERM ]

function [f,J] = dXdT_energetics(~,x,var,BX,K_BX,par,flag)

% mito membrane area per cell volume micron^{-1};
gamma  = 5.99;

%% LIST OF STATE VARIABLES (table 1 from Marzban 2020)
% 1 NAD_matrix
% 2 NADH_matrix
% 3 ADP_matrix
% 4 ATP_matrix
% 5 Pi_matrix
% 6 coQ_matrix
% 7 coQH2_matrix
% 8 ATP_c
% 9 ADP_c
% 10 Pi_c
% 11 phosphocreatine_c
% 12 creatine_c
% 13 AMP_c
% 14 cytocox_im
% 15 cytocred_im
% 16 Pi_im
% 17 ADP_im
% 18 ATP_im
% 19 AMP_im
% 20 h_matrix
% 21 m_matrix
% 22 k_matrix
% 23 h_c
% 24 m_c
% 25 k_c
% 26 h_im
% 29 m_im
% 28 k_im
% 29 DPsi_im_to_matrix

%% Tissue and cell composition constants
% All constants based on [https://doi.org/10.1152/ajpheart.00478.2003]
% with the mitochonrial content scaled by measured oxidative capacity
Ox_capacity = var(2);
VRegion_cytoplasm = 0.6801 + 0.2882*(1-Ox_capacity); % cytosolic volume fraction (mL cyto / mL cell) 
VRegion_im = 0.2882*Ox_capacity; % mito volume fraction (mL mito / mL cell)
VWater_cytoplasm = 0.8425; % cytosol water fraction (mL water / mL cytosol)
VWater_matrix = 0.6514; % mL matrix water / mito volume
VWater_im = 0.0724; % mL IMS water / mito volume

%% THERMODYNAMIC DATA
T = var(1);
RT = 8.314*(T+273.15)/1e3; % kJ  mol^{-1}
F = 0.096484; % kJ mol^{-1} mV^{-1}

%% STATE VARIABLES

% Concentrations of Reference Species
NAD_matrix = x(1);
NADH_matrix = x(2);
ADP_matrix = x(3);
ATP_matrix = x(4);
Pi_matrix = x(5);
coQ_matrix = x(6);
coQH2_matrix = x(7);
ATP_c = x(8);
ADP_c = x(9);
Pi_c = x(10);
phosphocreatine_c = x(11);
creatine_c = x(12);
AMP_c = x(13);
cytocox_im = x(14);
cytocred_im = x(15);
Pi_im = x(16);
ADP_im = x(17);
ATP_im = x(18);
AMP_im = x(19);
% Concentrations of H, Mg, and K
h_matrix = x(20);
m_matrix = x(21);
k_matrix = x(22);
h_c = x(23);
m_c = x(24);
k_c = x(25);
h_im = x(26);
m_im = x(27);
k_im = x(28);
% Membrane potential
DPsi_im_to_matrix = x(29);

%% DISSOCIATION CONSTANTS
% NAD_matrix
Kh(1) = Inf;
Km(1) = Inf;
Kk(1) = Inf;
% NADH_matrix
Kh(2) = Inf;
Km(2) = Inf;
Kk(2) = Inf;
% ADP_matrix
Kh(3) = 4.1057373286278706e-07;
Km(3) = 0.00071485288102034542;
Kk(3) = 0.1319048728526625;
% ATP_matrix
Kh(4) = 2.7566016511811682e-07;
Km(4) = 8.4303632255990519e-05;
Kk(4) = 0.09708512798373839;
% Pi_matrix
Kh(5) = 2.3075318323898138e-07;
Km(5) = 0.028149416570991275;
Kk(5) = 0.38034165601214515;
% coQ_matrix
Kh(6) = Inf;
Km(6) = Inf;
Kk(6) = Inf;
% coQH2_matrix
Kh(7) = Inf;
Km(7) = Inf;
Kk(7) = Inf;
% ATP_c
Kh(8) = 2.7566016511811682e-07;
Km(8) = 8.4303632255990519e-05;
Kk(8) = 0.09708512798373839;
% ADP_c
Kh(9) = 4.1057373286278706e-07;
Km(9) = 0.00071485288102034542;
Kk(9) = 0.1319048728526625;
% Pi_c
Kh(10) = 2.3075318323898138e-07;
Km(10) = 0.028149416570991275;
Kk(10) = 0.38034165601214515;
% phosphocreatine_c
Kh(11) = 4.0591821589741896e-05;
Km(11) = 0.043284385904395387;
Kk(11) = 0.58907947769630464;
% creatine_c
Kh(12) = 0.0023442288153199225;
Km(12) = Inf;
Kk(12) = Inf;
% AMP_c
Kh(13) = 6.2176978692908864e-07;
Km(13) = 0.019462203648493423;
Kk(13) = 0.23997936127262226;
% cytocox_im
Kh(14) = Inf;
Km(14) = Inf;
Kk(14) = Inf;
% cytocred_im
Kh(15) = Inf;
Km(15) = Inf;
Kk(15) = Inf;
% Pi_im
Kh(16) = 2.3075318323898138e-07;
Km(16) = 0.028149416570991275;
Kk(16) = 0.38034165601214515;
% ADP_im
Kh(17) = 4.1057373286278706e-07;
Km(17) = 0.00071485288102034542;
Kk(17) = 0.1319048728526625;
% ATP_im
Kh(18) = 2.7566016511811682e-07;
Km(18) = 8.4303632255990519e-05;
Kk(18) = 0.09708512798373839;
% AMP_im
Kh(19) = 6.2176978692908864e-07;
Km(19) = 0.019462203648493423;
Kk(19) = 0.23997936127262226;

%% BINDING POLYNOMIALS
% matrix species:
P( 1 ) = 1  + h_matrix/Kh(1) + m_matrix/Km(1) + k_matrix/Kk(1);
P( 2 ) = 1  + h_matrix/Kh(2) + m_matrix/Km(2) + k_matrix/Kk(2);
P( 3 ) = 1  + h_matrix/Kh(3) + m_matrix/Km(3) + k_matrix/Kk(3);
P( 4 ) = 1  + h_matrix/Kh(4) + m_matrix/Km(4) + k_matrix/Kk(4);
P( 5 ) = 1  + h_matrix/Kh(5) + m_matrix/Km(5) + k_matrix/Kk(5);
P( 6 ) = 1  + h_matrix/Kh(6) + m_matrix/Km(6) + k_matrix/Kk(6);
P( 7 ) = 1  + h_matrix/Kh(7) + m_matrix/Km(7) + k_matrix/Kk(7);
% cytosolic species:
P( 8 ) = 1  + h_c/Kh(8) + m_c/Km(8) + k_c/Kk(8);
P( 9 ) = 1  + h_c/Kh(9) + m_c/Km(9) + k_c/Kk(9);
P( 10 ) = 1  + h_c/Kh(10) + m_c/Km(10) + k_c/Kk(10);
P( 11 ) = 1  + h_c/Kh(11) + m_c/Km(11) + k_c/Kk(11);
P( 12 ) = 1  + h_c/Kh(12) + m_c/Km(12) + k_c/Kk(12);
P( 13 ) = 1  + h_c/Kh(13) + m_c/Km(13) + k_c/Kk(13);
% IMS species:
P( 14 ) = 1  + h_im/Kh(14) + m_im/Km(14) + k_im/Kk(14);
P( 15 ) = 1  + h_im/Kh(15) + m_im/Km(15) + k_im/Kk(15);
P( 16 ) = 1  + h_im/Kh(16) + m_im/Km(16) + k_im/Kk(16);
P( 17 ) = 1  + h_im/Kh(17) + m_im/Km(17) + k_im/Kk(17);
P( 18 ) = 1  + h_im/Kh(18) + m_im/Km(18) + k_im/Kk(18);
P( 19 ) = 1  + h_im/Kh(19) + m_im/Km(19) + k_im/Kk(19);

%% FLUX EQUATIONS (from Beard 2005)
%DH_matrix
x_DH=par(1);
r_DH=35.4;
k_DH=2.87e-2;
n_DH=1.4452;
nd=NAD_matrix;
ratio=(ADP_matrix*Pi_matrix/ATP_matrix)^n_DH;
%J_DH_matrix=x_DH*(r_DH*min(nd,2.58e-3)-NADH_matrix)/(1+(ATP_matrix*k_DH/ADP_matrix/Pi_matrix)^n_DH);
J_DH_matrix=x_DH*ratio*(r_DH*min(nd,2.58e-3)-NADH_matrix)/(ratio+k_DH^n_DH);
%SDH_BBV_matrix
x_SDH=par(2);
J_SDH_BBV_matrix=x_SDH*J_DH_matrix;
%ATPASE_cytoplasm
x_ATPase=par(3);
% Ki_ADP=2.41e-4;
dG0=4.5985;
Keq_ATPASE_cytoplasm=exp(-dG0/RT)*(P(9)*P(10)/P(8))/h_c;
dGrapp = -RT*log(Keq_ATPASE_cytoplasm);
dGrATPase = dGrapp + RT*log(ADP_c*Pi_c/ATP_c);
% Cytosolic ATP hydrolysis flux
% J_ATPASE_cytoplasm=x_ATPase/(1+ADP_c/Ki_ADP)*(1-ADP_c*Pi_c/ATP_c/Keq_ATPASE_cytoplasm);
J_ATPASE_cytoplasm = x_ATPase;
% Cytosolic creatine kinase flux
%CK_cytoplasm
x_CK=1e7;
K_CK=exp(50.78/RT);
ATP_c1=ATP_c*1/P(8);
%unboundspecies;
ADP_c1=ADP_c*1/P(9);
%unboundspecies;
J_CK_cytoplasm=x_CK*(K_CK*ADP_c1*phosphocreatine_c*h_c-ATP_c1*creatine_c);
%AK_cytoplasm
DGro_AK =-0.0036506;
Keq_AK_cytoplasm = exp(-DGro_AK/RT)*P(8)/P(9)^2*P(13);
x_AK=1e7;
% Cytosolic adenylate kinase flux
J_AK_cytoplasm=x_AK*(Keq_AK_cytoplasm*ADP_c*ADP_c-AMP_c*ATP_c);
%Modified,06/11/08
%F1F0ATPASE:im_to_matrix
k1_F1F0 =par(5);
dG0 = -4.5985;
nH = 8/3;
% Keq_F1F0ATPASE = exp(-(dG0 - nH*F*DPsi_im_to_matrix)/RT)*(h_im^3/h_matrix^2)*7.6406/2.31/1.4402;
Keq_F1F0ATPASE = exp(-(dG0 - nH*F*DPsi_im_to_matrix)/RT)*(h_im^3/h_matrix^2)*P(4)/P(3)/P(5);
J_F1F0ATPASE_im_to_matrix = k1_F1F0*(Keq_F1F0ATPASE*ADP_matrix*Pi_matrix-ATP_matrix);
%ETC1:im_to_matrix
ETC1Activity =par(6);
% Xcp0_NADH = 3.5e-3;
% Kn_NADH = 0.3e-3;
% compute free NADH ratio;
% NADH_free = (NADH_matrix - Kn_NADH - Xcp0_NADH + sqrt((NADH_matrix - Kn_NADH - Xcp0_NADH)^2 + 4*Kn_NADH *NADH_matrix))/2;
% beta_n=1;
Hp = h_im;
Hn = h_matrix;
NADH1 = NADH_matrix;
% NADH_free;
NAD1 = NAD_matrix;
Q1 = coQ_matrix;
QH21 = coQH2_matrix;
O2 = 200e-6; % constant oxygen concentration (M)
H2O2 = 0;
SO = 0;
% Update midpoint potentials from thermodynamic data
dGf_UQH2 = -19.150559;
dGf_UQ = 69.2309743;
dGf_NADH = 40.342360;
dGf_NAD = 19.230826;
dGf_SO = 12.4107669;
dGf_H2O2 = -132.491533;
dGf_O2 = 16.4;
CI_Em0_Q_QH2 = (dGf_UQH2 - dGf_UQ)/(-2*F);
% this is 18.67 mV less than value used for CI parameterization
CI_Em0_NADH = (dGf_NADH - dGf_NAD)/(-2*F);
% this is 2.3 mV higher than value used for CI parameterization
CI_Em0_SO = (dGf_SO - dGf_O2)/(-1*F);
% this is 8.3 mV higher than value used for CI parameterization
CI_Em0_H2O2 = (dGf_H2O2 - dGf_O2)/(-2*F);
% this is 8.4 mV less than value used for CI parameterization
% Update midpoint potentials for 1e- quinone couples (pH 0)
CI_Kstability = 10;
CI_Em0_Q_SQ = CI_Em0_Q_QH2 + RT/F/2*log(CI_Kstability*1e-14);
% mV
CI_Em0_SQ_QH2 = 2*CI_Em0_Q_QH2 - CI_Em0_Q_SQ;
% mV
% Binding Polynomials for Protonated States
CI_KiH1 = 4.0722e-08;
CI_KiH3 = 3.8958e-07;
PH1 = (1/(1+Hn/CI_KiH1));
PH2 = 1;
PH3 = (Hn/CI_KiH3/(1+Hn/CI_KiH3));
% Binding Polynomials for enzyme, substrates, products and regulators
% NADH binding constants
KdNADHo=4.6065e-05;
KdNADo=7.0545e-04;
KdNADHr=4.9867e-04;
KdNADr=1.1844e-05;
KdNADHrad=1;
KdNADrad=1.5358e-04;
% Quinone binding constants
KdQH2o = 0.1;
KdQo = 0.0175;
KdQH2r = 0.1;
KdQr = 0.0175;
% Binding polynomials
PNo = 1 + NADH1/KdNADHo + NAD1/KdNADo;
% oxidized binding constants
PNr = 1 + NADH1/KdNADHr + NAD1/KdNADr;
% reduced binding constants
PNrad = 1 + NADH1/KdNADHrad + NAD1/KdNADrad;
% reduced binding constants
PQr = 1 + QH21/KdQH2r + Q1/KdQr;
% reduced binding constants
PQo = 1 + QH21/KdQH2o + Q1/KdQo;
% oxidized binding constants
muH = F*DPsi_im_to_matrix+RT*log(Hp/Hn);
% proton chemical potential (kJ/mol)
% Compute pH Corrected Midpoint potentials
% NADH potential
% CI_Em0_FMN = 55.1446;
Em_NADH = CI_Em0_NADH + log(10)*RT/F/2*log10(Hn);
% FMN potentials
CI_Em0_FMN = 55.1446;
% CI_FMNred.pK = 7.0998;
CI_Em0_FMN2 = 86.7495;
CI_Em0_FMN1 = 23.5369;
CI_FMNrad = 7.9074;
CI_FMNred_pK = 7.0998;
Em_FMNred = CI_Em0_FMN - RT/2/F*log(Hn/10^-CI_FMNred_pK/(1+Hn/10^-CI_FMNred_pK)/Hn.^2);
% FMNred/FMN, pK1 of fully reduced FMN very high (>10)0)
Em_FMNrad = CI_Em0_FMN2 - RT/F*log((Hn/10^-(CI_FMNred_pK)/(1+Hn/10^-CI_FMNred_pK))/((Hn/10^-CI_FMNrad)/(1+Hn/10^-CI_FMNrad))/Hn);
% FMNred/FMNrad
Em_FMNox = CI_Em0_FMN1 - RT/F*log((Hn/10^-CI_FMNrad)/(1+Hn/10^-CI_FMNrad)/Hn);
%FMN/FMNrad
% N2 potential
CI_Em0_N2 = -90;
CI_N2_pKox = 6;
CI_N2_pKred = 8.5;
Em_N2 = CI_Em0_N2 - RT/F*log((Hn+10^-CI_N2_pKox)/(Hn+10^-CI_N2_pKred));
% superoxide and H2O2 potentials
Em_SO = CI_Em0_SO;
% Em0 stays the same for pH > 6, pKa of SO is ~4.8
Em_H2O2 = CI_Em0_H2O2 + log(10)*RT/F*log10(Hn);
% 1st pKa is ~11
% Quinone potentials
% Em_Q_QH2 = CI_Em0_Q_QH2 + 2*log(10)*RT/F/2*log10(Hn);
% assuming linked to N-side
Em_Q_SQ = CI_Em0_Q_SQ;
% pH independent
Em_SQ_QH2 = CI_Em0_SQ_QH2 + 2*log(10)*RT/F*log10(Hn);
% assuming linked to N-side
% State Transition Thermodynamics
% NADH-QH2 reductase related
K_NADHFMNred = exp((2*F*(Em_FMNred - Em_NADH))/RT);
% [Fr][NAD]/[F][NADH][H]
K_FMNradN2 = exp((F*(Em_N2 - Em_FMNrad))/RT);
% [Frad][N2r]/[Fr]/[N2]
K_FMNoxN2 = exp((F*(Em_N2 - Em_FMNox))/RT);
% [F][N2r]/[Frad]/[N2]
K_N2SQ = exp((F*(Em_Q_SQ - Em_N2))/RT);
% [N2][SQ]/[N2r][Q]
K_N2QH2 = exp((F*(Em_SQ_QH2 - Em_N2))/RT);
% [N2][QH2]/[N2r][SQ]
% superoxide and hydrogen peroxide
K_FMNredH2O2 = exp((2*F*(Em_H2O2 - Em_FMNred))/RT);
% [F][H2O2]/[Fr][O2]
K_FMNradO2 = exp((F*(Em_SO - Em_FMNrad))/RT);
% [Frad][O2.-]/[Fr][O2]
K_FMNoxO2 = exp((F*(Em_SO - Em_FMNox))/RT);
% [F][O2.-]/[Frad][O2]
K_N2O2 = exp((F*(Em_SO - Em_N2))/RT);
% [N2][O2.-]/[N2r][O2]
K_SQO2 = exp((F*(Em_SO - Em_Q_SQ))/RT);
% [Q][O2.-]/[SQ][O2]
% Substates via Rapid Equilibrium
% substate equilibrium constants
K1 = K_FMNradN2;
%[Frad][N2r]/[Fr][N2]
K2 = K_FMNoxN2;
% [F][N2r]/[Frad][N2]
KSQ = K_N2SQ*KdQr;
% [N2-SQ]/[N2r-Q], KdQr = [N2r][Q]/[N2r-Q]
KQ = (Q1/KdQr/PQr)*KSQ;
% [N2-SQ]/[N2r]
% 0e-
PS0 = 1;
% total
s0_F_N2 = 1/PS0;
% 1e-
PS1 = (1 + K2*(1 + KQ));
% total sum(cell2mat(struct2cell(s1)))
s1_Frad_N2 = 1/PS1;
s1_F_N2r = K2*s1_Frad_N2;
s1_F_N2_SQ = KQ*s1_F_N2r;
% 2e-
PS2 = (1 + K1*(1 + KQ*(1 + K2)));
% total sum(cell2mat(struct2cell(s2)))
s2_Fr_N2 = 1/PS2;
s2_Frad_N2r = K1*s2_Fr_N2;
s2_Frad_N2_SQ = KQ*s2_Frad_N2r;
s2_F_N2r_SQ = K2*s2_Frad_N2_SQ;
% 3e-
PS3 = (1 + KQ*(1 + K1));
% total sum(cell2mat(struct2cell(s3)))
s3_Fr_N2r = 1/PS3;
s3_Fr_N2_SQ = KQ*s3_Fr_N2r;
s3_Frad_N2r_SQ = K1*s3_Fr_N2_SQ;
% 4e-
PS4 = 1;
% total
s4_Fr_N2r_SQ = 1/PS4;

% State Transition Rates
% reverse NADH-Q reductase rates
CI_kfNADH_02 = 1.9642e+03;
CI_kfNADH_24 = 186.3167;
CI_kfNADH_13 = 4.6093;
CI_kfQ_20 = 5.8175e+03;
CI_kfQ_42 = 7.7767e+10;
CI_kfQ_31 = 8.6596;
krN0_02 = CI_kfNADH_02/K_NADHFMNred*KdNADr/KdNADHo*PH1*PH2;
krN0_24 = CI_kfNADH_24/K_NADHFMNred*KdNADr/KdNADHo*PH1*PH2;
krN0_13 = CI_kfNADH_13/K_NADHFMNred*KdNADr/KdNADHo*PH1*PH2;
krQ0_20 = CI_kfQ_20/K_N2QH2*KdQH2o*(PNr*PQo/PNo/PQr)*PH3;
krQ0_42 = CI_kfQ_42/K_N2QH2*KdQH2o*(PNr*PQo/PNo/PQr)*PH3;
krQ0_31 = CI_kfQ_31/K_N2QH2*KdQH2o*(PNr*PQo/PNo/PQr)*PH3;
% reverse ROS rates
% O2 and SQ
CI_kfSQSO_10 = 1.6667e+09;
CI_kfSQSO_21 = 0.0021;
CI_kfSQSO_21b = 0.0021;
CI_kfSQSO_32 = 1.1450e-05;
CI_kfSQSO_32b = 2.3953e-09;
CI_kfSQSO_43 = 2.1602e-07;
krSQSO_10 = CI_kfSQSO_10/K_SQO2;
krSQSO_21 = CI_kfSQSO_21/K_SQO2;
krSQSO_21b = CI_kfSQSO_21b/K_SQO2;
krSQSO_32 = CI_kfSQSO_32/K_SQO2;
krSQSO_32b = CI_kfSQSO_32b/K_SQO2;
krSQSO_43 = CI_kfSQSO_43/K_SQO2;
% O2 and N2
CI_kfN2rSO_10 = 0;
CI_kfN2rSO_21 = 0;
CI_kfN2rSO_21b = 0;
CI_kfN2rSO_32 = 0;
CI_kfN2rSO_32b = 0;
CI_kfN2rSO_43 = 0;
krN2rSO_10 = CI_kfN2rSO_10/K_N2O2;
krN2rSO_21 = CI_kfN2rSO_21/K_N2O2;
krN2rSO_21b = CI_kfN2rSO_21b/K_N2O2;
krN2rSO_32 = CI_kfN2rSO_32/K_N2O2;
krN2rSO_32b = CI_kfN2rSO_32b/K_N2O2;
krN2rSO_43 = CI_kfN2rSO_43/K_N2O2;
% O2 and FMNred
CI_kfFrSO_21 = 451900;
CI_kfFradSO_21b = 1.1943e-07;
CI_kfFrSO_32 = 5.3833e-04;
CI_kfFrSO_32b = 2.8333e-05;
CI_kfFrSO_43 = 1.4005e+05;
krFrSO_21 = CI_kfFrSO_21/K_FMNradO2;
krFradSO_21b = CI_kfFradSO_21b/K_FMNradO2;
krFrSO_32 = CI_kfFrSO_32/K_FMNradO2;
krFrSO_32b = CI_kfFrSO_32b/K_FMNradO2;
krFrSO_43 = CI_kfFrSO_43/K_FMNradO2;
% O2 and FMNrad
CI_kfFradSO_10 = 68055000;
CI_kfFradSO_21 = 4.0417e-11;
CI_kfFradSO_32 = 0.014;
krFradSO_10 = CI_kfFradSO_10/K_FMNoxO2;
krFradSO_21 = CI_kfFradSO_21/K_FMNoxO2;
krFradSO_32 = CI_kfFradSO_32/K_FMNoxO2;
% H2O2 and FMNred
CI_kfH_20 = 7.8243e-07;
CI_kfH_31 = 2.1228e+06;
CI_kfH_31b = 1.7637e+06;
CI_kfH_42 = 62.0583;
krH_20 = CI_kfH_20/K_FMNredH2O2;
krH_31 = CI_kfH_31/K_FMNredH2O2;
krH_31b = CI_kfH_31b/K_FMNredH2O2;
krH_42 = CI_kfH_42/K_FMNredH2O2;
% Nucleotide and pH effects
CI_kfH_20 = CI_kfH_20/PNr;
CI_kfH_42 = CI_kfH_42/PNr;
CI_kfH_31 = CI_kfH_31/PNr;
CI_kfH_31b = CI_kfH_31b/PNr;
CI_kfFrSO_21 = CI_kfFrSO_21/PNr/(1+Hn/10^-CI_FMNred_pK);
CI_kfFrSO_32 = CI_kfFrSO_32/PNr/(1+Hn/10^-CI_FMNred_pK);
CI_kfFrSO_32b = CI_kfFrSO_32b/PNr/(1+Hn/10^-CI_FMNred_pK);
CI_kfFrSO_43 = CI_kfFrSO_43/PNr/(1+Hn/10^-CI_FMNred_pK);
CI_kfFradSO_10 = CI_kfFradSO_10/PNrad/(1+Hn/10^-CI_FMNrad);
CI_kfFradSO_21 = CI_kfFradSO_21/PNrad/(1+Hn/10^-CI_FMNrad);
CI_kfFradSO_21b = CI_kfFradSO_21b/PNrad/(1+Hn/10^-CI_FMNrad);
CI_kfFradSO_32 = CI_kfFradSO_32/PNrad/(1+Hn/10^-CI_FMNrad);
CI_kfN2rSO_10 = CI_kfN2rSO_10/(1+Hn/10^-CI_N2_pKred);
CI_kfN2rSO_21 = CI_kfN2rSO_21/(1+Hn/10^-CI_N2_pKred);
CI_kfN2rSO_21b = CI_kfN2rSO_21b/(1+Hn/10^-CI_N2_pKred);
CI_kfN2rSO_32 = CI_kfN2rSO_32/(1+Hn/10^-CI_N2_pKred);
CI_kfN2rSO_32b = CI_kfN2rSO_32b/(1+Hn/10^-CI_N2_pKred);
CI_kfN2rSO_43 = CI_kfN2rSO_43/(1+Hn/10^-CI_N2_pKred);
% net forward rates % membrane potential dependence
dPsiNf = 1;
dPsiNr = 1;
CI_beta1 = 0.5;
dPsiQf = exp(-4*CI_beta1*muH/RT);
dPsiQr = exp(4*(1-CI_beta1)*muH/RT);
% NADH transitions
kfNADH_02 = CI_kfNADH_02*dPsiNf*NADH1/KdNADHo/PNo*(s0_F_N2)*PH1*PH2;
krNADH_02 = krN0_02*dPsiNr*NAD1/KdNADr/PNr*(s2_Fr_N2);
kfNADH_24 = CI_kfNADH_24*dPsiNf*NADH1/KdNADHo/PNo*(s2_F_N2r_SQ)*PH1*PH2;
krNADH_24 = krN0_24*dPsiNr*NAD1/KdNADr/PNr*(s4_Fr_N2r_SQ);
kfNADH_13 = CI_kfNADH_13*dPsiNf*NADH1/KdNADHo/PNo*(s1_F_N2r + s1_F_N2_SQ)*PH1*PH2;
krNADH_13 = krN0_13*dPsiNr*NAD1/KdNADr/PNr*(s3_Fr_N2r + s3_Fr_N2_SQ);
% QH2 transitions
kfQ_20 = CI_kfQ_20*dPsiQf*(s2_F_N2r_SQ)*PH3;
krQ_20 = krQ0_20*dPsiQr*QH21/KdQH2o/PQo*(s0_F_N2);
kfQ_42 = CI_kfQ_42*dPsiQf*(s4_Fr_N2r_SQ)*PH3;
krQ_42 = krQ0_42*dPsiQr*QH21/KdQH2o/PQo*(s2_Fr_N2);
kfQ_31 = CI_kfQ_31*dPsiQf*(s3_Frad_N2r_SQ)*PH3;
krQ_31 = krQ0_31*dPsiQr*QH21/KdQH2o/PQo*(s1_Frad_N2);
% H2O2 transitions
kfH2O2_20 = CI_kfH_20*O2*s2_Fr_N2;
krH2O2_20 = krH_20*H2O2*s0_F_N2;
kfH2O2_31 = CI_kfH_31*O2*s3_Fr_N2r + CI_kfH_31b*O2*s3_Fr_N2_SQ;
krH2O2_31 = krH_31*H2O2*s1_F_N2r + krH_31b*H2O2*s1_F_N2_SQ;
kfH2O2_42 = CI_kfH_42*O2*s4_Fr_N2r_SQ;
krH2O2_42 = krH_42*H2O2*(s2_F_N2r_SQ);
% SO transitions
kfSO_10 = O2*(CI_kfSQSO_10*s1_F_N2_SQ + CI_kfFradSO_10*s1_Frad_N2 + CI_kfN2rSO_10*s1_F_N2r);
krSO_10 = SO*s0_F_N2*(krSQSO_10*(Q1/KdQo/PQo) + krFradSO_10 + krN2rSO_10);
kfSO_21 = O2*((CI_kfSQSO_21*s2_Frad_N2_SQ + CI_kfSQSO_21b*s2_F_N2r_SQ) + CI_kfFrSO_21*s2_Fr_N2 + (CI_kfFradSO_21*s2_Frad_N2r + CI_kfFradSO_21b*s2_Frad_N2_SQ) + (CI_kfN2rSO_21*s2_Frad_N2r + CI_kfN2rSO_21b*s2_F_N2r_SQ));
krSO_21 = SO*(krFrSO_21*s1_Frad_N2 + krN2rSO_21*s1_Frad_N2 + krN2rSO_21b*s1_F_N2_SQ + krFradSO_21*s1_F_N2r + krFradSO_21b*s1_F_N2_SQ + (krSQSO_21*s1_Frad_N2*(Q1/KdQo/PQo) + krSQSO_21b*s1_F_N2r*(Q1/KdQr/PQr)));
kfSO_32 = O2*((CI_kfSQSO_32*s3_Fr_N2_SQ + CI_kfSQSO_32b*s3_Frad_N2r_SQ) + (CI_kfFrSO_32*s3_Fr_N2r + CI_kfFrSO_32b*s3_Fr_N2_SQ) + CI_kfFradSO_32*s3_Frad_N2r_SQ + (CI_kfN2rSO_32*s3_Fr_N2r + CI_kfN2rSO_32b*s3_Frad_N2r_SQ));
krSO_32 = SO*(krSQSO_32*s2_Fr_N2*(Q1/KdQo/PQo) + krSQSO_32b*s2_Frad_N2r*(Q1/KdQr/PQr) + krFradSO_32*s2_F_N2r_SQ + krFrSO_32*s2_Frad_N2r + krFrSO_32b*s2_Frad_N2_SQ + krN2rSO_32*s2_Fr_N2 + krN2rSO_32b*s2_Frad_N2_SQ);
kfSO_43 = O2*((CI_kfSQSO_43 + CI_kfN2rSO_43)*(s4_Fr_N2r_SQ) + CI_kfFrSO_43*s4_Fr_N2r_SQ);
krSO_43 = SO*(krSQSO_43*s3_Fr_N2r*(Q1/KdQr/PQr) + krFrSO_43*s3_Frad_N2r_SQ + krN2rSO_43*s3_Fr_N2_SQ);
% S0 <-> S2
k0_2 = kfNADH_02 + krQ_20 + krH2O2_20;
k2_0 = kfQ_20 + krNADH_02 + kfH2O2_20;
% S2 <-> S4
k2_4 = kfNADH_24 + krQ_42 + krH2O2_42;
k4_2 = kfQ_42 + krNADH_24 + kfH2O2_42;
% S1 <-> S3
k1_3 = kfNADH_13 + krQ_31 + krH2O2_31;
k3_1 = kfQ_31 + krNADH_13 + kfH2O2_31;
% S0 <-> S1 (superoxide)
k1_0 = kfSO_10;
k0_1 = krSO_10;
% S2 <-> S1 (superoxide)
k2_1 = kfSO_21;
k1_2 = krSO_21;
% S3 <-> S2 (superoxide)
k3_2 = kfSO_32;
k2_3 = krSO_32;
% S4 <-> S3 (superoxide)
k4_3 = kfSO_43;
k3_4 = krSO_43;
% Steady-State Fractional Occupancies (solved analytically)
S0=(k1_0*k2_0*k3_1*k4_2 + k1_0*k2_0*k3_1*k4_3 + k1_0*k2_0*k3_2*k4_2 + k1_0*k2_1*k3_1*k4_2 + k1_0*k2_0*k3_2*k4_3 + k1_0*k2_1*k3_1*k4_3 + k1_0*k2_1*k3_2*k4_2 + k1_2*k2_0*k3_1*k4_2 + k1_0*k2_0*k3_4*k4_2 + k1_0*k2_1*k3_2*k4_3 + k1_0*k2_3*k3_1*k4_2 + k1_2*k2_0*k3_1*k4_3 + k1_2*k2_0*k3_2*k4_2 + k1_0*k2_1*k3_4*k4_2 + k1_0*k2_3*k3_1*k4_3 + k1_2*k2_0*k3_2*k4_3 + k1_3*k2_0*k3_2*k4_2 + k1_0*k2_4*k3_1*k4_3 + k1_2*k2_0*k3_4*k4_2 + k1_3*k2_0*k3_2*k4_3 + k1_3*k2_0*k3_4*k4_2)/(k0_2*k1_0*k2_4*k3_1 + k0_1*k1_2*k2_4*k3_1 + k0_1*k1_3*k2_0*k3_4 + k0_2*k1_0*k2_4*k3_2 + k0_1*k1_2*k2_4*k3_2 + k0_1*k1_3*k2_1*k3_4 + k0_2*k1_0*k2_3*k3_4 + k0_2*k1_2*k2_4*k3_1 + k0_1*k1_2*k2_3*k3_4 + k0_1*k1_3*k2_4*k3_2 + k0_2*k1_0*k2_4*k3_4 + k0_2*k1_2*k2_4*k3_2 + k0_2*k1_3*k2_1*k3_4 + k0_1*k1_2*k2_4*k3_4 + k0_1*k1_3*k2_3*k3_4 + k0_2*k1_2*k2_3*k3_4 + k0_2*k1_3*k2_4*k3_2 + k0_1*k1_3*k2_4*k3_4 + k0_2*k1_2*k2_4*k3_4 + k0_2*k1_3*k2_3*k3_4 + k0_2*k1_3*k2_4*k3_4 + k0_1*k1_3*k2_0*k4_2 + k0_1*k1_3*k2_0*k4_3 + k0_1*k1_3*k2_1*k4_2 + k0_2*k1_0*k2_3*k4_2 + k0_1*k1_2*k2_3*k4_2 + k0_1*k1_3*k2_1*k4_3 + k0_2*k1_0*k2_3*k4_3 + k0_2*k1_3*k2_1*k4_2 + k0_1*k1_2*k2_3*k4_3 + k0_1*k1_3*k2_3*k4_2 + k0_2*k1_0*k2_4*k4_3 + k0_2*k1_2*k2_3*k4_2 + k0_2*k1_3*k2_1*k4_3 + k0_1*k1_2*k2_4*k4_3 + k0_1*k1_3*k2_3*k4_3 + k0_2*k1_2*k2_3*k4_3 + k0_2*k1_3*k2_3*k4_2 + k0_1*k1_3*k2_4*k4_3 + k0_2*k1_2*k2_4*k4_3 + k0_2*k1_3*k2_3*k4_3 + k0_2*k1_3*k2_4*k4_3 + k0_2*k1_0*k3_1*k4_2 + k0_1*k1_2*k3_1*k4_2 + k0_2*k1_0*k3_1*k4_3 + k0_2*k1_0*k3_2*k4_2 + k0_1*k1_2*k3_1*k4_3 + k0_1*k1_2*k3_2*k4_2 + k0_2*k1_0*k3_2*k4_3 + k0_2*k1_2*k3_1*k4_2 + k0_1*k1_2*k3_2*k4_3 + k0_1*k1_3*k3_2*k4_2 + k0_2*k1_0*k3_4*k4_2 + k0_2*k1_2*k3_1*k4_3 + k0_2*k1_2*k3_2*k4_2 + k0_1*k1_2*k3_4*k4_2 + k0_1*k1_3*k3_2*k4_3 + k0_2*k1_2*k3_2*k4_3 + k0_2*k1_3*k3_2*k4_2 + k0_1*k1_3*k3_4*k4_2 + k0_2*k1_2*k3_4*k4_2 + k0_2*k1_3*k3_2*k4_3 + k0_2*k1_3*k3_4*k4_2 + k0_1*k2_0*k3_1*k4_2 + k0_1*k2_0*k3_1*k4_3 + k0_1*k2_0*k3_2*k4_2 + k0_1*k2_1*k3_1*k4_2 + k0_1*k2_0*k3_2*k4_3 + k0_1*k2_1*k3_1*k4_3 + k0_1*k2_1*k3_2*k4_2 + k0_2*k2_1*k3_1*k4_2 + k0_1*k2_0*k3_4*k4_2 + k0_1*k2_1*k3_2*k4_3 + k0_1*k2_3*k3_1*k4_2 + k0_2*k2_1*k3_1*k4_3 + k0_2*k2_1*k3_2*k4_2 + k0_1*k2_1*k3_4*k4_2 + k0_1*k2_3*k3_1*k4_3 + k0_2*k2_1*k3_2*k4_3 + k0_2*k2_3*k3_1*k4_2 + k0_1*k2_4*k3_1*k4_3 + k0_2*k2_1*k3_4*k4_2 + k0_2*k2_3*k3_1*k4_3 + k0_2*k2_4*k3_1*k4_3 + k1_0*k2_0*k3_1*k4_2 + k1_0*k2_0*k3_1*k4_3 + k1_0*k2_0*k3_2*k4_2 + k1_0*k2_1*k3_1*k4_2 + k1_0*k2_0*k3_2*k4_3 + k1_0*k2_1*k3_1*k4_3 + k1_0*k2_1*k3_2*k4_2 + k1_2*k2_0*k3_1*k4_2 + k1_0*k2_0*k3_4*k4_2 + k1_0*k2_1*k3_2*k4_3 + k1_0*k2_3*k3_1*k4_2 + k1_2*k2_0*k3_1*k4_3 + k1_2*k2_0*k3_2*k4_2 + k1_0*k2_1*k3_4*k4_2 + k1_0*k2_3*k3_1*k4_3 + k1_2*k2_0*k3_2*k4_3 + k1_3*k2_0*k3_2*k4_2 + k1_0*k2_4*k3_1*k4_3 + k1_2*k2_0*k3_4*k4_2 + k1_3*k2_0*k3_2*k4_3 + k1_3*k2_0*k3_4*k4_2);
S1=(k0_1*k2_0*k3_1*k4_2 + k0_1*k2_0*k3_1*k4_3 + k0_1*k2_0*k3_2*k4_2 + k0_1*k2_1*k3_1*k4_2 + k0_1*k2_0*k3_2*k4_3 + k0_1*k2_1*k3_1*k4_3 + k0_1*k2_1*k3_2*k4_2 + k0_2*k2_1*k3_1*k4_2 + k0_1*k2_0*k3_4*k4_2 + k0_1*k2_1*k3_2*k4_3 + k0_1*k2_3*k3_1*k4_2 + k0_2*k2_1*k3_1*k4_3 + k0_2*k2_1*k3_2*k4_2 + k0_1*k2_1*k3_4*k4_2 + k0_1*k2_3*k3_1*k4_3 + k0_2*k2_1*k3_2*k4_3 + k0_2*k2_3*k3_1*k4_2 + k0_1*k2_4*k3_1*k4_3 + k0_2*k2_1*k3_4*k4_2 + k0_2*k2_3*k3_1*k4_3 + k0_2*k2_4*k3_1*k4_3)/(k0_2*k1_0*k2_4*k3_1 + k0_1*k1_2*k2_4*k3_1 + k0_1*k1_3*k2_0*k3_4 + k0_2*k1_0*k2_4*k3_2 + k0_1*k1_2*k2_4*k3_2 + k0_1*k1_3*k2_1*k3_4 + k0_2*k1_0*k2_3*k3_4 + k0_2*k1_2*k2_4*k3_1 + k0_1*k1_2*k2_3*k3_4 + k0_1*k1_3*k2_4*k3_2 + k0_2*k1_0*k2_4*k3_4 + k0_2*k1_2*k2_4*k3_2 + k0_2*k1_3*k2_1*k3_4 + k0_1*k1_2*k2_4*k3_4 + k0_1*k1_3*k2_3*k3_4 + k0_2*k1_2*k2_3*k3_4 + k0_2*k1_3*k2_4*k3_2 + k0_1*k1_3*k2_4*k3_4 + k0_2*k1_2*k2_4*k3_4 + k0_2*k1_3*k2_3*k3_4 + k0_2*k1_3*k2_4*k3_4 + k0_1*k1_3*k2_0*k4_2 + k0_1*k1_3*k2_0*k4_3 + k0_1*k1_3*k2_1*k4_2 + k0_2*k1_0*k2_3*k4_2 + k0_1*k1_2*k2_3*k4_2 + k0_1*k1_3*k2_1*k4_3 + k0_2*k1_0*k2_3*k4_3 + k0_2*k1_3*k2_1*k4_2 + k0_1*k1_2*k2_3*k4_3 + k0_1*k1_3*k2_3*k4_2 + k0_2*k1_0*k2_4*k4_3 + k0_2*k1_2*k2_3*k4_2 + k0_2*k1_3*k2_1*k4_3 + k0_1*k1_2*k2_4*k4_3 + k0_1*k1_3*k2_3*k4_3 + k0_2*k1_2*k2_3*k4_3 + k0_2*k1_3*k2_3*k4_2 + k0_1*k1_3*k2_4*k4_3 + k0_2*k1_2*k2_4*k4_3 + k0_2*k1_3*k2_3*k4_3 + k0_2*k1_3*k2_4*k4_3 + k0_2*k1_0*k3_1*k4_2 + k0_1*k1_2*k3_1*k4_2 + k0_2*k1_0*k3_1*k4_3 + k0_2*k1_0*k3_2*k4_2 + k0_1*k1_2*k3_1*k4_3 + k0_1*k1_2*k3_2*k4_2 + k0_2*k1_0*k3_2*k4_3 + k0_2*k1_2*k3_1*k4_2 + k0_1*k1_2*k3_2*k4_3 + k0_1*k1_3*k3_2*k4_2 + k0_2*k1_0*k3_4*k4_2 + k0_2*k1_2*k3_1*k4_3 + k0_2*k1_2*k3_2*k4_2 + k0_1*k1_2*k3_4*k4_2 + k0_1*k1_3*k3_2*k4_3 + k0_2*k1_2*k3_2*k4_3 + k0_2*k1_3*k3_2*k4_2 + k0_1*k1_3*k3_4*k4_2 + k0_2*k1_2*k3_4*k4_2 + k0_2*k1_3*k3_2*k4_3 + k0_2*k1_3*k3_4*k4_2 + k0_1*k2_0*k3_1*k4_2 + k0_1*k2_0*k3_1*k4_3 + k0_1*k2_0*k3_2*k4_2 + k0_1*k2_1*k3_1*k4_2 + k0_1*k2_0*k3_2*k4_3 + k0_1*k2_1*k3_1*k4_3 + k0_1*k2_1*k3_2*k4_2 + k0_2*k2_1*k3_1*k4_2 + k0_1*k2_0*k3_4*k4_2 + k0_1*k2_1*k3_2*k4_3 + k0_1*k2_3*k3_1*k4_2 + k0_2*k2_1*k3_1*k4_3 + k0_2*k2_1*k3_2*k4_2 + k0_1*k2_1*k3_4*k4_2 + k0_1*k2_3*k3_1*k4_3 + k0_2*k2_1*k3_2*k4_3 + k0_2*k2_3*k3_1*k4_2 + k0_1*k2_4*k3_1*k4_3 + k0_2*k2_1*k3_4*k4_2 + k0_2*k2_3*k3_1*k4_3 + k0_2*k2_4*k3_1*k4_3 + k1_0*k2_0*k3_1*k4_2 + k1_0*k2_0*k3_1*k4_3 + k1_0*k2_0*k3_2*k4_2 + k1_0*k2_1*k3_1*k4_2 + k1_0*k2_0*k3_2*k4_3 + k1_0*k2_1*k3_1*k4_3 + k1_0*k2_1*k3_2*k4_2 + k1_2*k2_0*k3_1*k4_2 + k1_0*k2_0*k3_4*k4_2 + k1_0*k2_1*k3_2*k4_3 + k1_0*k2_3*k3_1*k4_2 + k1_2*k2_0*k3_1*k4_3 + k1_2*k2_0*k3_2*k4_2 + k1_0*k2_1*k3_4*k4_2 + k1_0*k2_3*k3_1*k4_3 + k1_2*k2_0*k3_2*k4_3 + k1_3*k2_0*k3_2*k4_2 + k1_0*k2_4*k3_1*k4_3 + k1_2*k2_0*k3_4*k4_2 + k1_3*k2_0*k3_2*k4_3 + k1_3*k2_0*k3_4*k4_2);
S2=(k0_2*k1_0*k3_1*k4_2 + k0_1*k1_2*k3_1*k4_2 + k0_2*k1_0*k3_1*k4_3 + k0_2*k1_0*k3_2*k4_2 + k0_1*k1_2*k3_1*k4_3 + k0_1*k1_2*k3_2*k4_2 + k0_2*k1_0*k3_2*k4_3 + k0_2*k1_2*k3_1*k4_2 + k0_1*k1_2*k3_2*k4_3 + k0_1*k1_3*k3_2*k4_2 + k0_2*k1_0*k3_4*k4_2 + k0_2*k1_2*k3_1*k4_3 + k0_2*k1_2*k3_2*k4_2 + k0_1*k1_2*k3_4*k4_2 + k0_1*k1_3*k3_2*k4_3 + k0_2*k1_2*k3_2*k4_3 + k0_2*k1_3*k3_2*k4_2 + k0_1*k1_3*k3_4*k4_2 + k0_2*k1_2*k3_4*k4_2 + k0_2*k1_3*k3_2*k4_3 + k0_2*k1_3*k3_4*k4_2)/(k0_2*k1_0*k2_4*k3_1 + k0_1*k1_2*k2_4*k3_1 + k0_1*k1_3*k2_0*k3_4 + k0_2*k1_0*k2_4*k3_2 + k0_1*k1_2*k2_4*k3_2 + k0_1*k1_3*k2_1*k3_4 + k0_2*k1_0*k2_3*k3_4 + k0_2*k1_2*k2_4*k3_1 + k0_1*k1_2*k2_3*k3_4 + k0_1*k1_3*k2_4*k3_2 + k0_2*k1_0*k2_4*k3_4 + k0_2*k1_2*k2_4*k3_2 + k0_2*k1_3*k2_1*k3_4 + k0_1*k1_2*k2_4*k3_4 + k0_1*k1_3*k2_3*k3_4 + k0_2*k1_2*k2_3*k3_4 + k0_2*k1_3*k2_4*k3_2 + k0_1*k1_3*k2_4*k3_4 + k0_2*k1_2*k2_4*k3_4 + k0_2*k1_3*k2_3*k3_4 + k0_2*k1_3*k2_4*k3_4 + k0_1*k1_3*k2_0*k4_2 + k0_1*k1_3*k2_0*k4_3 + k0_1*k1_3*k2_1*k4_2 + k0_2*k1_0*k2_3*k4_2 + k0_1*k1_2*k2_3*k4_2 + k0_1*k1_3*k2_1*k4_3 + k0_2*k1_0*k2_3*k4_3 + k0_2*k1_3*k2_1*k4_2 + k0_1*k1_2*k2_3*k4_3 + k0_1*k1_3*k2_3*k4_2 + k0_2*k1_0*k2_4*k4_3 + k0_2*k1_2*k2_3*k4_2 + k0_2*k1_3*k2_1*k4_3 + k0_1*k1_2*k2_4*k4_3 + k0_1*k1_3*k2_3*k4_3 + k0_2*k1_2*k2_3*k4_3 + k0_2*k1_3*k2_3*k4_2 + k0_1*k1_3*k2_4*k4_3 + k0_2*k1_2*k2_4*k4_3 + k0_2*k1_3*k2_3*k4_3 + k0_2*k1_3*k2_4*k4_3 + k0_2*k1_0*k3_1*k4_2 + k0_1*k1_2*k3_1*k4_2 + k0_2*k1_0*k3_1*k4_3 + k0_2*k1_0*k3_2*k4_2 + k0_1*k1_2*k3_1*k4_3 + k0_1*k1_2*k3_2*k4_2 + k0_2*k1_0*k3_2*k4_3 + k0_2*k1_2*k3_1*k4_2 + k0_1*k1_2*k3_2*k4_3 + k0_1*k1_3*k3_2*k4_2 + k0_2*k1_0*k3_4*k4_2 + k0_2*k1_2*k3_1*k4_3 + k0_2*k1_2*k3_2*k4_2 + k0_1*k1_2*k3_4*k4_2 + k0_1*k1_3*k3_2*k4_3 + k0_2*k1_2*k3_2*k4_3 + k0_2*k1_3*k3_2*k4_2 + k0_1*k1_3*k3_4*k4_2 + k0_2*k1_2*k3_4*k4_2 + k0_2*k1_3*k3_2*k4_3 + k0_2*k1_3*k3_4*k4_2 + k0_1*k2_0*k3_1*k4_2 + k0_1*k2_0*k3_1*k4_3 + k0_1*k2_0*k3_2*k4_2 + k0_1*k2_1*k3_1*k4_2 + k0_1*k2_0*k3_2*k4_3 + k0_1*k2_1*k3_1*k4_3 + k0_1*k2_1*k3_2*k4_2 + k0_2*k2_1*k3_1*k4_2 + k0_1*k2_0*k3_4*k4_2 + k0_1*k2_1*k3_2*k4_3 + k0_1*k2_3*k3_1*k4_2 + k0_2*k2_1*k3_1*k4_3 + k0_2*k2_1*k3_2*k4_2 + k0_1*k2_1*k3_4*k4_2 + k0_1*k2_3*k3_1*k4_3 + k0_2*k2_1*k3_2*k4_3 + k0_2*k2_3*k3_1*k4_2 + k0_1*k2_4*k3_1*k4_3 + k0_2*k2_1*k3_4*k4_2 + k0_2*k2_3*k3_1*k4_3 + k0_2*k2_4*k3_1*k4_3 + k1_0*k2_0*k3_1*k4_2 + k1_0*k2_0*k3_1*k4_3 + k1_0*k2_0*k3_2*k4_2 + k1_0*k2_1*k3_1*k4_2 + k1_0*k2_0*k3_2*k4_3 + k1_0*k2_1*k3_1*k4_3 + k1_0*k2_1*k3_2*k4_2 + k1_2*k2_0*k3_1*k4_2 + k1_0*k2_0*k3_4*k4_2 + k1_0*k2_1*k3_2*k4_3 + k1_0*k2_3*k3_1*k4_2 + k1_2*k2_0*k3_1*k4_3 + k1_2*k2_0*k3_2*k4_2 + k1_0*k2_1*k3_4*k4_2 + k1_0*k2_3*k3_1*k4_3 + k1_2*k2_0*k3_2*k4_3 + k1_3*k2_0*k3_2*k4_2 + k1_0*k2_4*k3_1*k4_3 + k1_2*k2_0*k3_4*k4_2 + k1_3*k2_0*k3_2*k4_3 + k1_3*k2_0*k3_4*k4_2);
S3=(k0_1*k1_3*k2_0*k4_2 + k0_1*k1_3*k2_0*k4_3 + k0_1*k1_3*k2_1*k4_2 + k0_2*k1_0*k2_3*k4_2 + k0_1*k1_2*k2_3*k4_2 + k0_1*k1_3*k2_1*k4_3 + k0_2*k1_0*k2_3*k4_3 + k0_2*k1_3*k2_1*k4_2 + k0_1*k1_2*k2_3*k4_3 + k0_1*k1_3*k2_3*k4_2 + k0_2*k1_0*k2_4*k4_3 + k0_2*k1_2*k2_3*k4_2 + k0_2*k1_3*k2_1*k4_3 + k0_1*k1_2*k2_4*k4_3 + k0_1*k1_3*k2_3*k4_3 + k0_2*k1_2*k2_3*k4_3 + k0_2*k1_3*k2_3*k4_2 + k0_1*k1_3*k2_4*k4_3 + k0_2*k1_2*k2_4*k4_3 + k0_2*k1_3*k2_3*k4_3 + k0_2*k1_3*k2_4*k4_3)/(k0_2*k1_0*k2_4*k3_1 + k0_1*k1_2*k2_4*k3_1 + k0_1*k1_3*k2_0*k3_4 + k0_2*k1_0*k2_4*k3_2 + k0_1*k1_2*k2_4*k3_2 + k0_1*k1_3*k2_1*k3_4 + k0_2*k1_0*k2_3*k3_4 + k0_2*k1_2*k2_4*k3_1 + k0_1*k1_2*k2_3*k3_4 + k0_1*k1_3*k2_4*k3_2 + k0_2*k1_0*k2_4*k3_4 + k0_2*k1_2*k2_4*k3_2 + k0_2*k1_3*k2_1*k3_4 + k0_1*k1_2*k2_4*k3_4 + k0_1*k1_3*k2_3*k3_4 + k0_2*k1_2*k2_3*k3_4 + k0_2*k1_3*k2_4*k3_2 + k0_1*k1_3*k2_4*k3_4 + k0_2*k1_2*k2_4*k3_4 + k0_2*k1_3*k2_3*k3_4 + k0_2*k1_3*k2_4*k3_4 + k0_1*k1_3*k2_0*k4_2 + k0_1*k1_3*k2_0*k4_3 + k0_1*k1_3*k2_1*k4_2 + k0_2*k1_0*k2_3*k4_2 + k0_1*k1_2*k2_3*k4_2 + k0_1*k1_3*k2_1*k4_3 + k0_2*k1_0*k2_3*k4_3 + k0_2*k1_3*k2_1*k4_2 + k0_1*k1_2*k2_3*k4_3 + k0_1*k1_3*k2_3*k4_2 + k0_2*k1_0*k2_4*k4_3 + k0_2*k1_2*k2_3*k4_2 + k0_2*k1_3*k2_1*k4_3 + k0_1*k1_2*k2_4*k4_3 + k0_1*k1_3*k2_3*k4_3 + k0_2*k1_2*k2_3*k4_3 + k0_2*k1_3*k2_3*k4_2 + k0_1*k1_3*k2_4*k4_3 + k0_2*k1_2*k2_4*k4_3 + k0_2*k1_3*k2_3*k4_3 + k0_2*k1_3*k2_4*k4_3 + k0_2*k1_0*k3_1*k4_2 + k0_1*k1_2*k3_1*k4_2 + k0_2*k1_0*k3_1*k4_3 + k0_2*k1_0*k3_2*k4_2 + k0_1*k1_2*k3_1*k4_3 + k0_1*k1_2*k3_2*k4_2 + k0_2*k1_0*k3_2*k4_3 + k0_2*k1_2*k3_1*k4_2 + k0_1*k1_2*k3_2*k4_3 + k0_1*k1_3*k3_2*k4_2 + k0_2*k1_0*k3_4*k4_2 + k0_2*k1_2*k3_1*k4_3 + k0_2*k1_2*k3_2*k4_2 + k0_1*k1_2*k3_4*k4_2 + k0_1*k1_3*k3_2*k4_3 + k0_2*k1_2*k3_2*k4_3 + k0_2*k1_3*k3_2*k4_2 + k0_1*k1_3*k3_4*k4_2 + k0_2*k1_2*k3_4*k4_2 + k0_2*k1_3*k3_2*k4_3 + k0_2*k1_3*k3_4*k4_2 + k0_1*k2_0*k3_1*k4_2 + k0_1*k2_0*k3_1*k4_3 + k0_1*k2_0*k3_2*k4_2 + k0_1*k2_1*k3_1*k4_2 + k0_1*k2_0*k3_2*k4_3 + k0_1*k2_1*k3_1*k4_3 + k0_1*k2_1*k3_2*k4_2 + k0_2*k2_1*k3_1*k4_2 + k0_1*k2_0*k3_4*k4_2 + k0_1*k2_1*k3_2*k4_3 + k0_1*k2_3*k3_1*k4_2 + k0_2*k2_1*k3_1*k4_3 + k0_2*k2_1*k3_2*k4_2 + k0_1*k2_1*k3_4*k4_2 + k0_1*k2_3*k3_1*k4_3 + k0_2*k2_1*k3_2*k4_3 + k0_2*k2_3*k3_1*k4_2 + k0_1*k2_4*k3_1*k4_3 + k0_2*k2_1*k3_4*k4_2 + k0_2*k2_3*k3_1*k4_3 + k0_2*k2_4*k3_1*k4_3 + k1_0*k2_0*k3_1*k4_2 + k1_0*k2_0*k3_1*k4_3 + k1_0*k2_0*k3_2*k4_2 + k1_0*k2_1*k3_1*k4_2 + k1_0*k2_0*k3_2*k4_3 + k1_0*k2_1*k3_1*k4_3 + k1_0*k2_1*k3_2*k4_2 + k1_2*k2_0*k3_1*k4_2 + k1_0*k2_0*k3_4*k4_2 + k1_0*k2_1*k3_2*k4_3 + k1_0*k2_3*k3_1*k4_2 + k1_2*k2_0*k3_1*k4_3 + k1_2*k2_0*k3_2*k4_2 + k1_0*k2_1*k3_4*k4_2 + k1_0*k2_3*k3_1*k4_3 + k1_2*k2_0*k3_2*k4_3 + k1_3*k2_0*k3_2*k4_2 + k1_0*k2_4*k3_1*k4_3 + k1_2*k2_0*k3_4*k4_2 + k1_3*k2_0*k3_2*k4_3 + k1_3*k2_0*k3_4*k4_2);
S4 = 1 - S0 - S1 - S2 - S3;
% Turnover Rate
JNADH = ETC1Activity*(kfNADH_02*S0 - krNADH_02*S2 ...
+ kfNADH_24*S2 - krNADH_24*S4 ...
+ kfNADH_13*S1 - krNADH_13*S3);
% NADH oxidation releases 1 proton in matrix space
% JSO = ETC1Activity*(S1.*(kfSO_10 - krSO_21) + S2.*(kfSO_21 - krSO_32) + S3.*(kfSO_32 - krSO_43) + S4.*kfSO_43 - krSO_10.*S0);
% JH2O2 = ETC1Activity*(S2.*(kfH2O2_20 - krH2O2_42) + S3.*kfH2O2_31 + S4.*kfH2O2_42 - krH2O2_31.*S1 - krH2O2_20.*S0);
J_ETC1_im_to_matrix = JNADH;
%ETC3:im_to_matrix
ETC3_activity =par(7);
dPsi = DPsi_im_to_matrix;
Hp = h_im;
Hn = h_matrix;
Q1 = coQ_matrix;
QH21 = coQH2_matrix;
c2 = cytocred_im;
c3 = cytocox_im;
O2 = 200e-6; % constant oxygen concentration (M)
SO = 0;
% R = 8.314e-3;
% F = 0.0965;
% Set KDQH2 at Qi site
CIII_KQH2o = 8e-4;
CIII_KQi = 1e-3;
CIII_Kc3 = 1.1193e-06;
CIII_KQo = 5e-4;
CIII_Kc2 = 1.1666e-06;
CIII_KQH2i = CIII_KQH2o^2*CIII_KQi*CIII_Kc3^2/CIII_KQo^2/CIII_Kc2^2;
% MR constraint
% Update midpoint potentials from thermodynamic data
dGf_UQH2 = -19.1505595;
dGf_UQ = 69.2309743;
dGf_c2 = -27.586932;
dGf_c3 = -6.92309743;
dGf_SO = 12.41076695;
dGf_O2 = 16.4;
CIII_Em0_Q_QH2 = (dGf_UQH2 - dGf_UQ)/(-2*F);
CIII_Em0_c = (dGf_c2 - dGf_c3)/(-1*F);
CIII_Em0_SO = (dGf_SO - dGf_O2)/(-1*F);
% Update Q Thermodynamics
CIII_Kstabo = 1e-9;
CIII_Kstabi = 0.0078;
CIII_Em0_Q_SQo = CIII_Em0_Q_QH2 + RT/F/2*log(CIII_Kstabo*1e-14);
% mV
CIII_Em0_SQ_QH2o = 2*CIII_Em0_Q_QH2 - CIII_Em0_Q_SQo;
% mV (200 - 300)
CIII_Em0_Q_SQi = CIII_Em0_Q_QH2 + RT/F/2*log(CIII_Kstabi*1e-14);
% mV (assumes Kstabi at pH 7)
CIII_Em0_SQ_QH2i = 2*CIII_Em0_Q_QH2 - CIII_Em0_Q_SQi;
% mV (16-150)
% Binding Polynomials for Protonated States
CIII_pK_ISPox2 = 9.16;
CIII_pK_ISPox1 = 7.63;
CIII_pK_bLox = 5.9;
CIII_pK_bLred = 7.9;
CIII_pK_bHox = 5.7;
CIII_pK_bHred = 7.7;
CIII_pK_QH = 13.2;
CIII_pK_QH2 = 11.3;
CIII_pK_SO = 4.7;
P_ISP = (1 + Hp/10^-CIII_pK_ISPox2 + Hp^2/10^-CIII_pK_ISPox2/10^-CIII_pK_ISPox1);
% P_bLox = (1 + Hp/10^-CIII_pK_bLox);
% P_bLred = (1 + Hp/10^-CIII_pK_bLred);
% P_bHox = (1 + Hn/10^-CIII_pK_bHox);
% P_bHred = (1 + Hn/10^-CIII_pK_bHred);
P_QH2o = (1 + Hp/10^-CIII_pK_QH + Hp^2/10^-CIII_pK_QH/10^-CIII_pK_QH2);
P_QH2i = (1 + Hn/10^-CIII_pK_QH + Hn^2/10^-CIII_pK_QH/10^-CIII_pK_QH2);
P_SO = (1 + Hp/10^-CIII_pK_SO);
% Binding Polynomials for enzyme, substrates, products and regulators
% Qo-site
CIII_AA = 0;
CIII_KAA = 1e-10;
P_Qo = (1 + Q1/CIII_KQo + QH21/CIII_KQH2o);
% Cytc c-site
P_c = (1 + c2/CIII_Kc2 + c3/CIII_Kc3);
% Qi-site
P_Qi = (1 + Q1/CIII_KQi + QH21/CIII_KQH2i + CIII_AA/CIII_KAA);
% Midpoint potentials
% Qo-site
% CIII_Em0_ISP = 311;
CIII_Em0_bL = 39;
Em_c = CIII_Em0_c;
% pH independent
% Em_Q_QH2o = CIII_Em0_Q_QH2 - log(10)*RT/F/2*log10(Hp^2/10^-CIII_pK_QH2/10^-CIII_pK_QH/P_QH2o/Hp^2);
% assuming linked to P-side
Em_Q_SQo = CIII_Em0_Q_SQo;
% pH independent
Em_SQ_QH2o = CIII_Em0_SQ_QH2o - log(10)*RT/F*log10((Hp^2/10^-CIII_pK_QH2/10^-CIII_pK_QH)/P_QH2o/Hp^2);
% Em_ISP = CIII_Em0_ISP - log(10)*RT/F*log10(P_ISP/(Hp^2/10^-CIII_pK_ISPox2/10^-CIII_pK_ISPox1));
Em_bL = CIII_Em0_bL - log(10)*RT/F*log10((Hp+10^-CIII_pK_bLox)/(Hp+10^-CIII_pK_bLred));
% assuming linked to P-side (Izrailev et al. 1999, also Crofts)
% Qi-site
CIII_Em0_bH = 160;
Em_bH = CIII_Em0_bH - log(10)*RT/F*log10((Hn+10^-CIII_pK_bHox)/(Hn+10^-CIII_pK_bHred));
% assuming linked to N-side, log10((Hn/10^-CIII_pK_bHred/P_bHred)/(Hn/10^-CIII_pK_bHox/P_bHox))
% Em_Q_QH2i = CIII_Em0_Q_QH2 - log(10)*RT/F/2*log10(Hn^2/10^-CIII_pK_QH2/10^-CIII_pK_QH/P_QH2i/Hn^2);
% assuming linked to N-siden
Em_Q_SQi = CIII_Em0_Q_SQi;
% pH independent
Em_SQ_QH2i = CIII_Em0_SQ_QH2i - log(10)*RT/F*log10((Hn^2/10^-CIII_pK_QH2/10^-CIII_pK_QH)/P_QH2i/Hn^2);
% assuming linked to N-side
% Superoxide
Em_SO = CIII_Em0_SO - log(10)*RT/F*log10(1/P_SO/10^-CIII_pK_SO);
% State Transition Thermodynamics
% QH2o<->SQo
dGo = -F*(Em_c - Em_SQ_QH2o);
% dG1 = dGo + RT*log((E2*Q1*c2*Hp^2)/(QH21*c3*E1*1e-14)*(PE1*CIII_KQH2o*CIII_Kc3)/(PE2*CIII_KQo*CIII_Kc2));
% bHred<->bLred
dGi = -F*(Em_bH - Em_bL);
% dG2 = dGi + RT*log((E3)/(E2*Q1)*(PE2*CIII_KQi)/(PE3)) + F*dPsi;
% Substates via Rapid Equilibrium
% ISP protonation state
PISP = Hp/10^-CIII_pK_ISPox2*(1+Hp/10^-CIII_pK_ISPox1)/P_ISP;
% bL reduction by SQo
dGQo = -F*(Em_bL - Em_Q_SQo);
KQo = exp(-dGQo/RT);
% Q-bLred / SQ-bLox
% State 2 fractional species polynomial
% fQoa = 1;
% SQo
rQo = KQo*P_Qo*CIII_KQo/Q1;
% bL_red
PQo = 1 + rQo;
% SQo
% Q reduction by bHred
dGQi = -F*(Em_Q_SQi - Em_bH);
KQi = exp(-dGQi/RT);
% bLox-SQ / bLred-Q
% State 3 fractional species polynomial
% fQia = 1;
% bHred
rQi = KQi*Q1/CIII_KQi/P_Qi;
% bHox-SQi
PQi = 1 + rQi;
% bHred
% SQi reduction by bHred
CIII_beta2 = 0.5;
dGQi2 = -F*(Em_SQ_QH2i - Em_bH) + CIII_beta2*2*F*dPsi;
% protons included in Em_SQ_QH2i
KQi2 = exp(-dGQi2/RT);
% fQi2a = 1;
% bHred-SQi
rQi2 = KQi2*P_Qi*CIII_KQH2i/QH21;
% bHox
PQi2 = 1 + rQi2;
%bHred-SQi
% Superoxide production/consumption rates
% SO thermodynamics
CIII_kSO = 1600;
dG_SO = -F*(Em_SO - Em_Q_SQo);
kSOr = CIII_kSO*exp(dG_SO/RT);
% Reverse Rate Constants
CIII_k120 = 3.0728e+03;
CIII_k230 = 33100000;
CIII_k340 = 3.0728e+03;
CIII_k410 = 14370;
CIII_k450 = 3;
CIII_k520 = 14370;
CIII_k260 = 3.0728e+03;
CIII_k640 = 33100000;
k210 = CIII_k120*exp(dGo/RT);
k320 = CIII_k230*exp(dGi/RT);
k430 = CIII_k340*exp(dGo/RT);
k140 = CIII_k410*exp(dGi/RT);
k540 = CIII_k450*exp(dGo/RT);
k250 = CIII_k520*exp(dGi/RT);
k620 = CIII_k260*exp(dGo/RT);
k460 = CIII_k640*exp(dGi/RT);
% State Transition Rates
% Coupled QH2 oxidation /c3 reduction
kfQH2c3_12 = CIII_k120*QH21/CIII_KQH2o/P_Qo*c3/CIII_Kc3/P_c*PISP;
krQH2c3_12 = k210*c2/CIII_Kc2/P_c*PISP/PQo;
kfQH2c3_34 = CIII_k340*QH21/CIII_KQH2o/P_Qo*c3/CIII_Kc3/P_c*PISP;
krQH2c3_34 = k430*c2/CIII_Kc2/P_c*PISP/PQo;
kfQH2c3_45 = CIII_k450*QH21/CIII_KQH2o/P_Qo*c3/CIII_Kc3/P_c*PISP*rQo/PQo;
krQH2c3_45 = k540*c2/CIII_Kc2/P_c*PISP;
kfQH2c3_26 = CIII_k260*QH21/CIII_KQH2o/P_Qo*c3/CIII_Kc3/P_c*PISP*rQo/PQo;
krQH2c3_26 = k620*c2/CIII_Kc2/P_c*PISP;
% Q uptake at SQ site
CIII_beta1 = 0.5;
kfQ_23 = CIII_k230*rQo/PQo*rQi2/PQi2*exp(-CIII_beta1*F*dPsi/RT/2);
krQ_23 = k320/PQi*exp(CIII_beta1*F*dPsi/RT/2);
kfQ_64 = CIII_k640*rQi2/PQi2*exp(-CIII_beta1*F*dPsi/RT/2);
krQ_64 = k460/PQo/PQi*exp(CIII_beta1*F*dPsi/RT/2);
% QH2 regeneration
kfQH2_41 = CIII_k410*rQo/PQo*rQi/PQi*exp(-CIII_beta1*F*dPsi/RT/2);
krQH2_41 = k140/PQi2*exp(CIII_beta1*F*dPsi/RT/2);
kfQH2_52 = CIII_k520*rQi/PQi*exp(-CIII_beta1*F*dPsi/RT/2);
krQH2_52 = k250/PQo/PQi2*exp(CIII_beta1*F*dPsi/RT/2);
% SO production
kfSO_21 = CIII_kSO*O2/PQo;
krSO_21 = kSOr*SO*Q1/CIII_KQo/P_Qo;
kfSO_43 = CIII_kSO*O2/PQo;
krSO_43 = kSOr*SO*Q1/CIII_KQo/P_Qo;
kfSO_54 = CIII_kSO*O2;
krSO_54 = kSOr*SO*rQo/PQo*Q1/CIII_KQo/P_Qo;
kfSO_62 = CIII_kSO*O2;
krSO_62 = kSOr*SO*rQo/PQo*Q1/CIII_KQo/P_Qo;
% QH2o + c3 -> SQo + c2 + 2H+
k12 = kfQH2c3_12 + krSO_21;
k21 = krQH2c3_12 + kfSO_21;
% bL -> bH
k23 = kfQ_23;
k32 = krQ_23;
% QH2o + c3 -> SQo + c2 + 2H+
k34 = kfQH2c3_34 + krSO_43;
k43 = krQH2c3_34 + kfSO_43;
% bL -> bH
k41 = kfQH2_41;
k14 = krQH2_41;
% QH2o + c3 -> SQo + c2 + 2H+
k45 = kfQH2c3_45 + krSO_54;
k54 = krQH2c3_45 + kfSO_54;
% bL -> bH
k52 = kfQH2_52;
k25 = krQH2_52;
% QH2o + c3 -> SQo + c2 + 2H+
k26 = kfQH2c3_26 + krSO_62;
k62 = krQH2c3_26 + kfSO_62;
% bL -> bH
k64 = kfQ_64;
k46 = krQ_64;
% Steady-State Fractional Occupancies (solved analytically)
E1=(k21*k32*k41*k52*k62 + k21*k32*k41*k52*k64 + k21*k32*k41*k54*k62 + k21*k32*k43*k52*k62 + k21*k34*k41*k52*k62 + k21*k32*k41*k54*k64 + k21*k32*k43*k52*k64 + k21*k32*k43*k54*k62 + k21*k32*k45*k52*k62 + k21*k34*k41*k52*k64 + k21*k34*k41*k54*k62 + k23*k34*k41*k52*k62 + k21*k32*k46*k52*k62 + k21*k32*k43*k54*k64 + k21*k32*k45*k52*k64 + k21*k34*k41*k54*k64 + k21*k34*k45*k52*k62 + k23*k34*k41*k52*k64 + k23*k34*k41*k54*k62 + k25*k32*k41*k54*k62 + k21*k32*k46*k54*k62 + k21*k34*k46*k52*k62 + k26*k32*k41*k52*k64 + k21*k34*k45*k52*k64 + k23*k34*k41*k54*k64 + k25*k32*k41*k54*k64 + k25*k34*k41*k54*k62 + k21*k34*k46*k54*k62 + k26*k32*k41*k54*k64 + k26*k34*k41*k52*k64 + k25*k34*k41*k54*k64 + k26*k34*k41*k54*k64)/(k12*k26*k32*k41*k52 + k12*k26*k32*k41*k54 + k12*k26*k32*k43*k52 + k12*k26*k34*k41*k52 + k14*k21*k32*k46*k52 + k12*k23*k34*k46*k52 + k12*k26*k32*k43*k54 + k12*k26*k32*k45*k52 + k12*k26*k34*k41*k54 + k14*k21*k32*k46*k54 + k14*k21*k34*k46*k52 + k14*k26*k32*k43*k52 + k12*k26*k32*k46*k52 + k12*k23*k34*k46*k54 + k12*k25*k32*k46*k54 + k12*k26*k34*k45*k52 + k14*k21*k34*k46*k54 + k14*k23*k34*k46*k52 + k14*k26*k32*k43*k54 + k14*k26*k32*k45*k52 + k12*k26*k32*k46*k54 + k12*k26*k34*k46*k52 + k14*k26*k32*k46*k52 + k12*k25*k34*k46*k54 + k14*k23*k34*k46*k54 + k14*k25*k32*k46*k54 + k14*k26*k34*k45*k52 + k12*k25*k32*k41*k62 + k12*k26*k34*k46*k54 + k14*k26*k32*k46*k54 + k14*k26*k34*k46*k52 + k14*k25*k34*k46*k54 + k12*k25*k32*k41*k64 + k12*k25*k32*k43*k62 + k12*k25*k34*k41*k62 + k14*k21*k32*k45*k62 + k14*k26*k34*k46*k54 + k12*k23*k34*k45*k62 + k12*k25*k32*k43*k64 + k12*k25*k32*k45*k62 + k12*k25*k34*k41*k64 + k14*k21*k32*k45*k64 + k14*k21*k34*k45*k62 + k14*k25*k32*k43*k62 + k12*k25*k32*k46*k62 + k12*k23*k34*k45*k64 + k12*k25*k32*k45*k64 + k12*k25*k34*k45*k62 + k14*k21*k34*k45*k64 + k14*k23*k34*k45*k62 + k14*k25*k32*k43*k64 + k14*k25*k32*k45*k62 + k12*k25*k34*k46*k62 + k12*k26*k32*k45*k64 + k14*k25*k32*k46*k62 + k12*k25*k34*k45*k64 + k14*k23*k34*k45*k64 + k14*k25*k32*k45*k64 + k14*k25*k34*k45*k62 + k12*k26*k34*k45*k64 + k14*k21*k32*k52*k62 + k14*k25*k34*k46*k62 + k14*k26*k32*k45*k64 + k14*k25*k34*k45*k64 + k12*k23*k34*k52*k62 + k14*k21*k32*k52*k64 + k14*k21*k32*k54*k62 + k14*k21*k34*k52*k62 + k14*k26*k34*k45*k64 + k12*k23*k34*k52*k64 + k12*k23*k34*k54*k62 + k12*k25*k32*k54*k62 + k14*k21*k32*k54*k64 + k14*k21*k34*k52*k64 + k14*k21*k34*k54*k62 + k14*k23*k34*k52*k62 + k12*k26*k32*k52*k64 + k12*k23*k34*k54*k64 + k12*k25*k32*k54*k64 + k12*k25*k34*k54*k62 + k14*k21*k34*k54*k64 + k14*k23*k34*k52*k64 + k14*k23*k34*k54*k62 + k14*k25*k32*k54*k62 + k12*k26*k32*k54*k64 + k12*k26*k34*k52*k64 + k14*k26*k32*k52*k64 + k12*k25*k34*k54*k64 + k14*k23*k34*k54*k64 + k14*k25*k32*k54*k64 + k14*k25*k34*k54*k62 + k12*k23*k41*k52*k62 + k12*k26*k34*k54*k64 + k14*k26*k32*k54*k64 + k14*k26*k34*k52*k64 + k14*k25*k34*k54*k64 + k12*k23*k41*k52*k64 + k12*k23*k41*k54*k62 + k12*k23*k43*k52*k62 + k14*k21*k43*k52*k62 + k14*k26*k34*k54*k64 + k12*k23*k41*k54*k64 + k12*k23*k43*k52*k64 + k12*k23*k43*k54*k62 + k12*k23*k45*k52*k62 + k14*k21*k43*k52*k64 + k14*k21*k43*k54*k62 + k14*k23*k43*k52*k62 + k12*k23*k46*k52*k62 + k12*k23*k43*k54*k64 + k12*k23*k45*k52*k64 + k12*k25*k43*k54*k62 + k14*k21*k43*k54*k64 + k14*k23*k43*k52*k64 + k14*k23*k43*k54*k62 + k14*k23*k45*k52*k62 + k12*k23*k46*k54*k62 + k12*k26*k43*k52*k64 + k14*k23*k46*k52*k62 + k12*k25*k43*k54*k64 + k14*k23*k43*k54*k64 + k14*k23*k45*k52*k64 + k14*k25*k43*k54*k62 + k12*k26*k43*k54*k64 + k12*k32*k41*k52*k62 + k14*k23*k46*k54*k62 + k14*k26*k43*k52*k64 + k14*k25*k43*k54*k64 + k12*k32*k41*k52*k64 + k12*k32*k41*k54*k62 + k12*k32*k43*k52*k62 + k12*k34*k41*k52*k62 + k14*k26*k43*k54*k64 + k12*k32*k41*k54*k64 + k12*k32*k43*k52*k64 + k12*k32*k43*k54*k62 + k12*k32*k45*k52*k62 + k12*k34*k41*k52*k64 + k12*k34*k41*k54*k62 + k14*k32*k43*k52*k62 + k12*k32*k46*k52*k62 + k12*k32*k43*k54*k64 + k12*k32*k45*k52*k64 + k12*k34*k41*k54*k64 + k12*k34*k45*k52*k62 + k14*k32*k43*k52*k64 + k14*k32*k43*k54*k62 + k14*k32*k45*k52*k62 + k12*k32*k46*k54*k62 + k12*k34*k46*k52*k62 + k14*k32*k46*k52*k62 + k12*k34*k45*k52*k64 + k14*k32*k43*k54*k64 + k14*k32*k45*k52*k64 + k14*k34*k45*k52*k62 + k12*k34*k46*k54*k62 + k14*k32*k46*k54*k62 + k14*k34*k46*k52*k62 + k21*k32*k41*k52*k62 + k14*k34*k45*k52*k64 + k14*k34*k46*k54*k62 + k21*k32*k41*k52*k64 + k21*k32*k41*k54*k62 + k21*k32*k43*k52*k62 + k21*k34*k41*k52*k62 + k21*k32*k41*k54*k64 + k21*k32*k43*k52*k64 + k21*k32*k43*k54*k62 + k21*k32*k45*k52*k62 + k21*k34*k41*k52*k64 + k21*k34*k41*k54*k62 + k23*k34*k41*k52*k62 + k21*k32*k46*k52*k62 + k21*k32*k43*k54*k64 + k21*k32*k45*k52*k64 + k21*k34*k41*k54*k64 + k21*k34*k45*k52*k62 + k23*k34*k41*k52*k64 + k23*k34*k41*k54*k62 + k25*k32*k41*k54*k62 + k21*k32*k46*k54*k62 + k21*k34*k46*k52*k62 + k26*k32*k41*k52*k64 + k21*k34*k45*k52*k64 + k23*k34*k41*k54*k64 + k25*k32*k41*k54*k64 + k25*k34*k41*k54*k62 + k21*k34*k46*k54*k62 + k26*k32*k41*k54*k64 + k26*k34*k41*k52*k64 + k25*k34*k41*k54*k64 + k26*k34*k41*k54*k64);
E2=(k12*k32*k41*k52*k62 + k12*k32*k41*k52*k64 + k12*k32*k41*k54*k62 + k12*k32*k43*k52*k62 + k12*k34*k41*k52*k62 + k12*k32*k41*k54*k64 + k12*k32*k43*k52*k64 + k12*k32*k43*k54*k62 + k12*k32*k45*k52*k62 + k12*k34*k41*k52*k64 + k12*k34*k41*k54*k62 + k14*k32*k43*k52*k62 + k12*k32*k46*k52*k62 + k12*k32*k43*k54*k64 + k12*k32*k45*k52*k64 + k12*k34*k41*k54*k64 + k12*k34*k45*k52*k62 + k14*k32*k43*k52*k64 + k14*k32*k43*k54*k62 + k14*k32*k45*k52*k62 + k12*k32*k46*k54*k62 + k12*k34*k46*k52*k62 + k14*k32*k46*k52*k62 + k12*k34*k45*k52*k64 + k14*k32*k43*k54*k64 + k14*k32*k45*k52*k64 + k14*k34*k45*k52*k62 + k12*k34*k46*k54*k62 + k14*k32*k46*k54*k62 + k14*k34*k46*k52*k62 + k14*k34*k45*k52*k64 + k14*k34*k46*k54*k62)/(k12*k26*k32*k41*k52 + k12*k26*k32*k41*k54 + k12*k26*k32*k43*k52 + k12*k26*k34*k41*k52 + k14*k21*k32*k46*k52 + k12*k23*k34*k46*k52 + k12*k26*k32*k43*k54 + k12*k26*k32*k45*k52 + k12*k26*k34*k41*k54 + k14*k21*k32*k46*k54 + k14*k21*k34*k46*k52 + k14*k26*k32*k43*k52 + k12*k26*k32*k46*k52 + k12*k23*k34*k46*k54 + k12*k25*k32*k46*k54 + k12*k26*k34*k45*k52 + k14*k21*k34*k46*k54 + k14*k23*k34*k46*k52 + k14*k26*k32*k43*k54 + k14*k26*k32*k45*k52 + k12*k26*k32*k46*k54 + k12*k26*k34*k46*k52 + k14*k26*k32*k46*k52 + k12*k25*k34*k46*k54 + k14*k23*k34*k46*k54 + k14*k25*k32*k46*k54 + k14*k26*k34*k45*k52 + k12*k25*k32*k41*k62 + k12*k26*k34*k46*k54 + k14*k26*k32*k46*k54 + k14*k26*k34*k46*k52 + k14*k25*k34*k46*k54 + k12*k25*k32*k41*k64 + k12*k25*k32*k43*k62 + k12*k25*k34*k41*k62 + k14*k21*k32*k45*k62 + k14*k26*k34*k46*k54 + k12*k23*k34*k45*k62 + k12*k25*k32*k43*k64 + k12*k25*k32*k45*k62 + k12*k25*k34*k41*k64 + k14*k21*k32*k45*k64 + k14*k21*k34*k45*k62 + k14*k25*k32*k43*k62 + k12*k25*k32*k46*k62 + k12*k23*k34*k45*k64 + k12*k25*k32*k45*k64 + k12*k25*k34*k45*k62 + k14*k21*k34*k45*k64 + k14*k23*k34*k45*k62 + k14*k25*k32*k43*k64 + k14*k25*k32*k45*k62 + k12*k25*k34*k46*k62 + k12*k26*k32*k45*k64 + k14*k25*k32*k46*k62 + k12*k25*k34*k45*k64 + k14*k23*k34*k45*k64 + k14*k25*k32*k45*k64 + k14*k25*k34*k45*k62 + k12*k26*k34*k45*k64 + k14*k21*k32*k52*k62 + k14*k25*k34*k46*k62 + k14*k26*k32*k45*k64 + k14*k25*k34*k45*k64 + k12*k23*k34*k52*k62 + k14*k21*k32*k52*k64 + k14*k21*k32*k54*k62 + k14*k21*k34*k52*k62 + k14*k26*k34*k45*k64 + k12*k23*k34*k52*k64 + k12*k23*k34*k54*k62 + k12*k25*k32*k54*k62 + k14*k21*k32*k54*k64 + k14*k21*k34*k52*k64 + k14*k21*k34*k54*k62 + k14*k23*k34*k52*k62 + k12*k26*k32*k52*k64 + k12*k23*k34*k54*k64 + k12*k25*k32*k54*k64 + k12*k25*k34*k54*k62 + k14*k21*k34*k54*k64 + k14*k23*k34*k52*k64 + k14*k23*k34*k54*k62 + k14*k25*k32*k54*k62 + k12*k26*k32*k54*k64 + k12*k26*k34*k52*k64 + k14*k26*k32*k52*k64 + k12*k25*k34*k54*k64 + k14*k23*k34*k54*k64 + k14*k25*k32*k54*k64 + k14*k25*k34*k54*k62 + k12*k23*k41*k52*k62 + k12*k26*k34*k54*k64 + k14*k26*k32*k54*k64 + k14*k26*k34*k52*k64 + k14*k25*k34*k54*k64 + k12*k23*k41*k52*k64 + k12*k23*k41*k54*k62 + k12*k23*k43*k52*k62 + k14*k21*k43*k52*k62 + k14*k26*k34*k54*k64 + k12*k23*k41*k54*k64 + k12*k23*k43*k52*k64 + k12*k23*k43*k54*k62 + k12*k23*k45*k52*k62 + k14*k21*k43*k52*k64 + k14*k21*k43*k54*k62 + k14*k23*k43*k52*k62 + k12*k23*k46*k52*k62 + k12*k23*k43*k54*k64 + k12*k23*k45*k52*k64 + k12*k25*k43*k54*k62 + k14*k21*k43*k54*k64 + k14*k23*k43*k52*k64 + k14*k23*k43*k54*k62 + k14*k23*k45*k52*k62 + k12*k23*k46*k54*k62 + k12*k26*k43*k52*k64 + k14*k23*k46*k52*k62 + k12*k25*k43*k54*k64 + k14*k23*k43*k54*k64 + k14*k23*k45*k52*k64 + k14*k25*k43*k54*k62 + k12*k26*k43*k54*k64 + k12*k32*k41*k52*k62 + k14*k23*k46*k54*k62 + k14*k26*k43*k52*k64 + k14*k25*k43*k54*k64 + k12*k32*k41*k52*k64 + k12*k32*k41*k54*k62 + k12*k32*k43*k52*k62 + k12*k34*k41*k52*k62 + k14*k26*k43*k54*k64 + k12*k32*k41*k54*k64 + k12*k32*k43*k52*k64 + k12*k32*k43*k54*k62 + k12*k32*k45*k52*k62 + k12*k34*k41*k52*k64 + k12*k34*k41*k54*k62 + k14*k32*k43*k52*k62 + k12*k32*k46*k52*k62 + k12*k32*k43*k54*k64 + k12*k32*k45*k52*k64 + k12*k34*k41*k54*k64 + k12*k34*k45*k52*k62 + k14*k32*k43*k52*k64 + k14*k32*k43*k54*k62 + k14*k32*k45*k52*k62 + k12*k32*k46*k54*k62 + k12*k34*k46*k52*k62 + k14*k32*k46*k52*k62 + k12*k34*k45*k52*k64 + k14*k32*k43*k54*k64 + k14*k32*k45*k52*k64 + k14*k34*k45*k52*k62 + k12*k34*k46*k54*k62 + k14*k32*k46*k54*k62 + k14*k34*k46*k52*k62 + k21*k32*k41*k52*k62 + k14*k34*k45*k52*k64 + k14*k34*k46*k54*k62 + k21*k32*k41*k52*k64 + k21*k32*k41*k54*k62 + k21*k32*k43*k52*k62 + k21*k34*k41*k52*k62 + k21*k32*k41*k54*k64 + k21*k32*k43*k52*k64 + k21*k32*k43*k54*k62 + k21*k32*k45*k52*k62 + k21*k34*k41*k52*k64 + k21*k34*k41*k54*k62 + k23*k34*k41*k52*k62 + k21*k32*k46*k52*k62 + k21*k32*k43*k54*k64 + k21*k32*k45*k52*k64 + k21*k34*k41*k54*k64 + k21*k34*k45*k52*k62 + k23*k34*k41*k52*k64 + k23*k34*k41*k54*k62 + k25*k32*k41*k54*k62 + k21*k32*k46*k54*k62 + k21*k34*k46*k52*k62 + k26*k32*k41*k52*k64 + k21*k34*k45*k52*k64 + k23*k34*k41*k54*k64 + k25*k32*k41*k54*k64 + k25*k34*k41*k54*k62 + k21*k34*k46*k54*k62 + k26*k32*k41*k54*k64 + k26*k34*k41*k52*k64 + k25*k34*k41*k54*k64 + k26*k34*k41*k54*k64);
E3=(k12*k23*k41*k52*k62 + k12*k23*k41*k52*k64 + k12*k23*k41*k54*k62 + k12*k23*k43*k52*k62 + k14*k21*k43*k52*k62 + k12*k23*k41*k54*k64 + k12*k23*k43*k52*k64 + k12*k23*k43*k54*k62 + k12*k23*k45*k52*k62 + k14*k21*k43*k52*k64 + k14*k21*k43*k54*k62 + k14*k23*k43*k52*k62 + k12*k23*k46*k52*k62 + k12*k23*k43*k54*k64 + k12*k23*k45*k52*k64 + k12*k25*k43*k54*k62 + k14*k21*k43*k54*k64 + k14*k23*k43*k52*k64 + k14*k23*k43*k54*k62 + k14*k23*k45*k52*k62 + k12*k23*k46*k54*k62 + k12*k26*k43*k52*k64 + k14*k23*k46*k52*k62 + k12*k25*k43*k54*k64 + k14*k23*k43*k54*k64 + k14*k23*k45*k52*k64 + k14*k25*k43*k54*k62 + k12*k26*k43*k54*k64 + k14*k23*k46*k54*k62 + k14*k26*k43*k52*k64 + k14*k25*k43*k54*k64 + k14*k26*k43*k54*k64)/(k12*k26*k32*k41*k52 + k12*k26*k32*k41*k54 + k12*k26*k32*k43*k52 + k12*k26*k34*k41*k52 + k14*k21*k32*k46*k52 + k12*k23*k34*k46*k52 + k12*k26*k32*k43*k54 + k12*k26*k32*k45*k52 + k12*k26*k34*k41*k54 + k14*k21*k32*k46*k54 + k14*k21*k34*k46*k52 + k14*k26*k32*k43*k52 + k12*k26*k32*k46*k52 + k12*k23*k34*k46*k54 + k12*k25*k32*k46*k54 + k12*k26*k34*k45*k52 + k14*k21*k34*k46*k54 + k14*k23*k34*k46*k52 + k14*k26*k32*k43*k54 + k14*k26*k32*k45*k52 + k12*k26*k32*k46*k54 + k12*k26*k34*k46*k52 + k14*k26*k32*k46*k52 + k12*k25*k34*k46*k54 + k14*k23*k34*k46*k54 + k14*k25*k32*k46*k54 + k14*k26*k34*k45*k52 + k12*k25*k32*k41*k62 + k12*k26*k34*k46*k54 + k14*k26*k32*k46*k54 + k14*k26*k34*k46*k52 + k14*k25*k34*k46*k54 + k12*k25*k32*k41*k64 + k12*k25*k32*k43*k62 + k12*k25*k34*k41*k62 + k14*k21*k32*k45*k62 + k14*k26*k34*k46*k54 + k12*k23*k34*k45*k62 + k12*k25*k32*k43*k64 + k12*k25*k32*k45*k62 + k12*k25*k34*k41*k64 + k14*k21*k32*k45*k64 + k14*k21*k34*k45*k62 + k14*k25*k32*k43*k62 + k12*k25*k32*k46*k62 + k12*k23*k34*k45*k64 + k12*k25*k32*k45*k64 + k12*k25*k34*k45*k62 + k14*k21*k34*k45*k64 + k14*k23*k34*k45*k62 + k14*k25*k32*k43*k64 + k14*k25*k32*k45*k62 + k12*k25*k34*k46*k62 + k12*k26*k32*k45*k64 + k14*k25*k32*k46*k62 + k12*k25*k34*k45*k64 + k14*k23*k34*k45*k64 + k14*k25*k32*k45*k64 + k14*k25*k34*k45*k62 + k12*k26*k34*k45*k64 + k14*k21*k32*k52*k62 + k14*k25*k34*k46*k62 + k14*k26*k32*k45*k64 + k14*k25*k34*k45*k64 + k12*k23*k34*k52*k62 + k14*k21*k32*k52*k64 + k14*k21*k32*k54*k62 + k14*k21*k34*k52*k62 + k14*k26*k34*k45*k64 + k12*k23*k34*k52*k64 + k12*k23*k34*k54*k62 + k12*k25*k32*k54*k62 + k14*k21*k32*k54*k64 + k14*k21*k34*k52*k64 + k14*k21*k34*k54*k62 + k14*k23*k34*k52*k62 + k12*k26*k32*k52*k64 + k12*k23*k34*k54*k64 + k12*k25*k32*k54*k64 + k12*k25*k34*k54*k62 + k14*k21*k34*k54*k64 + k14*k23*k34*k52*k64 + k14*k23*k34*k54*k62 + k14*k25*k32*k54*k62 + k12*k26*k32*k54*k64 + k12*k26*k34*k52*k64 + k14*k26*k32*k52*k64 + k12*k25*k34*k54*k64 + k14*k23*k34*k54*k64 + k14*k25*k32*k54*k64 + k14*k25*k34*k54*k62 + k12*k23*k41*k52*k62 + k12*k26*k34*k54*k64 + k14*k26*k32*k54*k64 + k14*k26*k34*k52*k64 + k14*k25*k34*k54*k64 + k12*k23*k41*k52*k64 + k12*k23*k41*k54*k62 + k12*k23*k43*k52*k62 + k14*k21*k43*k52*k62 + k14*k26*k34*k54*k64 + k12*k23*k41*k54*k64 + k12*k23*k43*k52*k64 + k12*k23*k43*k54*k62 + k12*k23*k45*k52*k62 + k14*k21*k43*k52*k64 + k14*k21*k43*k54*k62 + k14*k23*k43*k52*k62 + k12*k23*k46*k52*k62 + k12*k23*k43*k54*k64 + k12*k23*k45*k52*k64 + k12*k25*k43*k54*k62 + k14*k21*k43*k54*k64 + k14*k23*k43*k52*k64 + k14*k23*k43*k54*k62 + k14*k23*k45*k52*k62 + k12*k23*k46*k54*k62 + k12*k26*k43*k52*k64 + k14*k23*k46*k52*k62 + k12*k25*k43*k54*k64 + k14*k23*k43*k54*k64 + k14*k23*k45*k52*k64 + k14*k25*k43*k54*k62 + k12*k26*k43*k54*k64 + k12*k32*k41*k52*k62 + k14*k23*k46*k54*k62 + k14*k26*k43*k52*k64 + k14*k25*k43*k54*k64 + k12*k32*k41*k52*k64 + k12*k32*k41*k54*k62 + k12*k32*k43*k52*k62 + k12*k34*k41*k52*k62 + k14*k26*k43*k54*k64 + k12*k32*k41*k54*k64 + k12*k32*k43*k52*k64 + k12*k32*k43*k54*k62 + k12*k32*k45*k52*k62 + k12*k34*k41*k52*k64 + k12*k34*k41*k54*k62 + k14*k32*k43*k52*k62 + k12*k32*k46*k52*k62 + k12*k32*k43*k54*k64 + k12*k32*k45*k52*k64 + k12*k34*k41*k54*k64 + k12*k34*k45*k52*k62 + k14*k32*k43*k52*k64 + k14*k32*k43*k54*k62 + k14*k32*k45*k52*k62 + k12*k32*k46*k54*k62 + k12*k34*k46*k52*k62 + k14*k32*k46*k52*k62 + k12*k34*k45*k52*k64 + k14*k32*k43*k54*k64 + k14*k32*k45*k52*k64 + k14*k34*k45*k52*k62 + k12*k34*k46*k54*k62 + k14*k32*k46*k54*k62 + k14*k34*k46*k52*k62 + k21*k32*k41*k52*k62 + k14*k34*k45*k52*k64 + k14*k34*k46*k54*k62 + k21*k32*k41*k52*k64 + k21*k32*k41*k54*k62 + k21*k32*k43*k52*k62 + k21*k34*k41*k52*k62 + k21*k32*k41*k54*k64 + k21*k32*k43*k52*k64 + k21*k32*k43*k54*k62 + k21*k32*k45*k52*k62 + k21*k34*k41*k52*k64 + k21*k34*k41*k54*k62 + k23*k34*k41*k52*k62 + k21*k32*k46*k52*k62 + k21*k32*k43*k54*k64 + k21*k32*k45*k52*k64 + k21*k34*k41*k54*k64 + k21*k34*k45*k52*k62 + k23*k34*k41*k52*k64 + k23*k34*k41*k54*k62 + k25*k32*k41*k54*k62 + k21*k32*k46*k54*k62 + k21*k34*k46*k52*k62 + k26*k32*k41*k52*k64 + k21*k34*k45*k52*k64 + k23*k34*k41*k54*k64 + k25*k32*k41*k54*k64 + k25*k34*k41*k54*k62 + k21*k34*k46*k54*k62 + k26*k32*k41*k54*k64 + k26*k34*k41*k52*k64 + k25*k34*k41*k54*k64 + k26*k34*k41*k54*k64);
E4=(k14*k21*k32*k52*k62 + k12*k23*k34*k52*k62 + k14*k21*k32*k52*k64 + k14*k21*k32*k54*k62 + k14*k21*k34*k52*k62 + k12*k23*k34*k52*k64 + k12*k23*k34*k54*k62 + k12*k25*k32*k54*k62 + k14*k21*k32*k54*k64 + k14*k21*k34*k52*k64 + k14*k21*k34*k54*k62 + k14*k23*k34*k52*k62 + k12*k26*k32*k52*k64 + k12*k23*k34*k54*k64 + k12*k25*k32*k54*k64 + k12*k25*k34*k54*k62 + k14*k21*k34*k54*k64 + k14*k23*k34*k52*k64 + k14*k23*k34*k54*k62 + k14*k25*k32*k54*k62 + k12*k26*k32*k54*k64 + k12*k26*k34*k52*k64 + k14*k26*k32*k52*k64 + k12*k25*k34*k54*k64 + k14*k23*k34*k54*k64 + k14*k25*k32*k54*k64 + k14*k25*k34*k54*k62 + k12*k26*k34*k54*k64 + k14*k26*k32*k54*k64 + k14*k26*k34*k52*k64 + k14*k25*k34*k54*k64 + k14*k26*k34*k54*k64)/(k12*k26*k32*k41*k52 + k12*k26*k32*k41*k54 + k12*k26*k32*k43*k52 + k12*k26*k34*k41*k52 + k14*k21*k32*k46*k52 + k12*k23*k34*k46*k52 + k12*k26*k32*k43*k54 + k12*k26*k32*k45*k52 + k12*k26*k34*k41*k54 + k14*k21*k32*k46*k54 + k14*k21*k34*k46*k52 + k14*k26*k32*k43*k52 + k12*k26*k32*k46*k52 + k12*k23*k34*k46*k54 + k12*k25*k32*k46*k54 + k12*k26*k34*k45*k52 + k14*k21*k34*k46*k54 + k14*k23*k34*k46*k52 + k14*k26*k32*k43*k54 + k14*k26*k32*k45*k52 + k12*k26*k32*k46*k54 + k12*k26*k34*k46*k52 + k14*k26*k32*k46*k52 + k12*k25*k34*k46*k54 + k14*k23*k34*k46*k54 + k14*k25*k32*k46*k54 + k14*k26*k34*k45*k52 + k12*k25*k32*k41*k62 + k12*k26*k34*k46*k54 + k14*k26*k32*k46*k54 + k14*k26*k34*k46*k52 + k14*k25*k34*k46*k54 + k12*k25*k32*k41*k64 + k12*k25*k32*k43*k62 + k12*k25*k34*k41*k62 + k14*k21*k32*k45*k62 + k14*k26*k34*k46*k54 + k12*k23*k34*k45*k62 + k12*k25*k32*k43*k64 + k12*k25*k32*k45*k62 + k12*k25*k34*k41*k64 + k14*k21*k32*k45*k64 + k14*k21*k34*k45*k62 + k14*k25*k32*k43*k62 + k12*k25*k32*k46*k62 + k12*k23*k34*k45*k64 + k12*k25*k32*k45*k64 + k12*k25*k34*k45*k62 + k14*k21*k34*k45*k64 + k14*k23*k34*k45*k62 + k14*k25*k32*k43*k64 + k14*k25*k32*k45*k62 + k12*k25*k34*k46*k62 + k12*k26*k32*k45*k64 + k14*k25*k32*k46*k62 + k12*k25*k34*k45*k64 + k14*k23*k34*k45*k64 + k14*k25*k32*k45*k64 + k14*k25*k34*k45*k62 + k12*k26*k34*k45*k64 + k14*k21*k32*k52*k62 + k14*k25*k34*k46*k62 + k14*k26*k32*k45*k64 + k14*k25*k34*k45*k64 + k12*k23*k34*k52*k62 + k14*k21*k32*k52*k64 + k14*k21*k32*k54*k62 + k14*k21*k34*k52*k62 + k14*k26*k34*k45*k64 + k12*k23*k34*k52*k64 + k12*k23*k34*k54*k62 + k12*k25*k32*k54*k62 + k14*k21*k32*k54*k64 + k14*k21*k34*k52*k64 + k14*k21*k34*k54*k62 + k14*k23*k34*k52*k62 + k12*k26*k32*k52*k64 + k12*k23*k34*k54*k64 + k12*k25*k32*k54*k64 + k12*k25*k34*k54*k62 + k14*k21*k34*k54*k64 + k14*k23*k34*k52*k64 + k14*k23*k34*k54*k62 + k14*k25*k32*k54*k62 + k12*k26*k32*k54*k64 + k12*k26*k34*k52*k64 + k14*k26*k32*k52*k64 + k12*k25*k34*k54*k64 + k14*k23*k34*k54*k64 + k14*k25*k32*k54*k64 + k14*k25*k34*k54*k62 + k12*k23*k41*k52*k62 + k12*k26*k34*k54*k64 + k14*k26*k32*k54*k64 + k14*k26*k34*k52*k64 + k14*k25*k34*k54*k64 + k12*k23*k41*k52*k64 + k12*k23*k41*k54*k62 + k12*k23*k43*k52*k62 + k14*k21*k43*k52*k62 + k14*k26*k34*k54*k64 + k12*k23*k41*k54*k64 + k12*k23*k43*k52*k64 + k12*k23*k43*k54*k62 + k12*k23*k45*k52*k62 + k14*k21*k43*k52*k64 + k14*k21*k43*k54*k62 + k14*k23*k43*k52*k62 + k12*k23*k46*k52*k62 + k12*k23*k43*k54*k64 + k12*k23*k45*k52*k64 + k12*k25*k43*k54*k62 + k14*k21*k43*k54*k64 + k14*k23*k43*k52*k64 + k14*k23*k43*k54*k62 + k14*k23*k45*k52*k62 + k12*k23*k46*k54*k62 + k12*k26*k43*k52*k64 + k14*k23*k46*k52*k62 + k12*k25*k43*k54*k64 + k14*k23*k43*k54*k64 + k14*k23*k45*k52*k64 + k14*k25*k43*k54*k62 + k12*k26*k43*k54*k64 + k12*k32*k41*k52*k62 + k14*k23*k46*k54*k62 + k14*k26*k43*k52*k64 + k14*k25*k43*k54*k64 + k12*k32*k41*k52*k64 + k12*k32*k41*k54*k62 + k12*k32*k43*k52*k62 + k12*k34*k41*k52*k62 + k14*k26*k43*k54*k64 + k12*k32*k41*k54*k64 + k12*k32*k43*k52*k64 + k12*k32*k43*k54*k62 + k12*k32*k45*k52*k62 + k12*k34*k41*k52*k64 + k12*k34*k41*k54*k62 + k14*k32*k43*k52*k62 + k12*k32*k46*k52*k62 + k12*k32*k43*k54*k64 + k12*k32*k45*k52*k64 + k12*k34*k41*k54*k64 + k12*k34*k45*k52*k62 + k14*k32*k43*k52*k64 + k14*k32*k43*k54*k62 + k14*k32*k45*k52*k62 + k12*k32*k46*k54*k62 + k12*k34*k46*k52*k62 + k14*k32*k46*k52*k62 + k12*k34*k45*k52*k64 + k14*k32*k43*k54*k64 + k14*k32*k45*k52*k64 + k14*k34*k45*k52*k62 + k12*k34*k46*k54*k62 + k14*k32*k46*k54*k62 + k14*k34*k46*k52*k62 + k21*k32*k41*k52*k62 + k14*k34*k45*k52*k64 + k14*k34*k46*k54*k62 + k21*k32*k41*k52*k64 + k21*k32*k41*k54*k62 + k21*k32*k43*k52*k62 + k21*k34*k41*k52*k62 + k21*k32*k41*k54*k64 + k21*k32*k43*k52*k64 + k21*k32*k43*k54*k62 + k21*k32*k45*k52*k62 + k21*k34*k41*k52*k64 + k21*k34*k41*k54*k62 + k23*k34*k41*k52*k62 + k21*k32*k46*k52*k62 + k21*k32*k43*k54*k64 + k21*k32*k45*k52*k64 + k21*k34*k41*k54*k64 + k21*k34*k45*k52*k62 + k23*k34*k41*k52*k64 + k23*k34*k41*k54*k62 + k25*k32*k41*k54*k62 + k21*k32*k46*k54*k62 + k21*k34*k46*k52*k62 + k26*k32*k41*k52*k64 + k21*k34*k45*k52*k64 + k23*k34*k41*k54*k64 + k25*k32*k41*k54*k64 + k25*k34*k41*k54*k62 + k21*k34*k46*k54*k62 + k26*k32*k41*k54*k64 + k26*k34*k41*k52*k64 + k25*k34*k41*k54*k64 + k26*k34*k41*k54*k64);
E5=(k12*k25*k32*k41*k62 + k12*k25*k32*k41*k64 + k12*k25*k32*k43*k62 + k12*k25*k34*k41*k62 + k14*k21*k32*k45*k62 + k12*k23*k34*k45*k62 + k12*k25*k32*k43*k64 + k12*k25*k32*k45*k62 + k12*k25*k34*k41*k64 + k14*k21*k32*k45*k64 + k14*k21*k34*k45*k62 + k14*k25*k32*k43*k62 + k12*k25*k32*k46*k62 + k12*k23*k34*k45*k64 + k12*k25*k32*k45*k64 + k12*k25*k34*k45*k62 + k14*k21*k34*k45*k64 + k14*k23*k34*k45*k62 + k14*k25*k32*k43*k64 + k14*k25*k32*k45*k62 + k12*k25*k34*k46*k62 + k12*k26*k32*k45*k64 + k14*k25*k32*k46*k62 + k12*k25*k34*k45*k64 + k14*k23*k34*k45*k64 + k14*k25*k32*k45*k64 + k14*k25*k34*k45*k62 + k12*k26*k34*k45*k64 + k14*k25*k34*k46*k62 + k14*k26*k32*k45*k64 + k14*k25*k34*k45*k64 + k14*k26*k34*k45*k64)/(k12*k26*k32*k41*k52 + k12*k26*k32*k41*k54 + k12*k26*k32*k43*k52 + k12*k26*k34*k41*k52 + k14*k21*k32*k46*k52 + k12*k23*k34*k46*k52 + k12*k26*k32*k43*k54 + k12*k26*k32*k45*k52 + k12*k26*k34*k41*k54 + k14*k21*k32*k46*k54 + k14*k21*k34*k46*k52 + k14*k26*k32*k43*k52 + k12*k26*k32*k46*k52 + k12*k23*k34*k46*k54 + k12*k25*k32*k46*k54 + k12*k26*k34*k45*k52 + k14*k21*k34*k46*k54 + k14*k23*k34*k46*k52 + k14*k26*k32*k43*k54 + k14*k26*k32*k45*k52 + k12*k26*k32*k46*k54 + k12*k26*k34*k46*k52 + k14*k26*k32*k46*k52 + k12*k25*k34*k46*k54 + k14*k23*k34*k46*k54 + k14*k25*k32*k46*k54 + k14*k26*k34*k45*k52 + k12*k25*k32*k41*k62 + k12*k26*k34*k46*k54 + k14*k26*k32*k46*k54 + k14*k26*k34*k46*k52 + k14*k25*k34*k46*k54 + k12*k25*k32*k41*k64 + k12*k25*k32*k43*k62 + k12*k25*k34*k41*k62 + k14*k21*k32*k45*k62 + k14*k26*k34*k46*k54 + k12*k23*k34*k45*k62 + k12*k25*k32*k43*k64 + k12*k25*k32*k45*k62 + k12*k25*k34*k41*k64 + k14*k21*k32*k45*k64 + k14*k21*k34*k45*k62 + k14*k25*k32*k43*k62 + k12*k25*k32*k46*k62 + k12*k23*k34*k45*k64 + k12*k25*k32*k45*k64 + k12*k25*k34*k45*k62 + k14*k21*k34*k45*k64 + k14*k23*k34*k45*k62 + k14*k25*k32*k43*k64 + k14*k25*k32*k45*k62 + k12*k25*k34*k46*k62 + k12*k26*k32*k45*k64 + k14*k25*k32*k46*k62 + k12*k25*k34*k45*k64 + k14*k23*k34*k45*k64 + k14*k25*k32*k45*k64 + k14*k25*k34*k45*k62 + k12*k26*k34*k45*k64 + k14*k21*k32*k52*k62 + k14*k25*k34*k46*k62 + k14*k26*k32*k45*k64 + k14*k25*k34*k45*k64 + k12*k23*k34*k52*k62 + k14*k21*k32*k52*k64 + k14*k21*k32*k54*k62 + k14*k21*k34*k52*k62 + k14*k26*k34*k45*k64 + k12*k23*k34*k52*k64 + k12*k23*k34*k54*k62 + k12*k25*k32*k54*k62 + k14*k21*k32*k54*k64 + k14*k21*k34*k52*k64 + k14*k21*k34*k54*k62 + k14*k23*k34*k52*k62 + k12*k26*k32*k52*k64 + k12*k23*k34*k54*k64 + k12*k25*k32*k54*k64 + k12*k25*k34*k54*k62 + k14*k21*k34*k54*k64 + k14*k23*k34*k52*k64 + k14*k23*k34*k54*k62 + k14*k25*k32*k54*k62 + k12*k26*k32*k54*k64 + k12*k26*k34*k52*k64 + k14*k26*k32*k52*k64 + k12*k25*k34*k54*k64 + k14*k23*k34*k54*k64 + k14*k25*k32*k54*k64 + k14*k25*k34*k54*k62 + k12*k23*k41*k52*k62 + k12*k26*k34*k54*k64 + k14*k26*k32*k54*k64 + k14*k26*k34*k52*k64 + k14*k25*k34*k54*k64 + k12*k23*k41*k52*k64 + k12*k23*k41*k54*k62 + k12*k23*k43*k52*k62 + k14*k21*k43*k52*k62 + k14*k26*k34*k54*k64 + k12*k23*k41*k54*k64 + k12*k23*k43*k52*k64 + k12*k23*k43*k54*k62 + k12*k23*k45*k52*k62 + k14*k21*k43*k52*k64 + k14*k21*k43*k54*k62 + k14*k23*k43*k52*k62 + k12*k23*k46*k52*k62 + k12*k23*k43*k54*k64 + k12*k23*k45*k52*k64 + k12*k25*k43*k54*k62 + k14*k21*k43*k54*k64 + k14*k23*k43*k52*k64 + k14*k23*k43*k54*k62 + k14*k23*k45*k52*k62 + k12*k23*k46*k54*k62 + k12*k26*k43*k52*k64 + k14*k23*k46*k52*k62 + k12*k25*k43*k54*k64 + k14*k23*k43*k54*k64 + k14*k23*k45*k52*k64 + k14*k25*k43*k54*k62 + k12*k26*k43*k54*k64 + k12*k32*k41*k52*k62 + k14*k23*k46*k54*k62 + k14*k26*k43*k52*k64 + k14*k25*k43*k54*k64 + k12*k32*k41*k52*k64 + k12*k32*k41*k54*k62 + k12*k32*k43*k52*k62 + k12*k34*k41*k52*k62 + k14*k26*k43*k54*k64 + k12*k32*k41*k54*k64 + k12*k32*k43*k52*k64 + k12*k32*k43*k54*k62 + k12*k32*k45*k52*k62 + k12*k34*k41*k52*k64 + k12*k34*k41*k54*k62 + k14*k32*k43*k52*k62 + k12*k32*k46*k52*k62 + k12*k32*k43*k54*k64 + k12*k32*k45*k52*k64 + k12*k34*k41*k54*k64 + k12*k34*k45*k52*k62 + k14*k32*k43*k52*k64 + k14*k32*k43*k54*k62 + k14*k32*k45*k52*k62 + k12*k32*k46*k54*k62 + k12*k34*k46*k52*k62 + k14*k32*k46*k52*k62 + k12*k34*k45*k52*k64 + k14*k32*k43*k54*k64 + k14*k32*k45*k52*k64 + k14*k34*k45*k52*k62 + k12*k34*k46*k54*k62 + k14*k32*k46*k54*k62 + k14*k34*k46*k52*k62 + k21*k32*k41*k52*k62 + k14*k34*k45*k52*k64 + k14*k34*k46*k54*k62 + k21*k32*k41*k52*k64 + k21*k32*k41*k54*k62 + k21*k32*k43*k52*k62 + k21*k34*k41*k52*k62 + k21*k32*k41*k54*k64 + k21*k32*k43*k52*k64 + k21*k32*k43*k54*k62 + k21*k32*k45*k52*k62 + k21*k34*k41*k52*k64 + k21*k34*k41*k54*k62 + k23*k34*k41*k52*k62 + k21*k32*k46*k52*k62 + k21*k32*k43*k54*k64 + k21*k32*k45*k52*k64 + k21*k34*k41*k54*k64 + k21*k34*k45*k52*k62 + k23*k34*k41*k52*k64 + k23*k34*k41*k54*k62 + k25*k32*k41*k54*k62 + k21*k32*k46*k54*k62 + k21*k34*k46*k52*k62 + k26*k32*k41*k52*k64 + k21*k34*k45*k52*k64 + k23*k34*k41*k54*k64 + k25*k32*k41*k54*k64 + k25*k34*k41*k54*k62 + k21*k34*k46*k54*k62 + k26*k32*k41*k54*k64 + k26*k34*k41*k52*k64 + k25*k34*k41*k54*k64 + k26*k34*k41*k54*k64);
% E6=(k12*k26*k32*k41*k52 + k12*k26*k32*k41*k54 + k12*k26*k32*k43*k52 + k12*k26*k34*k41*k52 + k14*k21*k32*k46*k52 + k12*k23*k34*k46*k52 + k12*k26*k32*k43*k54 + k12*k26*k32*k45*k52 + k12*k26*k34*k41*k54 + k14*k21*k32*k46*k54 + k14*k21*k34*k46*k52 + k14*k26*k32*k43*k52 + k12*k26*k32*k46*k52 + k12*k23*k34*k46*k54 + k12*k25*k32*k46*k54 + k12*k26*k34*k45*k52 + k14*k21*k34*k46*k54 + k14*k23*k34*k46*k52 + k14*k26*k32*k43*k54 + k14*k26*k32*k45*k52 + k12*k26*k32*k46*k54 + k12*k26*k34*k46*k52 + k14*k26*k32*k46*k52 + k12*k25*k34*k46*k54 + k14*k23*k34*k46*k54 + k14*k25*k32*k46*k54 + k14*k26*k34*k45*k52 + k12*k26*k34*k46*k54 + k14*k26*k32*k46*k54 + k14*k26*k34*k46*k52 + k14*k25*k34*k46*k54 + k14*k26*k34*k46*k54)/(k12*k26*k32*k41*k52 + k12*k26*k32*k41*k54 + k12*k26*k32*k43*k52 + k12*k26*k34*k41*k52 + k14*k21*k32*k46*k52 + k12*k23*k34*k46*k52 + k12*k26*k32*k43*k54 + k12*k26*k32*k45*k52 + k12*k26*k34*k41*k54 + k14*k21*k32*k46*k54 + k14*k21*k34*k46*k52 + k14*k26*k32*k43*k52 + k12*k26*k32*k46*k52 + k12*k23*k34*k46*k54 + k12*k25*k32*k46*k54 + k12*k26*k34*k45*k52 + k14*k21*k34*k46*k54 + k14*k23*k34*k46*k52 + k14*k26*k32*k43*k54 + k14*k26*k32*k45*k52 + k12*k26*k32*k46*k54 + k12*k26*k34*k46*k52 + k14*k26*k32*k46*k52 + k12*k25*k34*k46*k54 + k14*k23*k34*k46*k54 + k14*k25*k32*k46*k54 + k14*k26*k34*k45*k52 + k12*k25*k32*k41*k62 + k12*k26*k34*k46*k54 + k14*k26*k32*k46*k54 + k14*k26*k34*k46*k52 + k14*k25*k34*k46*k54 + k12*k25*k32*k41*k64 + k12*k25*k32*k43*k62 + k12*k25*k34*k41*k62 + k14*k21*k32*k45*k62 + k14*k26*k34*k46*k54 + k12*k23*k34*k45*k62 + k12*k25*k32*k43*k64 + k12*k25*k32*k45*k62 + k12*k25*k34*k41*k64 + k14*k21*k32*k45*k64 + k14*k21*k34*k45*k62 + k14*k25*k32*k43*k62 + k12*k25*k32*k46*k62 + k12*k23*k34*k45*k64 + k12*k25*k32*k45*k64 + k12*k25*k34*k45*k62 + k14*k21*k34*k45*k64 + k14*k23*k34*k45*k62 + k14*k25*k32*k43*k64 + k14*k25*k32*k45*k62 + k12*k25*k34*k46*k62 + k12*k26*k32*k45*k64 + k14*k25*k32*k46*k62 + k12*k25*k34*k45*k64 + k14*k23*k34*k45*k64 + k14*k25*k32*k45*k64 + k14*k25*k34*k45*k62 + k12*k26*k34*k45*k64 + k14*k21*k32*k52*k62 + k14*k25*k34*k46*k62 + k14*k26*k32*k45*k64 + k14*k25*k34*k45*k64 + k12*k23*k34*k52*k62 + k14*k21*k32*k52*k64 + k14*k21*k32*k54*k62 + k14*k21*k34*k52*k62 + k14*k26*k34*k45*k64 + k12*k23*k34*k52*k64 + k12*k23*k34*k54*k62 + k12*k25*k32*k54*k62 + k14*k21*k32*k54*k64 + k14*k21*k34*k52*k64 + k14*k21*k34*k54*k62 + k14*k23*k34*k52*k62 + k12*k26*k32*k52*k64 + k12*k23*k34*k54*k64 + k12*k25*k32*k54*k64 + k12*k25*k34*k54*k62 + k14*k21*k34*k54*k64 + k14*k23*k34*k52*k64 + k14*k23*k34*k54*k62 + k14*k25*k32*k54*k62 + k12*k26*k32*k54*k64 + k12*k26*k34*k52*k64 + k14*k26*k32*k52*k64 + k12*k25*k34*k54*k64 + k14*k23*k34*k54*k64 + k14*k25*k32*k54*k64 + k14*k25*k34*k54*k62 + k12*k23*k41*k52*k62 + k12*k26*k34*k54*k64 + k14*k26*k32*k54*k64 + k14*k26*k34*k52*k64 + k14*k25*k34*k54*k64 + k12*k23*k41*k52*k64 + k12*k23*k41*k54*k62 + k12*k23*k43*k52*k62 + k14*k21*k43*k52*k62 + k14*k26*k34*k54*k64 + k12*k23*k41*k54*k64 + k12*k23*k43*k52*k64 + k12*k23*k43*k54*k62 + k12*k23*k45*k52*k62 + k14*k21*k43*k52*k64 + k14*k21*k43*k54*k62 + k14*k23*k43*k52*k62 + k12*k23*k46*k52*k62 + k12*k23*k43*k54*k64 + k12*k23*k45*k52*k64 + k12*k25*k43*k54*k62 + k14*k21*k43*k54*k64 + k14*k23*k43*k52*k64 + k14*k23*k43*k54*k62 + k14*k23*k45*k52*k62 + k12*k23*k46*k54*k62 + k12*k26*k43*k52*k64 + k14*k23*k46*k52*k62 + k12*k25*k43*k54*k64 + k14*k23*k43*k54*k64 + k14*k23*k45*k52*k64 + k14*k25*k43*k54*k62 + k12*k26*k43*k54*k64 + k12*k32*k41*k52*k62 + k14*k23*k46*k54*k62 + k14*k26*k43*k52*k64 + k14*k25*k43*k54*k64 + k12*k32*k41*k52*k64 + k12*k32*k41*k54*k62 + k12*k32*k43*k52*k62 + k12*k34*k41*k52*k62 + k14*k26*k43*k54*k64 + k12*k32*k41*k54*k64 + k12*k32*k43*k52*k64 + k12*k32*k43*k54*k62 + k12*k32*k45*k52*k62 + k12*k34*k41*k52*k64 + k12*k34*k41*k54*k62 + k14*k32*k43*k52*k62 + k12*k32*k46*k52*k62 + k12*k32*k43*k54*k64 + k12*k32*k45*k52*k64 + k12*k34*k41*k54*k64 + k12*k34*k45*k52*k62 + k14*k32*k43*k52*k64 + k14*k32*k43*k54*k62 + k14*k32*k45*k52*k62 + k12*k32*k46*k54*k62 + k12*k34*k46*k52*k62 + k14*k32*k46*k52*k62 + k12*k34*k45*k52*k64 + k14*k32*k43*k54*k64 + k14*k32*k45*k52*k64 + k14*k34*k45*k52*k62 + k12*k34*k46*k54*k62 + k14*k32*k46*k54*k62 + k14*k34*k46*k52*k62 + k21*k32*k41*k52*k62 + k14*k34*k45*k52*k64 + k14*k34*k46*k54*k62 + k21*k32*k41*k52*k64 + k21*k32*k41*k54*k62 + k21*k32*k43*k52*k62 + k21*k34*k41*k52*k62 + k21*k32*k41*k54*k64 + k21*k32*k43*k52*k64 + k21*k32*k43*k54*k62 + k21*k32*k45*k52*k62 + k21*k34*k41*k52*k64 + k21*k34*k41*k54*k62 + k23*k34*k41*k52*k62 + k21*k32*k46*k52*k62 + k21*k32*k43*k54*k64 + k21*k32*k45*k52*k64 + k21*k34*k41*k54*k64 + k21*k34*k45*k52*k62 + k23*k34*k41*k52*k64 + k23*k34*k41*k54*k62 + k25*k32*k41*k54*k62 + k21*k32*k46*k54*k62 + k21*k34*k46*k52*k62 + k26*k32*k41*k52*k64 + k21*k34*k45*k52*k64 + k23*k34*k41*k54*k64 + k25*k32*k41*k54*k64 + k25*k34*k41*k54*k62 + k21*k34*k46*k54*k62 + k26*k32*k41*k54*k64 + k26*k34*k41*k52*k64 + k25*k34*k41*k54*k64 + k26*k34*k41*k54*k64);
E6 = 1 - E1 - E2 - E3 - E4 - E5;
% Net Turnover Rates
JQH2 = ETC3_activity*((kfQH2c3_12*E1 - krQH2c3_12*E2) + (kfQH2c3_34*E3 - krQH2c3_34*E4) + (kfQH2c3_45*E4 - krQH2c3_45*E5) + (kfQH2c3_26*E2 - krQH2c3_26*E6) - (kfQH2_41*E4 - krQH2_41*E1) - (kfQH2_52*E5 - krQH2_52*E2));
% QH2 consumption
JSO = ETC3_activity*((kfSO_21*E2 - krSO_21*E1) + (kfSO_43*E4 - krSO_43*E3) + (kfSO_54*E5 - krSO_54*E4) + (kfSO_62*E6 - krSO_62*E2));
% superoxide production
% for numerical stability, enforce strict stoichiometric coupling
J_ETC3_im_to_matrix = (2*JQH2 - JSO);
%ETC4:im_to_matrix
ETC4_activity = par(8);
c3 = cytocox_im; % M
c2 = cytocred_im; % M
O2 = 200e-6; % constant oxygen concentration (M)
% Hx = h_im; % M
% Hi = h_matrix; % M
dPsi = DPsi_im_to_matrix; % mV
KM = 1.6115e-04; % M
beta = 6.6054e-06; % unitless
DGro_ETC4 =-202.524;
Keq_ETC4 = exp(-DGro_ETC4/RT)* exp( -(4*F*DPsi_im_to_matrix)/RT)*P(14)^2/P(15)^2*h_matrix^4/h_im^2;
J_ETC4_im_to_matrix = 1e-0*ETC4_activity*(O2./(O2+1E-6)).*(c2.^2./(c2.^2 + KM^2*(1+beta*exp(dPsi*2*F/RT)))).*(1 - (c3.^2./c2.^2./O2.^.5)./Keq_ETC4);
% if O2 < 1e-12 | c2 < 1e-9
% J_ETC4_im_to_matrix = 0;
% end
%HLEAK:im_to_matrix
x_HLE =par(9);
FRT = DPsi_im_to_matrix*F/RT/2;
J_HLEAK_im_to_matrix = x_HLE*(h_im*exp(FRT) - h_matrix*exp(-FRT));
%PIH:im_to_matrix
a = Pi_im*(h_im/Kh(16))/P(16);
p = Pi_matrix*(h_matrix/Kh(5))/P(5);
x_PIH =par(10);
k_PIH = 1.61e-3;
J_PIH_im_to_matrix = (x_PIH/k_PIH)*(h_im*a - h_matrix*p)/(1+a/k_PIH)/(1+k_PIH);
%ANT:im_to_matrix
x_ANT =par(11);
ADP_i1 = ADP_im/P(17);
% ADP^3-;
ATP_i1 = ATP_im/P(18);
% ATP^4-;
ADP_x1 = ADP_matrix/P(3);
% ADP^3-;
ATP_x1 = ATP_matrix/P(4);
% ATP^4-;
del_D = 0.0167;
del_T = 0.0699;
k2_ANT = 9.54/60;
% = 1.59e-1
k3_ANT = 30.05/60;
% = 5.01e-1
K_D_o_ANT = 38.89e-6;
K_T_o_ANT = 56.05e-6;
A = +0.2829;
B = -0.2086;
C = +0.2372;
fi = F*DPsi_im_to_matrix/RT;
k2_ANT_fi = k2_ANT*exp((A*(-3)+B*(-4)+C)*fi);
k3_ANT_fi = k3_ANT*exp((A*(-4)+B*(-3)+C)*fi);
K_D_o_ANT_fi = K_D_o_ANT*exp(3*del_D*fi);
K_T_o_ANT_fi = K_T_o_ANT*exp(4*del_T*fi);
q = k3_ANT_fi*K_D_o_ANT_fi*exp(fi)/(k2_ANT_fi*K_T_o_ANT_fi);
term2 = k2_ANT_fi*ATP_x1.*ADP_i1*q/K_D_o_ANT_fi ;
term3 = k3_ANT_fi.*ADP_x1.*ATP_i1/K_T_o_ANT_fi;
num = term2 - term3;
den = (1 + ATP_i1/K_T_o_ANT_fi + ADP_i1/K_D_o_ANT_fi)*(ADP_x1 + ATP_x1*q);
J_ANT_im_to_matrix = x_ANT*num/den;
% x_ANT'(in the paper) = x_ANT/7.2679e-003*(0.70e-1);
%KH:im_to_matrix
k1_KH =par(12);
J_KH_im_to_matrix = k1_KH*(k_im*h_matrix - k_matrix*h_im);
%AMPPERM:cytoplasm_to_im
x_AMPPERM =par(13);
J_AMPPERM_cytoplasm_to_im = gamma * x_AMPPERM * (AMP_c - AMP_im);
%ADPPERM:cytoplasm_to_im
x_ADPPERM =par(14);
J_ADPPERM_cytoplasm_to_im = gamma * x_ADPPERM * (ADP_c - ADP_im);
%ATPPERM:cytoplasm_to_im
x_ATPPERM =par(15);
J_ATPPERM_cytoplasm_to_im = gamma * x_ATPPERM * (ATP_c - ATP_im);
%PIPERM:cytoplasm_to_im
x_PIPERM =par(16);
J_PIPERM_cytoplasm_to_im = gamma * x_PIPERM * (Pi_c - Pi_im);
%HPERM:cytoplasm_to_im
x_HPERM =par(17);
J_HPERM_cytoplasm_to_im = x_HPERM * (h_c - h_im);
%KPERM:cytoplasm_to_im
x_KPERM =par(18);
J_KPERM_cytoplasm_to_im = x_KPERM * (k_c - k_im);
%MPERM:cytoplasm_to_im
x_MPERM =par(19);
J_MPERM_cytoplasm_to_im = x_MPERM * (m_c - m_im);

%% REACTANT TIME DERIVATIVES
f(1,:) = ( 0  - 1*J_DH_matrix + 1*J_ETC1_im_to_matrix ) / VWater_matrix; % NAD_matrix
f(2,:) = ( 0  + 1*J_DH_matrix - 1*J_ETC1_im_to_matrix ) / VWater_matrix; % NADH_matrix
f(3,:) = ( 0  - 1*J_F1F0ATPASE_im_to_matrix + 1*J_ANT_im_to_matrix ) / VWater_matrix; % ADP_matrix
f(4,:) = ( 0  + 1*J_F1F0ATPASE_im_to_matrix - 1*J_ANT_im_to_matrix ) / VWater_matrix; % ATP_matrix
f(5,:) = ( 0  - 1*J_F1F0ATPASE_im_to_matrix + 1*J_PIH_im_to_matrix ) / VWater_matrix; % Pi_matrix
f(6,:) = ( 0  - 1*J_SDH_BBV_matrix - 1*J_ETC1_im_to_matrix + 1*J_ETC3_im_to_matrix ) / VWater_matrix; % coQ_matrix
f(7,:) = ( 0  + 1*J_SDH_BBV_matrix + 1*J_ETC1_im_to_matrix - 1*J_ETC3_im_to_matrix ) / VWater_matrix; % coQH2_matrix
f(8,:) = ( 0  - 1*J_ATPASE_cytoplasm + 1*J_CK_cytoplasm + 1*J_AK_cytoplasm - 1*J_ATPPERM_cytoplasm_to_im/VRegion_cytoplasm*VRegion_im ) / VWater_cytoplasm; % ATP_c
f(9,:) = ( 0  + 1*J_ATPASE_cytoplasm - 1*J_CK_cytoplasm - 2*J_AK_cytoplasm - 1*J_ADPPERM_cytoplasm_to_im/VRegion_cytoplasm*VRegion_im ) / VWater_cytoplasm; % ADP_c
f(10,:) = ( 0  + 1*J_ATPASE_cytoplasm - 1*J_PIPERM_cytoplasm_to_im/VRegion_cytoplasm*VRegion_im ) / VWater_cytoplasm; % Pi_c
f(11,:) = ( 0  - 1*J_CK_cytoplasm ) / VWater_cytoplasm; % phosphocreatine_c
f(12,:) = ( 0  + 1*J_CK_cytoplasm ) / VWater_cytoplasm; % creatine_c
f(13,:) = ( 0  + 1*J_AK_cytoplasm - 1*J_AMPPERM_cytoplasm_to_im/VRegion_cytoplasm*VRegion_im ) / VWater_cytoplasm; % AMP_c
f(14,:) = ( 0  - 2*J_ETC3_im_to_matrix + 2*J_ETC4_im_to_matrix ) / VWater_im; % cytocox_im
f(15,:) = ( 0  + 2*J_ETC3_im_to_matrix - 2*J_ETC4_im_to_matrix ) / VWater_im; % cytocred_im
f(16,:) = ( 0  - 1*J_PIH_im_to_matrix + 1*J_PIPERM_cytoplasm_to_im ) / VWater_im; % Pi_im
f(17,:) = ( 0  - 1*J_ANT_im_to_matrix + 1*J_ADPPERM_cytoplasm_to_im ) / VWater_im; % ADP_im
f(18,:) = ( 0  + 1*J_ANT_im_to_matrix + 1*J_ATPPERM_cytoplasm_to_im ) / VWater_im; % ATP_im
f(19,:) = ( 0  + 1*J_AMPPERM_cytoplasm_to_im ) / VWater_im; % AMP_im

%% ION EQUATIONS
% COMPARTMENT matrix:
ii = [1   2   3   4   5   6   7]; % Indices of SVs in compartment matrix
% PARTIAL DERIVATIVES
pHBpK = -sum( (h_matrix*x(ii)'./Kh(ii))./(Kk(ii).*P(ii).^2) );
pHBpM = -sum( (h_matrix*x(ii)'./Kh(ii))./(Km(ii).*P(ii).^2) );
pHBpH = +sum( (1+m_matrix./Km(ii)+k_matrix./Kk(ii)).*x(ii)'./(Kh(ii).*P(ii).^2) );
pMBpH = -sum( (m_matrix*x(ii)'./Km(ii))./(Kh(ii).*P(ii).^2) );
pMBpK = -sum( (m_matrix*x(ii)'./Km(ii))./(Kk(ii).*P(ii).^2) );
pMBpM = +sum( (1+h_matrix./Kh(ii)+k_matrix./Kk(ii)).*x(ii)'./(Km(ii).*P(ii).^2) );
pKBpH = -sum( (k_matrix*x(ii)'./Kk(ii))./(Kh(ii).*P(ii).^2) );
pKBpM = -sum( (k_matrix*x(ii)'./Kk(ii))./(Km(ii).*P(ii).^2) );
pKBpK = +sum( (1+h_matrix./Kh(ii)+m_matrix./Km(ii)).*x(ii)'./(Kk(ii).*P(ii).^2) );
% PHIs
J_H = (0 + 1*J_DH_matrix + 0*J_SDH_BBV_matrix + 1.66667*J_F1F0ATPASE_im_to_matrix + 0*J_F1F0ATPASE_im_to_matrix - 5*J_ETC1_im_to_matrix - 2*J_ETC3_im_to_matrix - 4*J_ETC4_im_to_matrix + 1*J_HLEAK_im_to_matrix + 2*J_PIH_im_to_matrix - 1*J_KH_im_to_matrix) / VWater_matrix;
J_M = (0) / VWater_matrix;
J_K = (0 + 1*J_KH_im_to_matrix) / VWater_matrix;
Phi_H = J_H - sum( h_matrix*f(ii)'./(Kh(ii).*P(ii)) );
Phi_M = J_M - sum( m_matrix*f(ii)'./(Km(ii).*P(ii)) );
Phi_K = J_K - sum( k_matrix*f(ii)'./(Kk(ii).*P(ii)) );
% ALPHAs
% aH = 1 + pHBpH;
aM = 1 + pMBpM;
aK = 1 + pKBpK;
% ADDITIONAL BUFFER for [H+]
aH = 1 + pHBpH + BX(1)/K_BX(1)/(1+h_matrix/K_BX(1))^2; % M
% DENOMINATOR
D = aH*pKBpM*pMBpK + aK*pHBpM*pMBpH + aM*pHBpK*pKBpH - ...
    aM*aK*aH - pHBpK*pKBpM*pMBpH - pHBpM*pMBpK*pKBpH;
% DERIVATIVES for H,Mg,K
f(20,:) =  ( (pKBpM*pMBpK - aM*aK)*Phi_H + ...
            (aK*pHBpM - pHBpK*pKBpM)*Phi_M + ...
            (aM*pHBpK - pHBpM*pMBpK)*Phi_K ) / D;
f(21,:) =  ( (aK*pMBpH - pKBpH*pMBpK)*Phi_H + ...
            (pKBpH*pHBpK - aH*aK)*Phi_M + ...
            (aH*pMBpK - pHBpK*pMBpH)*Phi_K ) / D;
f(22,:) =  ( (aM*pKBpH - pKBpM*pMBpH)*Phi_H + ...
            (aH*pKBpM - pKBpH*pHBpM)*Phi_M + ...
            (pMBpH*pHBpM - aH*aM)*Phi_K ) / D;
% COMPARTMENT cytoplasm:
ii = [8   9  10  11  12  13]; % Indices of SVs in compartment cytoplasm
% PARTIAL DERIVATIVES
pHBpK = -sum( (h_c*x(ii)'./Kh(ii))./(Kk(ii).*P(ii).^2) );
pHBpM = -sum( (h_c*x(ii)'./Kh(ii))./(Km(ii).*P(ii).^2) );
pHBpH = +sum( (1+m_c./Km(ii)+k_c./Kk(ii)).*x(ii)'./(Kh(ii).*P(ii).^2) );
pMBpH = -sum( (m_c*x(ii)'./Km(ii))./(Kh(ii).*P(ii).^2) );
pMBpK = -sum( (m_c*x(ii)'./Km(ii))./(Kk(ii).*P(ii).^2) );
pMBpM = +sum( (1+h_c./Kh(ii)+k_c./Kk(ii)).*x(ii)'./(Km(ii).*P(ii).^2) );
pKBpH = -sum( (k_c*x(ii)'./Kk(ii))./(Kh(ii).*P(ii).^2) );
pKBpM = -sum( (k_c*x(ii)'./Kk(ii))./(Km(ii).*P(ii).^2) );
pKBpK = +sum( (1+h_c./Kh(ii)+m_c./Km(ii)).*x(ii)'./(Kk(ii).*P(ii).^2) );
% PHIs
J_H = (0 + 1*J_ATPASE_cytoplasm - 1*J_CK_cytoplasm + 0*J_AK_cytoplasm - 1*J_HPERM_cytoplasm_to_im/VRegion_cytoplasm*VRegion_im) / VWater_cytoplasm;
J_M = (0 - 1*J_MPERM_cytoplasm_to_im/VRegion_cytoplasm*VRegion_im) / VWater_cytoplasm;
J_K = (0 - 1*J_KPERM_cytoplasm_to_im/VRegion_cytoplasm*VRegion_im) / VWater_cytoplasm;
Phi_H = J_H - sum( h_c*f(ii)'./(Kh(ii).*P(ii)) );
Phi_M = J_M - sum( m_c*f(ii)'./(Km(ii).*P(ii)) );
Phi_K = J_K - sum( k_c*f(ii)'./(Kk(ii).*P(ii)) );
% ALPHAs
% aH = 1 + pHBpH;
aM = 1 + pMBpM;
aK = 1 + pKBpK;
% ADDITIONAL BUFFER for [H+]
aH = 1 + pHBpH + BX(2)/K_BX(2)/(1+h_c/K_BX(2))^2; % M
% DENOMINATOR
D = aH*pKBpM*pMBpK + aK*pHBpM*pMBpH + aM*pHBpK*pKBpH - ...
    aM*aK*aH - pHBpK*pKBpM*pMBpH - pHBpM*pMBpK*pKBpH;
% DERIVATIVES for H,Mg,K
f(23,:) =  ( (pKBpM*pMBpK - aM*aK)*Phi_H + ...
            (aK*pHBpM - pHBpK*pKBpM)*Phi_M + ...
            (aM*pHBpK - pHBpM*pMBpK)*Phi_K ) / D;
f(24,:) =  ( (aK*pMBpH - pKBpH*pMBpK)*Phi_H + ...
            (pKBpH*pHBpK - aH*aK)*Phi_M + ...
            (aH*pMBpK - pHBpK*pMBpH)*Phi_K ) / D;
f(25,:) =  ( (aM*pKBpH - pKBpM*pMBpH)*Phi_H + ...
            (aH*pKBpM - pKBpH*pHBpM)*Phi_M + ...
            (pMBpH*pHBpM - aH*aM)*Phi_K ) / D;
% COMPARTMENT im:
ii = [14  15  16  17  18  19]; % Indices of SVs in compartment im
% PARTIAL DERIVATIVES
pHBpK = -sum( (h_im*x(ii)'./Kh(ii))./(Kk(ii).*P(ii).^2) );
pHBpM = -sum( (h_im*x(ii)'./Kh(ii))./(Km(ii).*P(ii).^2) );
pHBpH = +sum( (1+m_im./Km(ii)+k_im./Kk(ii)).*x(ii)'./(Kh(ii).*P(ii).^2) );
pMBpH = -sum( (m_im*x(ii)'./Km(ii))./(Kh(ii).*P(ii).^2) );
pMBpK = -sum( (m_im*x(ii)'./Km(ii))./(Kk(ii).*P(ii).^2) );
pMBpM = +sum( (1+h_im./Kh(ii)+k_im./Kk(ii)).*x(ii)'./(Km(ii).*P(ii).^2) );
pKBpH = -sum( (k_im*x(ii)'./Kk(ii))./(Kh(ii).*P(ii).^2) );
pKBpM = -sum( (k_im*x(ii)'./Kk(ii))./(Km(ii).*P(ii).^2) );
pKBpK = +sum( (1+h_im./Kh(ii)+m_im./Km(ii)).*x(ii)'./(Kk(ii).*P(ii).^2) );
% PHIs
J_H = (0 - 2.66667*J_F1F0ATPASE_im_to_matrix + 4*J_ETC1_im_to_matrix + 4*J_ETC3_im_to_matrix + 2*J_ETC4_im_to_matrix - 1*J_HLEAK_im_to_matrix - 2*J_PIH_im_to_matrix + 1*J_KH_im_to_matrix + 1*J_HPERM_cytoplasm_to_im) / VWater_im;
J_M = (0 + 1*J_MPERM_cytoplasm_to_im) / VWater_im;
J_K = (0 - 1*J_KH_im_to_matrix + 1*J_KPERM_cytoplasm_to_im) / VWater_im;
Phi_H = J_H - sum( h_im*f(ii)'./(Kh(ii).*P(ii)) );
Phi_M = J_M - sum( m_im*f(ii)'./(Km(ii).*P(ii)) );
Phi_K = J_K - sum( k_im*f(ii)'./(Kk(ii).*P(ii)) );
% ALPHAs
% aH = 1 + pHBpH;
aM = 1 + pMBpM;f
aK = 1 + pKBpK;
% ADDITIONAL BUFFER for [H+]
aH = 1 + pHBpH + BX(3)/K_BX(3)/(1+h_im/K_BX(3))^2; % M
% DENOMINATOR
D = aH*pKBpM*pMBpK + aK*pHBpM*pMBpH + aM*pHBpK*pKBpH - ...
    aM*aK*aH - pHBpK*pKBpM*pMBpH - pHBpM*pMBpK*pKBpH;
% DERIVATIVES for H,Mg,K
f(26,:) =  ( (pKBpM*pMBpK - aM*aK)*Phi_H + ...
            (aK*pHBpM - pHBpK*pKBpM)*Phi_M + ...
            (aM*pHBpK - pHBpM*pMBpK)*Phi_K ) / D;
f(27,:) =  ( (aK*pMBpH - pKBpH*pMBpK)*Phi_H + ...
            (pKBpH*pHBpK - aH*aK)*Phi_M + ...
            (aH*pMBpK - pHBpK*pMBpH)*Phi_K ) / D;
f(28,:) =  ( (aM*pKBpH - pKBpM*pMBpH)*Phi_H + ...
            (aH*pKBpM - pKBpH*pHBpM)*Phi_M + ...
            (pMBpH*pHBpM - aH*aM)*Phi_K ) / D;

%% ELECTROPHYS EQUATIONS
C_im_to_matrix = par(20);
f(29) = ( 0 - 2.66667*J_F1F0ATPASE_im_to_matrix  + 4*J_ETC1_im_to_matrix  + 2*J_ETC3_im_to_matrix  + 4*J_ETC4_im_to_matrix  - 1*J_HLEAK_im_to_matrix  - 1*J_ANT_im_to_matrix ) / C_im_to_matrix;

%% FLUX VECTOR:
J = [ J_DH_matrix J_SDH_BBV_matrix J_ATPASE_cytoplasm J_CK_cytoplasm J_AK_cytoplasm J_F1F0ATPASE_im_to_matrix J_ETC1_im_to_matrix J_ETC3_im_to_matrix J_ETC4_im_to_matrix J_HLEAK_im_to_matrix J_PIH_im_to_matrix J_ANT_im_to_matrix J_KH_im_to_matrix J_AMPPERM_cytoplasm_to_im J_ADPPERM_cytoplasm_to_im J_ATPPERM_cytoplasm_to_im J_PIPERM_cytoplasm_to_im J_HPERM_cytoplasm_to_im J_KPERM_cytoplasm_to_im J_MPERM_cytoplasm_to_im];


%% Calculations to output [MgATP]c, [MgADP]c, and [Pi]c

% MgATP
unbound_ATP_c = ATP_c/P(8);
MgATP_c = unbound_ATP_c*m_c/Km(8);

% MgADP
unbound_ADP_c = ADP_c/P(9);
MgADP_c = unbound_ADP_c*m_c/Km(9);

%fPi_c
unbound_Pi_c = Pi_c/P(10);
% MgPi_c = unbound_Pi_c*m_c/Km(10);
fPi_c = unbound_Pi_c;% - MgPi_c;

if flag == 1
  f = f(:);
else
  f = [MgATP_c,MgADP_c,fPi_c,dGrATPase];
end
