function f = TrisegEquations(x,Vw_LV,Vw_SEP,Vw_RV,SL_LV,SL_SEP,SL_RV, V_LV, V_RV,Amref_LV,Amref_SEP,Amref_RV)

xm_LV  = x(1); % LV heart geometry variable, cm
xm_SEP = x(2); % septum heart geometry variable, cm
xm_RV  = x(3); % RV heart geometry variable, cm
ym     = x(4); % Heart geometry variable, cm

Lsref     = 1.9; % Resting SL, micron
Kse       = 50000; % series element elastance, mmHg/micron (Changed to match the value in Tewari's code) (9/5 BM)

%% ventricular mechanics
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
SLo_LV = Lsref*exp(epsf_LV); SLo_SEP = Lsref*exp(epsf_SEP); SLo_RV = Lsref*exp(epsf_RV);

% Total forces
sigmaf_LV = Kse*(SLo_LV - SL_LV);
sigmaf_SEP = Kse*(SLo_SEP - SL_SEP);
sigmaf_RV = Kse*(SLo_RV - SL_RV);
% equilibrium of forces at junction circle

Tm_LV = (Vw_LV*sigmaf_LV/(2*Am_LV))*(1 + (z_LV^2)/3 + (z_LV^4)/5);
Tm_SEP = (Vw_SEP*sigmaf_SEP/(2*Am_SEP))*(1 + (z_SEP^2)/3 + (z_SEP^4)/5);
Tm_RV = (Vw_RV*sigmaf_RV/(2*Am_RV))*(1 + (z_RV^2)/3 + (z_RV^4)/5);

Tx_LV = Tm_LV*2*xm_LV*ym/(xm_LV^2 + ym^2);
Tx_SEP = Tm_SEP*2*xm_SEP*ym/(xm_SEP^2 + ym^2);
Tx_RV = Tm_RV*2*xm_RV*ym/(xm_RV^2 + ym^2);

Ty_LV = Tm_LV*(-xm_LV^2 + ym^2)/(xm_LV^2 + ym^2);
Ty_SEP = Tm_SEP*(-xm_SEP^2 + ym^2)/(xm_SEP^2 + ym^2);
Ty_RV = Tm_RV*(-xm_RV^2 + ym^2)/(xm_RV^2 + ym^2);

f(1) = (-V_LV - 0.5*Vw_LV - 0.5*Vw_SEP + Vm_SEP - Vm_LV)/V_LV; %xm_LV
f(2) = (Tx_LV + Tx_SEP + Tx_RV); % xm_SEP
f(3) = (V_RV + 0.5*Vw_RV + 0.5*Vw_SEP + Vm_SEP - Vm_RV)/V_LV; %xm_RV
f(4) = (Ty_LV + Ty_SEP + Ty_RV);  %ym



