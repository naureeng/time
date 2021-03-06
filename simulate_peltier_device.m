function T = simulate_peltier_device
% T = simulate_peltier_device
%  Calculates steady-state temperatures generated by the Peltier device
%  Output (T) is the steady state temperature relative to baseline
%  Arrays of values can be assigned to any of the parameters
%
% by Dmitriy Aronov and Michale S. Fee
% "Analyzing the dynamics of brain circuits with temperature:
% design and implementation of a miniature thermoelectric device"
% Journal of Neuroscience Methods (2011)




% SEMICONDUCTOR PROPERTIES
% Seebeck coefficient (V/K)
alpha = 148e-6;
% Electrical resistivity (Ohm�m)
rho = 1.74e-5;
% Thermal conductivity (W/m�K)
kappa = 1.38;
% Geometric constant (m)
G = 1.7e-3;
% Number of semiconductor junctions
N = 1;

% HEAT SINK PROPERTIES
% Thermal conductance of the body-coupled heat sink (W/K)
K_B = 7e-3;
% Thermal conductance of the convective heat sink (W/K)
K_C = 5e-3;
% Temperature of the convective fluid relative to body (K)
deltaT_f = -18;

% PROBE PROPERTIES
% Diameter of the probe (m)
d = 250e-6;
% Length of the segment above the brain (m)
l0 = 1e-3;
% Length of the insulated segment (m)
l1 = 1.5e-3;
% Length of the un-insulated segment
l2 = 1e-3;
% Thermal conductivity of material (W/m�K)
kappa_P = 429;
% Number of probes
N_P = 2;

% INSULATION PROPERTIES
% Thermal conductivity of the cold-plate insulation
kappa_cold = 0.033;
% Area of the cold-plate insulation (m�)
A_cold = 5e-6;
% Thickness of the cold-plate insulation (m)
t_cold = 1e-3;
% Calculate thermal conductance of cold-plate insulation (W/K)
K_I = A_cold.*kappa_cold./t_cold;
% Outer diameter of the first insulation layer (m)
do = 370.84e-6;
% Thermal conductivity of the first insulation layer (W/m�K)
kappa1 = 0.025;
% Thickness of the second insulation layer (m)
t_poly = 25.4e-6;
% Thermal conductivity of the second insulation layer (W/m�K)
kappa2 = 0.12;

% BRAIN PROPERTIES
% Length constant (m)
lambda = 1.59e-3;
% Thermal conductivity (W/m�K)
kappa_B = 0.54;
% Baseline temperature (K)
T_B = 314;

% ELECTRIC CURRENT (A)
I = linspace(-10,1,1000);

% LOCATION
% Distance from probe center (m)
r = 0;
% Vectical position relative to brain surface (m)
y = -2e-3;


% DEVICE SIMULATION
% Calculate surface thermal resistances of the probe segments
rm2 = 1./(kappa_B./lambda.*besselk(1,(d./2)./lambda) ...
    ./besselk(0,(d./2)./lambda));
r1 = 1./(pi.*do.*kappa_B./lambda.*besselk(1,(do./2)./lambda) ...
    ./besselk(0,(do./2)./lambda));
r_air = log((do-2.*t_poly)./d)./(2.*pi.*kappa1);
r_poly = log(do./(do-2.*t_poly))./(2.*pi.*kappa2);
rm1 = (r_air + r_poly + r1).*(pi.*d);

% Calculate "electrotonic" lengths of the probe segments
lambda1 = sqrt(d.*rm1.*kappa_P./4);
lambda2 = sqrt(d.*rm2.*kappa_P./4);

% Calculate coefficients of the cable equation solution
A = 1./(cosh(l1./lambda1) + ...
    lambda1./lambda2.*tanh(l2./lambda2).*sinh(l1./lambda1));
B = 1./(lambda2./lambda1.*cosh(l1./lambda1)./tanh(l2./lambda2) + ...
    sinh(l1./lambda1));
C = 1./(cosh(l1./lambda1).*cosh(l2./lambda2) + ...
    lambda1./lambda2.*sinh(l1./lambda1).*sinh(l2./lambda2));

% Calculate the thermal resistance of the entire probe
R0 = 4.*l0./(pi.*d.^2.*kappa_P);
Rfull = (4.*lambda1./(pi.*d.^2.*kappa_P))./ ...
    (A.*sinh(l1./lambda1)+B.*cosh(l1./lambda1));
Rprobe = R0 + Rfull;

% Calculate the full thermal load resistance on the cold plate
Kl = N_P./Rprobe + K_I;

% Calculate steady state temperature of the cold plate
Tc = ((2.*N.*kappa.*G + K_C + K_B + 2.*N.*alpha.*I) .* ...
    (N.*rho.*I.^2./G + 2.*N.*alpha.*T_B.*I) + 2.*N.*kappa.*G .* ...
    (N.*rho.*I.^2./G - 2.*N.*alpha.*T_B.*I + deltaT_f.*K_C)) ./ ...
    ((2.*N.*kappa.*G + Kl - 2.*N.*alpha.*I) .* ...
    (2.*N.*kappa.*G + K_C + K_B + 2.*N.*alpha.*I) - (2.*N.*kappa.*G).^2);

% Calculate steady state temperatures at all probe segments
Tsurf = Tc.*Rfull./Rprobe;
T0 = Tsurf+(Tc-Tsurf).*y./t_cold;
T1 = Tsurf.*(A.*cosh((l1+y)./lambda1) + B.*sinh((l1+y)./lambda1));
T2 = Tsurf.*C.*cosh((l2+l1+y)./lambda2);

% Assign temperatures to corresponding vertical positions
T = nan(size(Tc));
if numel(y)==1
    y = y*ones(size(Tc));
end
f = find(y>=0 & y<=t_cold);
T(f) = T0(f);
f = find(y<0 & -y<=l1);
T(f) = T1(f);
f = find(-y>l1 & -y<=l1+l2);
T(f) = T2(f);

% Calculate temperature as a function of the radius
if numel(r)==1
    r = r*ones(size(Tc));
end
if numel(T)==1
    T = T*ones(size(r));
end
f = find(y<0 & -y<=l1 & r>do./2);
T(f) = T(f).*besselk(0,r(f)./lambda)./besselk(0,(do./2)./lambda).* ...
    (r1./(r_air+r_poly+r1));
f = find(-y>l1 & -y<=l1+l2 & r>d./2);
T(f) = T(f).*besselk(0,r(f)./lambda)./besselk(0,(d./2)./lambda);