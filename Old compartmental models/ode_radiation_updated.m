%ode_radiation but now the effect only has an inward protective element,
%rather than an inward and outwards protective element. Now the model only
%means that a certain percentage of the population will be able to have the
%benefits of the radiation while they are inside rather than even after
%they leave. This is more accurate and a substantial improvement upon the
%previous model.

%note that the 



function dydt = ode_radiation(beta)

%gotta define the constants in here
%beta = 0.1986
N = 328.2e6;
eps = 38.4156/112.846;
pre_infec = 5.1;
sigma = 1/pre_infec;
eta = 0.5;
alpha = 0.5;
phi = 0.025;
gamma_A = 1/7;
gamma_I = 1/7;
gamma_H = 1/14;
delta = 0.015;


S = y(1);
E = y(2);
I = y(3);
A = y(4);
H = y(5);
R = y(6);
D = y(7);
S_R = y(8);
E_R = y(9);
I_R = y(10);
A_R = y(11);
H_R = y(12);
R_R = y(13);
D_R = y(14);

dS = -beta.*(I+eta.*A).*S./N - beta.*((1-0).*I_R+(1-0).*eta.*A_R).*S./N;
dE = beta.*(I+eta.*A).*S./N + beta.*((1-0)*I_R + (1-0)*eta*A_R).*S/N- sigma.*E;
dI = alpha.*sigma.*E-phi.*I-gamma_I.*I;
dA = (1-alpha).*sigma.*E-gamma_A.*A;
dH = phi.*I-delta.*H-gamma_H*H;
dR = gamma_I*I+gamma_A*A + gamma_H*H;
dD = delta*H;

dS_R = -beta.*(1-eps)*(I+eta*A).*S_R./N-beta.*(1-eps)*((1-0).*I_R + (1-0)*eta.*A_R)*S_R/N;
dE_R = beta.*(1-eps)*(I+eta*A).*S_R./N + beta.*(1-eps)*((1-0)*I_R+(1-0)*eta*A_R).*S_R./N - sigma.*E_R;
dI_R = alpha.*sigma.*E_R-phi.*I_R-gamma_I.*I_R;
dA_R = (1-alpha).*sigma.*E_R-gamma_A.*A_R;
dH_R = phi.*I_R-delta.*H_R-gamma_H*H_R;
dR_R = gamma_I*I_R+gamma_A*A_R + gamma_H*H_R;
dD_R = delta*H_R;

dydt = [dS;dE;dI;dA;dH;dR;dD;dS_R;dE_R;dI_R;dA_R;dH_R;dR_R;dD_R];
end