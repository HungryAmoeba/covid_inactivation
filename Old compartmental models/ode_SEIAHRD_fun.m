function dydt = ode_SEIAHRD_fun(t,y,beta)

%gotta define the constants in here
beta = 0.5/50;
pre_infec = 5.1;
sigma = 1/pre_infec;
eta = 0.5;
alpha = 0.5;
phi = 0.025;
gamma_A = 1/7;
gamma_I = 1/7;
gamma_H = 1/14;
delta = 0.015;
N = 145.93e6;

S = y(1);
E = y(2);
I = y(3);
A = y(4);
H = y(5);
R = y(6);
D = y(7);

dS = -beta.*(I+eta.*A).*S./N;
dE = beta.*(I+eta.*A).*S./N-sigma.*E;
dI = alpha.*sigma.*E-phi.*I-gamma_I.*I;
dA = (1-alpha).*sigma.*E-gamma_A.*A;
dH = phi.*I-delta.*H-gamma_H*H;
dR = gamma_I*I+gamma_A*A + gamma_H*H;
dD = delta*H;

dydt = [dS;dE;dI;dA;dH;dR;dD];
end