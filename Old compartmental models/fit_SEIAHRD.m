global beta eta sigma alpha phi delta gamma_A gamma_I gamma_H epsilon N
global S0 E0 I0 A0 H0 R0 D0 SI0 EI0 II0 AI0 HI0 RI0 DI0

global C
global init

close all
clearvars; clc;

% define some static parameters

eta = 0.5;
sigma = 1/5.1;
alpha = 0.5;
delta = .015;
epsilon = 38.4156/112.846;

% imports data from downloaded Johns Hopkins
M = csvread("US_covid_nums.csv");

CC = flip(M(1:173,1));
US_deaths = flip(M(1:173, 2));
US_recovered = flip(M(1:173, 3));


init = false;
for ii = 1:length(CC)
    C = CC(1:ii);
    
    b = parest();
    fprintf('%7.3f\n',beta)
%     STORE VALUES SOMEHOW
end

tspan = 0:2*length(C);
ic = [S0 E0 I0 A0 H0 R0 D0 SI0 EI0 II0 AI0 HI0 RI0 DI0]';
opts = [];

%simulate
[t,y] = ode45(@(t,z)ode_radiation_updated(), tspan, ic, opts);


function b = iniguess()
    global beta phi gamma_A gamma_I gamma_H 
    global S0 E0 I0 A0 H0 R0 D0 SI0 EI0 II0 AI0 HI0 RI0 DI0
    global C
    global init
    
    b = zeros(7,1);
    if ~init
        beta = .1986;
        phi = .025;
        gamma_A = 1/7;
        gamma_I = 1/7;
        gamma_H = 1/14;
        
        S0 = 328.2e6;
        E0 = 10;
        I0 = C(1);
        A0 = 0;
        H0 = 0;
        R0 = 0;
        D0 = 0;
        SI0 = 0;
        EI0 = 0;
        II0 = 0;
        AI0 = 0;
        HI0 = 0;
        RI0 = 0;
        DI0 = 0;
        init = true;     
    end
    b(1) = beta;
    b(2) = phi;
    b(3) = gamma_A;
    b(4) = gamma_I;
    b(5) = gamma_H;
    b(6) = S0;
    b(7) = E0;
    
end

function dydt = ode_radiation(~,y)

%gotta define the constants in here
%beta = 0.1986
global beta eta sigma alpha phi delta gamma_A gamma_I gamma_H epsilon N

N = S + E + I + A + H + R + D + S_R + E_R + I_R + A_R + H_R + R_R + D_R;

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

function b = parest
    global beta phi gamma_A gamma_I gamma_H N
    global S0 E0 I0 A0 H0 R0 D0 SI0 EI0 II0 AI0 HI0 RI0 DI0
    maxiter = 20000;
    maxfun =  20000;
    b0 = iniguess();
    
    options = optimset('Display', 'off', 'MaxIter', maxiter, ...
        'MaxFunEvals', maxfun);
    [b, fmin, flag] = fminsearch(@fun, b0, options);
    fun(b);
    beta = b(1);
    phi = b(2);
    gamma_A = b(3);
    gamma_I = b(4);
    gamma_H = b(5);
    S0 = b(6);
    E0 = b(7);
    N = S0 + E0 + I0 + A0 + H0 + R0 + D0 + SI0 + EI0 + II0 + AI0 + HI0 + RI0 + DI0;
end

function f = fun(b)
    global beta eta sigma alpha phi delta gamma_A gamma_I gamma_H epsilon N
    global S0 E0 I0 A0 H0 R0 D0 SI0 EI0 II0 AI0 HI0 RI0 DI0
    global C Ca
    
     beta = b(1);
    phi = b(2);
    gamma_A = b(3);
    gamma_I = b(4);
    gamma_H = b(5);
    S0 = b(6);
    E0 = b(7);
    
    tend = length(C);
    
    tspan = 0:tend-1;
    ic = [S0 E0 I0 A0 H0 R0 D0 SI0 EI0 II0 AI0 HI0 RI0 DI0]';
    
    %solve ODE
    try
        [tsol, zsol] = ode45(@SIR, tspan, ic);
    catch 
        f = NaN;
        return
        
    end
    
    if length(tsol) ~= length(tspan)
        f = NaN;return
    end
    
    Ca = zsol(:,3) + zsol(:,5) + zsol(:,6) + zsol(:,7) + zsol(:,10) + zsol(:,12) + zsol(:,13) + zsol(:,14);
    f = norm(C-Ca);
end
