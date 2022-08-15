function dydt = ode_fun_simple(t,y,beta)

Death = 0.034;
Pre_infec = 5.2;
f = 1/Pre_infec;
Duration = 14;
r = 1/Duration;
S = y(1);
E = y(2);
I = y(3);

dS = -beta*I.*S;
dE = beta*I.*S-f.*E;
dI = f*E - r*I - Death*I;
dR = r*I;
dD = Death*I;

dydt=[dS;dE;dI;dR;dD];

end