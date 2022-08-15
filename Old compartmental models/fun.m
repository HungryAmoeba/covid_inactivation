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