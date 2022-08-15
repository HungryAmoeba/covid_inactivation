
global beta gamma N % parameters 
global S0 I0 R0 % initial values (normalized S0=1)

global C
global init 
close all

% get data
[CC, date0] = getData();

fprintf('%12s %3s %10s %10s %10s %7s %7s %7s %7s\n', ...
   'Date', 'Day', 'N', 'Sinf', 'Rinf', 'beta', 'gamma', 'R0', 'R2' );
init = false;
for ii = 28:length(CC)
    
    % get data
    C = CC(1:ii);
    
    % estimate model parameters
    b = parest();
    
    % final number
    Rinf = Rmax();
    Sinf = S0*exp(-beta/gamma/N*(Rinf - R0));
    
    %calculate R2
    tspan = 0:length(C) - 1;   % final time
    ic = [S0 I0 R0]'; %initial conditions
    opts = [];   % no options set
    [~,z] = ode45(@SIR, tspan, ic, opts);
    
    % pretty sure this is calcuation of cumulative cases
    z = z(:,2) + z(:,3); 
    
    zbar = sum(C)/length(C);
    SStot = sum((C-zbar).^2);
    SSres = sum((C-z').^2);
    R2 = 1 - SSres/SStot;

    fprintf('%12s %3d %10d %10d %10d %7.3f %7.3f %7.3f %7.3f\n',...
        datestr(date0+ceil(length(C) - 1)), ceil(length(C)), round(N,0), ...
        round(Sinf,0),round(Rinf,0),beta,gamma,beta/gamma,R2)
end

%set parameters
tspan = 0:2*length(C);  %final time
ic = [S0 I0 R0]';   % initial conditions
opts = [];

% simulate
[t,z] = ode45(@SIR, tspan, ic, opts);

%plot results
figure
hold on
plot(t,z,'LineWidth', 2)
legend('Susceptible','Infected', 'Recovered', ...
    'Location',  'best', 'FontSize', 12)
xlabel('Day after jan 16 2020')
ylabel('Cases')
grid on
hold off
shg

% plot comparison
figure
hold on
plot(t, (z(:,2)+z(:,3))','k','LineWidth',2)
scatter(1:length(C), C, 50, 'filled')
legend('Predicted','Actua', ...
    'Location', 'best','FontSize',12)
xlabel('Day after jan 16 2020')
ylabel('Cases')
grid on
hold off
shg

% save to global
ta = t;

function dzdt = SIR(~,z)
%SIR model
    global beta gamma N
    S = z(1);
    I = z(2);
    R = z(3);
    N = S + I + R;
    dzdt = [-beta*I*S/N; beta*I*S/N - gamma*I; gamma*I];
end

function b = iniguess()
%iniguess obtains initial guess and defines beta, gamma, S0, I0 and R0
%values.
    global beta gamma
    global S0 I0 R0
    global C
    global init
    if ~init
        beta = 2.897;
        gamma = 2.689;
        S0 = 1e8;
        I0 = C(1);
        R0 = 0;
        init = true;
    end
    b(1) = beta;
    b(2) = gamma;
    b(3) = S0;
end

function b = parest
%PAREST is parameter estimation.
% fminsearch takes fun, original b0 estimation as inputs, and chooses b0
% values to minimize the function (norm between predicted and expected
% values)
    global beta gamma N
    global S0 I0 R0
    maxiter = 20000;
    maxfun = 20000;
    
    % this step just assigns base values to beta, gamma, S0, I0 and R0.
    b0 = iniguess();
    options = optimset('Display', 'off', 'MaxIter', maxiter, ...
        'MaxFunEvals', maxfun);
    [b, fmin, flag] = fminsearch(@fun, b0, options);
    fun(b); % obtain x ped
    beta = b(1);
    gamma = b(2);
    S0 = b(3);
    N = S0 + I0 + R0;
end

function f = fun(b)
%FUN Optimization function

%fun(b) runs on the given SIR equations and returns f, the difference
%between the calculated and actual number of cases as calculated by 
%the norm of the vector

    global beta gamma S0 I0 R0 
    global C Ca
    
    % set parameters
    beta = b(1);
    gamma = b(2);
    S0 = b(3);
    
    tend = length(C);
    
    tspan = 0:tend-1; %time interval
    ic = [S0 I0 R0]';
    
    %solve ODE
    try
            [tsol, zsol] = ode45(@SIR, tspan, ic);
    catch 
        f = NaN;
        return        
    end
    
    %check if calculation time equals sample time
    if length(tsol) ~= length(tspan)
        f = NaN;
        return
    end
    
    Ca = (zsol(:,2)+zsol(:,3))'; % calculated number of cases
    f = norm(C - Ca);
    
end


function r = Rmax()
%FSOLVE calculate number of recovered individuals after t=inf
    global N S0 beta gamma R0
    RN = beta/gamma;
    r = fzero(@f, [0,S0]);
    
    function z = f(x)
    z = x - (N-S0*exp(-RN*(x-R0)/N));
end
    

end

function [C, date0] = getData()
 date0=datenum('2020/01/16');
C = [ 45 62 121 198 291 440 580 845 1317 2015 2800 4581 6058 7813 9821 11948 14551 17389 20628 24553 28276 31439 34876 37552 40553 43099 44919 60326 64437 67100 69197 71329 73332 75198 75700 76676 77673 78651 79400 80088]'; 
 end