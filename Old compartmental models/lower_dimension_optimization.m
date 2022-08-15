
global US_confirmed;
global days_since_jan_21;

M = csvread("US_covid_nums.csv");

[days_since_jan_21,~] = size(M);



US_confirmed = M(:,1);
maxiter = 20000;
maxfun =  80000;
full_betas = ones(days_since_jan_21,1)*.1988;

for iter = 1:100
options = optimset('Display', 'off', 'MaxIter', maxiter, ...
         'MaxFunEvals', maxfun,'PlotFcns',@optimplotfval);

     %random initialization. could speed up convergance. 
optimum_betas = fminsearch(@err_pred, full_betas, options);
full_betas = optimum_betas;
end

tspan = 0:days_since_jan_21-1;
f = full_betas(tspan+1);

N = 328.2e6;
y0 = [N-10,9,1,0,0,0,0,0,0,0,0,0,0,0];
    
[t,y]= ode45(@(t,y)ode_radiation_time_dep(t,y,tspan,f),tspan,y0);
pred_cases = y(:,3) + y(:,5) + y(:,6) + y(:,7) + ...
            y(:,10) + y(:,12) + y(:,13) + y(:,14);
plot(pred_cases);
hold on;
plot(US_confirmed);

function MSE = err_pred(betas)

global US_confirmed days_since_jan_21;




tspan = 0:days_since_jan_21-1;
f = betas(round((tspan+1)/20)+1);

N = 328.2e6;
y0 = [N-10,9,1,0,0,0,0,0,0,0,0,0,0,0];
    
[t,y]= ode45(@(t,y)ode_radiation_time_dep(t,y,tspan,f),tspan,y0);
pred_cases = y(:,3) + y(:,5) + y(:,6) + y(:,7) + ...
            y(:,10) + y(:,12) + y(:,13) + y(:,14);

        MSE = sum(sum((pred_cases-US_confirmed).^2))/size(pred_cases,1);
end

