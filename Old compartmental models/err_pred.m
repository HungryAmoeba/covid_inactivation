function MSE = err_pred(betas)

M = csvread("US_covid_nums.csv");

US_confirmed = M(:,1);
US_deaths = M(:, 2);


[days_since_jan_21,~] = size(M);

tspan = 0:days_since_jan_21-1;
f = betas(tspan+1);

N = 328.2e6;
y0 = [N-10,9,1,0,0,0,0,0,0,0,0,0,0,0];
    
[t,y]= ode45(@(t,y)ode_radiation_time_dep(t,y,tspan,f),tspan,y0);
pred_cases = y(:,3) + y(:,5) + y(:,6) + y(:,7) + ...
            y(:,10) + y(:,12) + y(:,13) + y(:,14);

        MSE = sum(sum((pred_cases-US_confirmed).^2))/size(pred_cases,1);
end