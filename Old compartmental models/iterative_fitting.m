[days_since_jan_21,~] = size(M);

M = csvread("US_covid_nums.csv");

US_confirmed = M(:,1);

N = 328.2e6;
y0 = [N-10,9,1,0,0,0,0,0,0,0,0,0,0,0];

[days_since_jan_21,~] = size(M);
tspan = 0:days_since_jan_21-1;

for iteration = 1:100000
    f = full_betas(tspan+1);
    
    [t,y]= ode45(@(t,y)ode_radiation_time_dep(t,y,tspan,f),tspan,y0);
    pred_cases = y(:,3) + y(:,5) + y(:,6) + y(:,7) + ...
        y(:,10) + y(:,12) + y(:,13) + y(:,14);
    
    for iter = 1:days_since_jan_21 -1
        % if predicted outcome is less, raise beta
        if pred_cases(iter+1) < US_confirmed(iter+1)
            full_betas(iter) = full_betas(iter) - .00005;
        end
        
        % if pred outcome is greater, lower beta
        if pred_cases(iter+1) > US_confirmed(iter+1)
            full_betas(iter) = full_betas(iter) + .00005;
        end
        
    end
end