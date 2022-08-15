% This program is going to try and calculate the best beta up until the end
% of each time series. Beginning values (t<20 will be discarded).

clearvars; clc;

% imports data from downloaded NYT
% https://github.com/nytimes/covid-19-data/blob/master/us.csv

M = csvread("US_covid_nums.csv");

US_confirmed = M(:,1);
US_deaths = M(:, 2);


[days_since_jan_21,~] = size(M);

beta_guess = .1986;
N = 328.2e6;
y0 = [N-10,9,1,0,0,0,0,0,0,0,0,0,0,0];

fitted_betas = zeros(days_since_jan_21, 1);

lower_b = .001;
upper_b = .75;
b_step = .001;

 for time = 1:days_since_jan_21
        tspan = 0:time;
        best_sq_error = 10000000000000;
        
        for test_beta = lower_b:b_step:upper_b
            [t,y]= ode45(@(t,y)ode_radiation(t,y,eps,test_beta),tspan,y0);
            
            %sums together all I, H, R and D
            pred_cases = y(:,3) + y(:,5) + y(:,6) + y(:,7) + ...
            y(:,10) + y(:,12) + y(:,13) + y(:,14);
            pred_cases = pred_cases(2:end, 1);
            % calculates sq error
            Sq_error = sum((US_confirmed(1:time,1)- ...
            pred_cases).^2)/10000; % divide by 100 to avoid big num
            
            % replaces best Sq error
            if (Sq_error < best_sq_error)
            fitted_betas(time) = test_beta;
            best_sq_error = Sq_error;
            end
        end
 end