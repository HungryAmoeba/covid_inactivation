
% imports data from downloaded NYT
% https://github.com/nytimes/covid-19-data/blob/master/us.csv

M = csvread("US_covid_nums.csv");

US_confirmed = M(:,1);
US_deaths = M(:, 2);


[days_since_jan_21,~] = size(M);




maxiter = 100000;
maxfun =  100000;

    
options = optimset('Display', 'off', 'MaxIter', maxiter, ...
         'MaxFunEvals', maxfun);

[optimum_betas, low_MSE] = fminsearch(@err_pred, rand(182,1)/3, options);

