

M = csvread("US_covid_nums.csv");
[days_since_jan_21,~] = size(M);
US_confirmed = M(:,1);

% fits for beta gamma and delta constants
b_min = -.1225;
b_0 = .6810;
r = .0122;
gamma_0 = .0149;
gamma_1 = 1.8971;
T_g = 26.4170;
delta_0 = -.2912;
delta_1 = .2001;
T_d = -2.5443;
alpha = .5027;
phi = .025;

func_inputs(1) = b_min;
func_inputs(2) = b_0;
func_inputs(3) = r;
func_inputs(4) = gamma_0;
func_inputs(5) = gamma_1;
func_inputs(6) = T_g;
func_inputs(7) = delta_0;
func_inputs(8) = delta_1;
func_inputs(9) = T_d;
func_inputs(10) = alpha;
func_inputs(11) = phi;


tspan = 0:days_since_jan_21-1;
% define beta, gamma and delta
fbeta = b_min + (b_0 - b_min)*exp(-r*(tspan+1));

% for gamma function, optimize gamma_0, gamma_1 and T_g
fgamma = gamma_0 + exp(-gamma_1*(tspan+1 +T_g));

% for delta function, optimize delta_0, delta_1 and T_d
fdelta = delta_0 + exp(-delta_1*(tspan+T_d));

% Half intervention



% half intervention
N = 328.2e6;
eps = 38.4156/112.846;
y0 = [N/2-10,9,1,0,0,0,0,N/2,0,0,0,0,0,0];

[t,y]= ode45(@(t,y)ode_tspan(t,y,tspan, fbeta,fgamma, ...
    fdelta,alpha,phi,eps),tspan,y0);

pred_cases = y(:,3) + y(:,5) + y(:,6) + y(:,7) + ...
            y(:,10) + y(:,12) + y(:,13) + y(:,14);

% no intervention
y0 = [N-10,9,1,0,0,0,0,0,0,0,0,0,0,0];
[t,y]= ode45(@(t,y)ode_tspan(t,y,tspan, fbeta,fgamma, ...
    fdelta,alpha,phi,eps),tspan,y0);

pred_no_int = y(:,3) + y(:,5) + y(:,6) + y(:,7) + ...
            y(:,10) + y(:,12) + y(:,13) + y(:,14);



% Plot results
plot(t,pred_cases,'LineWidth',1.5,'MarkerSize',18);
hold on;
plot(t,pred_no_int,'LineWidth',1.5,'MarkerSize',18);
hold on;
plot(t, US_confirmed, 'r*');
legend('Half exposure','No exposure','Confirmed','Location','best')
xlabel('Days after Jan 22, 2020')
ylabel('Confirmed cases')
title('Predicted spread half intervention')
grid on;
grid minor;


% By percent intervention
peakhospitalnum = zeros(101,1);
peakhospitalday = zeros(101,1);

totaldeaths = zeros(101,1);

for percent_rad = 0:100
    
    y0 = [N*(100-percent_rad)/100-10,9,1,0,0,0,0,N*percent_rad/100-2,2,0,0,0,0,0];

    [t,y]= ode45(@(t,y)ode_tspan(t,y,tspan, fbeta,fgamma, ...
    fdelta,alpha,phi,eps),tspan,y0);
    
    hosp = y(:,5) + y(:,12);
    [peakhospitalnum(percent_rad+1),peakhospitalday(percent_rad+1)]=max(hosp);
    
    deaths(percent_rad + 1) = y(end,7) + y(end,14);
end

%Plot results
figure
plot(0:100, peakhospitalnum, 'LineWidth',1.5,'MarkerSize',18);
xlabel('Percent individuals exposed to RF')
ylabel('Population')
legend('Peak hospitalization')
hold off;

figure
plot(0:100, peakhospitalday, 'LineWidth',1.5,'MarkerSize',18);
xlabel('Percent individuals exposed to RF')
ylabel('Days after march 12')
legend('Peak hospitalization day')
hold off;

figure
plot(0:100, deaths, 'LineWidth',1.5,'MarkerSize',18);
xlabel('Percent individuals exposed to RF')
ylabel('Deaths')
legend('Total deaths')
hold off;

% Plot peak hospitalizations
hosp = zeros(days_since_jan_21,6);

for iter = 1:6
    percent_rad = (iter-1)*20;
    
    y0 = [N*(100-percent_rad)/100-10,9,1,0,0,0,0,N*percent_rad/100-2,2,0,0,0,0,0];
    [t,y]= ode45(@(t,y)ode_tspan(t,y,tspan, fbeta,fgamma, ...
    fdelta,alpha,phi,eps),tspan,y0);
    
    hosp(:,iter) = y(:,5) + y(:,12);
end

figure
plot(tspan, hosp, 'LineWidth', 1.5,  'MarkerSize', 18);
legend('0','20', '40','60','80','100')
xlabel('Days after Jan 22')
ylabel('Population')
title('Hospitalizations')

% Impact of varying time of inactivation
% The base area under the curve is 112.846

deaths_matrix = zeros(101,25);
peakhospitalnums_matrix = zeros(101,25);
peakhospitaldays_matrix = zeros(101,25);

% for inactivation_time = 5:25
%     root = 1000/61*lambertw(61*inactivation_time/1000);
%     fun = @(x) 6.8836*exp(-.061*x)-6.8836*x/inactivation_time;
%     area = integral(fun,0,root);
%     eps = area/112.846
% end

for percent_rad = 0:100
    
    
    for inactivation_time = 1:25
        
        y0 = [N*(100-percent_rad)/100-10,9,1,0,0,0,0,N*percent_rad/100-2,2,0,0,0,0,0];
        
        root = 1000/61*lambertw(61*inactivation_time/1000);
        fun = @(x) 6.8836*exp(-.061*x)-6.8836*x/inactivation_time;
        area = integral(fun,0,root);
        eps = area/112.846;
        
        [t,y]= ode45(@(t,y)ode_tspan(t,y,tspan, fbeta,fgamma, ...
    fdelta,alpha,phi,eps),tspan,y0);
    
    
    hosp = y(:,5) + y(:,12);
    [peakhospitalnums,peakhospitaldays]=max(hosp);
    peakhospitalnums_matrix(percent_rad+1, inactivation_time) = peakhospitalnums;
    peakhospitaldays_matrix(percent_rad+1, inactivation_time) = peakhospitaldays;
    
    deaths = y(end,7) + y(end,14);
    deaths_matrix(percent_rad+1, inactivation_time) = deaths;
    
    end
    

end

%Plot data

grid off; 
colormap jet;

figure
heatmap(peakhospitaldays_matrix)
title( ...
    'Day of peak hospitalization')
xlabel("Time for inactivation")
ylabel("Percent of population exposed")

figure
grid off; 
colormap jet;
heatmap(deaths_matrix)
title('Deaths')
xlabel("Time for inactivation")
ylabel("Percent of population exposed")

figure
grid off; 
colormap jet;
heatmap(peakhospitalnums_matrix)
title('Peak hospitalization')
xlabel("Time for inactivation")
ylabel("Percent of population exposed")


% differential equations for model
function dydt = ode_tspan(t,y,tspan,fbeta,fgamma,fdelta,alpha,phi,eps)

%gotta define the constants in here

N = 328.2e6;

pre_infec = 5.1;
sigma = 1/pre_infec;
eta = 0.5;

gamma = interp1(tspan,fgamma,t);
delta = interp1(tspan, fdelta, t);
beta = interp1(tspan,fbeta,t);

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
dI = alpha.*sigma.*E-phi.*I-gamma.*I;
dA = (1-alpha).*sigma.*E-gamma.*A;
dH = phi.*I-delta.*H-gamma*H;
dR = gamma*I+gamma*A + gamma*H;
dD = delta*H;

dS_R = -beta.*(1-eps)*(I+eta*A).*S_R./N-beta.*(1-eps)*((1-0).*I_R + (1-0)*eta.*A_R)*S_R/N;
dE_R = beta.*(1-eps)*(I+eta*A).*S_R./N + beta.*(1-eps)*((1-0)*I_R+(1-0)*eta*A_R).*S_R./N - sigma.*E_R;
dI_R = alpha.*sigma.*E_R-phi.*I_R-gamma.*I_R;
dA_R = (1-alpha).*sigma.*E_R-gamma.*A_R;
dH_R = phi.*I_R-delta.*H_R-gamma*H_R;
dR_R = gamma*I_R+gamma*A_R + gamma*H_R;
dD_R = delta*H_R;

dydt = [dS;dE;dI;dA;dH;dR;dD;dS_R;dE_R;dI_R;dA_R;dH_R;dR_R;dD_R];
end