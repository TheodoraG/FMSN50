clear all;
close all;
clc;

load powercurve_D236.mat;
%P(6.5) % power for the wind speed 6.5 m/s

%% 2.a
wind_power = 0:0.1:35;
power_speed = P(wind_power);

figure;
plot(wind_power, power_speed,'LineWidth',2.5)
title("Power curve")
xlabel("Wind speed v [m/s]");
ylabel("Output power P(v) [W]");

figure;
pdf_wind = wblpdf(wind_power, 1.9, 11.7);
plot(pdf_wind,'LineWidth',2.5)
title("Probability density function")

figure;
cdf_wind = wblcdf(wind_power,1.9, 11.7);
plot(cdf_wind,'LineWidth',2.5)
title("Cumulative distribution function")

%basic Monte Carlo - implementation from the lecture
lambda = [11.7, 10.7, 10.1, 8.8, 8.6, 8.9, 8.6, 8.9, 10.0, 10.9, 11.7, 11.7];
k = [2.0, 2.0, 2.0, 1.9, 1.9, 1.9, 1.9, 1.9, 2.0, 1.9, 2.0, 2.0];

N = 10000;
step =100;
tau_N =zeros(100,12);
LB =zeros(100,12);
UB =zeros(100,12);
WB =zeros(100,12);
%for each month
for i = 1:12
    l= 1;
    %for different number of samples
    for j = 10:step:N
        %generate some random Weibull distrib numbers
        X = wblrnd(lambda(i),k(i), [1 j]);
        %compute the function phi as P
        %the expectation for each month
        tau_N(l,i) = mean(P(X));
        
        %compute confidence interval
        LB(l,i) = tau_N(l,i) - abs(norminv(0.975))*std(P(X))/(sqrt(j));
        UB(l,i) = tau_N(l,i) + abs(norminv(0.975))*std(P(X))/(sqrt(j));
        WB(l,i) = UB(l,i) - LB(l,i);
        l=l+1;
    end    
end

figure(4);
hold on
plot(10:step:N,tau_N(:,1),'LineWidth',2.5,'color','r')
plot(10:step:N,UB(:,1),'--','LineWidth',1.5,'color','g')
plot(10:step:N,LB(:,1),'--','LineWidth',1.5,'color','g')

%truncated version
a = 3;
b = 30;
R = rand(1,N);
tau_Ntr = zeros(100,12);
LB_tr = zeros(100,12);
UB_tr = zeros(100,12);
WB_tr = zeros(100,12);
%for each month
for ii = 1:12
    ll = 1;    
    Fa = wblcdf(a,lambda(ii),k(ii));
    Fb = wblcdf(b,lambda(ii),k(ii)); 
    %for different number of samples
    for jj = 10:step:N
        %compute the inverse according to prob 1b
        v = Fa + (Fb-Fa)*rand(1,jj);
        %compensation
        const  = Fb-Fa;
        est= wblinv(v,lambda(ii),k(ii))*const;
        tau_Ntr(ll,ii) = mean(P(est));
        
        %compute confidence interval
        LB_tr(ll,ii) = tau_Ntr(ll,ii) - abs(norminv(0.975))*std(P(est))/(sqrt(jj));
        UB_tr(ll,ii) = tau_Ntr(ll,ii) + abs(norminv(0.975))*std(P(est))/(sqrt(jj));
        WB_tr(ll,ii) = UB_tr(ll,ii) - LB_tr(ll,ii);
        ll=ll+1;
    end
end

figure(4);
hold on
plot(10:step:N,tau_Ntr(:,1),'LineWidth',2.5,'color','k')
plot(10:step:N,UB_tr(:,1),'--','LineWidth',1.5,'color','b')
plot(10:step:N,LB_tr(:,1),'--','LineWidth',1.5,'color','b')


figure(5);
hold on
plot(10:step:N,tau_N(:,1),'LineWidth',2.5,'color','r')
hold on
plot(10:step:N,UB(:,1),'--','LineWidth',1.5,'color','g')
hold on
plot(10:step:N,LB(:,1),'--','LineWidth',1.5,'color','g')
hold on
plot(10:step:N,tau_Ntr(:,1),'LineWidth',2.5,'color','k')
hold on
plot(10:step:N,UB_tr(:,1),'--','LineWidth',1.5,'color','b')
hold on
plot(10:step:N,LB_tr(:,1),'--','LineWidth',1.5,'color','b')
hold off
legend("Standard MC","UB - standard MC","LB - standard MC","Truncated MC",...
    "UB - truncated MC","LB - truncated MC")
xlabel("Number of samples");
ylabel("Estimate");


%% 2b
%control variate
tau_N_CV = zeros(100,12);
LB_CV = zeros(100,12);
UB_CV = zeros(100,12);
%for each month
for i = 1:12
    l= 1;
    %for different number of samples
    for j = 10:step:N
        %generate some random Weibull distrib numbers
        X = wblrnd(lambda(i),k(i), [1 j]);
        %the control variate is V
        Y = wblrnd(lambda(i),k(i), [1 j]);
        %take the expectation of the control variate
        Expect_Y = gamma(1+(1/k(i)))*lambda(i);
        c = cov(P(X), Y'); 
        %compute the optimal beta
        beta = -c(1,2)./var(Y);
        Z = P(X) + beta.*(Y' - Expect_Y);
        tau_N_CV(l,i) = mean(Z);
       
        %compute confidence interval
        LB_CV(l,i) = tau_N_CV(l,i) - abs(norminv(0.975))*std(Z)/(sqrt(j));
        UB_CV(l,i) = tau_N_CV(l,i) + abs(norminv(0.975))*std(Z)/(sqrt(j));
        l=l+1;
    end
end

figure(4);
hold on
plot(10:step:N,tau_N_CV(:,1),'LineWidth',2.5,'color','#7E2F8E')
plot(10:step:N,UB_CV(:,1),'--','LineWidth',1.5,'color','#EDB120')
plot(10:step:N,LB_CV(:,1),'--','LineWidth',1.5,'color','#EDB120')
%legend("Standard","Truncated","Control variate")

figure(6);
hold on
plot(10:step:N,tau_N(:,1),'LineWidth',2.5,'color','r')
hold on
plot(10:step:N,UB(:,1),'--','LineWidth',1.5,'color','g')
hold on
plot(10:step:N,LB(:,1),'--','LineWidth',1.5,'color','g')
hold on
plot(10:step:N,tau_Ntr(:,1),'LineWidth',2.5,'color','k')
hold on
plot(10:step:N,UB_tr(:,1),'--','LineWidth',1.5,'color','b')
hold on
plot(10:step:N,LB_tr(:,1),'--','LineWidth',1.5,'color','b')
hold on
plot(10:step:N,tau_N_CV(:,1),'LineWidth',2.5,'color','#7E2F8E')
hold on
plot(10:step:N,UB_CV(:,1),'--','LineWidth',1.5,'color','#EDB120')
hold on
plot(10:step:N,LB_CV(:,1),'--','LineWidth',1.5,'color','#EDB120')
hold off
legend("Standard MC","UB - standard MC","LB - standard MC","Truncated MC",...
    "UB - truncated MC","LB - truncated MC","Control variate",...
    "UB - control variate ","LB - control variate")
xlabel("Number of samples");
ylabel("Estimate");


%% 2.c 
%importance sampling
tau_IS = zeros(100,12);
LB_IS = zeros(100,12);
UB_IS = zeros(100,12);
WB_IS = zeros(100,12);

%find mean and sigma for g distribution
mu = 11;
sigma =7;

x_new  = 0:0.5:35;
g_0 = normpdf(x_new,mu,sigma);
f_0 = wblpdf(x_new, lambda(1), k(1));
phi_0 = 1e-8*P(x_new);
omega_0 = f_0./g_0;

figure()
hold on
plot(x_new,g_0,'LineWidth',2.5)
hold on
plot(x_new,f_0,'LineWidth',2.5)
hold on
plot(x_new,phi_0,'LineWidth',2.5)
hold on
plot(x_new,phi_0.*omega_0','LineWidth',2.5)
hold off
legend("g","f","phi","phi*omega")

%for each month
for i3 = 1:12   
    l =1;
    %for different number of samples
    for j3 = 10:step:N
        %generate some normal distribution numbers
        X_N = sigma*randn(1,N)+mu;
        %f is the weibull distributon
        f = wblpdf(X_N, lambda(i3), k(i3));
        %let g be the normal distribution
        g = normpdf(X_N,mu,sigma);
        %compute the weight function
        omega = f./g;
        
        tau_IS(l,i3) = mean(P(X_N).*omega');
        %compute confidence interval
        LB_IS(l,i3) = tau_IS(l,i3) - abs(norminv(0.975))*std(P(X_N).*omega')/(sqrt(j3));
        UB_IS(l,i3) = tau_IS(l,i3) + abs(norminv(0.975))*std(P(X_N).*omega')/(sqrt(j3));
        WB_IS(l,i3) = UB_IS(l,i3) - LB_IS(l,i3);
        l = l + 1;
    end
end

figure(4);
hold on
plot(10:step:N,tau_IS(:,1),'LineWidth',2.5,'color','m')
plot(10:step:N,UB_IS(:,1),'--','LineWidth',1.5,'color','#4DBEEE')
plot(10:step:N,LB_IS(:,1),'--','LineWidth',1.5,'color','#4DBEEE')
%legend("Standard","Truncated","Control variate","Importance sampling")

figure(8);
hold on
plot(10:step:N,tau_N(:,1),'LineWidth',2.5,'color','r')
hold on
plot(10:step:N,UB(:,1),'--','LineWidth',1.5,'color','g')
hold on
plot(10:step:N,LB(:,1),'--','LineWidth',1.5,'color','g')
hold on
plot(10:step:N,tau_Ntr(:,1),'LineWidth',2.5,'color','k')
hold on
plot(10:step:N,UB_tr(:,1),'--','LineWidth',1.5,'color','b')
hold on
plot(10:step:N,LB_tr(:,1),'--','LineWidth',1.5,'color','b')
hold on
plot(10:step:N,tau_N_CV(:,1),'LineWidth',2.5,'color','#7E2F8E')
hold on
plot(10:step:N,UB_CV(:,1),'--','LineWidth',1.5,'color','#EDB120')
hold on
plot(10:step:N,LB_CV(:,1),'--','LineWidth',1.5,'color','#EDB120')
hold on
plot(10:step:N,tau_IS(:,1),'LineWidth',2.5,'color','m')
hold on
plot(10:step:N,UB_IS(:,1),'--','LineWidth',1.5,'color','#4DBEEE')
hold on
plot(10:step:N,LB_IS(:,1),'--','LineWidth',1.5,'color','#4DBEEE')

hold off
legend("Standard MC","UB - standard MC","LB - standard MC","Truncated MC",...
    "UB - truncated MC","LB - truncated MC","Control variate",...
    "UB - control variate ","LB - control variate", "Importance sampling",...
    "UB - importance sampling ","LB - importance sampling")
xlabel("Number of samples");
ylabel("Estimate");

%% 2.d 
%antithetic sampling
tau_AS = zeros(100,12);
LB_AS = zeros(100,12);
UB_AS = zeros(100,12);
WB_AS = zeros(100,12);
%for each month
for ii = 1:12
    ll = 1;    
    Fa = wblcdf(3.0,lambda(ii),k(ii));
    Fb = wblcdf(30.0,lambda(ii),k(ii));  
    
    %for dfferent number of samples
    for jj = 10:step:N
        %compute the V variables for each interval
        R = rand(1,round(jj/2));
        V1 = Fa + (Fb - Fa)*R;
        V2 = Fa + (Fb - Fa)*(1-R);
        %compensation constant
        const  = Fb - Fa;
        V1_inv = wblinv(V1, lambda(ii),k(ii))*const;
        V2_inv = wblinv(V2, lambda(ii),k(ii))*const;
        %check if the covariance is negative (check for January for
        %example)
        if ll == 100 && ii ==1
            Cov1 = sign(cov(V1_inv,V2_inv));
        end
        W = (P(V1_inv)+P(V2_inv))/2;
        tau_AS(ll,ii) = mean(W);
        
        %compute confidence interval
        LB_AS(ll,ii) = tau_AS(ll,ii) - abs(norminv(0.975))*std(W)/(sqrt(jj));
        UB_AS(ll,ii) = tau_AS(ll,ii) + abs(norminv(0.975))*std(W)/(sqrt(jj));
        WB_AS(ll,ii) = UB_AS(ll,ii) - LB_AS(ll,ii);
        ll=ll+1;
    end
end

figure(4);
hold on
plot(10:step:N,tau_AS(:,1),'LineWidth',2.5,'color','c')
plot(10:step:N,UB_AS(:,1),'--','LineWidth',1.5,'color','b')
plot(10:step:N,LB_AS(:,1),'--','LineWidth',1.5,'color','b')
hold off
legend("Standard MC","UB - standard MC","LB - standard MC","Truncated MC",...
    "UB - truncated MC","LB - truncated MC",...
    "Control variate","UB - Control variate","LB - control variate",...
    "Importance sampling","UB - Importance sampling","LB - importance sampling",...
    "Antithetic sampling","UB - antithetic sampling","LB - antithetic sampling")
xlabel("Number of samples");
ylabel("Estimate");

%% 2e
size = 10:100:10000;
%for each month
for  i =1:12
    %for different number of samples
    for j = 1:length(size)
        %generate Weibull numbers
        wind_power = wblrnd(lambda(i),k(i),[1 size(j)]);
        %compute the power
        power_speed = P(wind_power);
        %find the number of postive powers
        positive_powers = length(find(power_speed~=0));
        %compute the probability
        probability(i,j) = positive_powers/size(j);
    end
end
figure()
plot(10:100:10000,probability,'LineWidth',2.5)
title("Estimation of P>0")
legend("Jan","Feb","March","April","May","June","July","Aug","Sep","Oct","Nov","Dec");
xlabel("Number of samples");
ylabel("Probability");

%% 2f
m = 3;
ro = 1.225;
d = 236;
const = 0.5*ro*pi*d^2/4;

%for each month
for i =1:12
	%compute P total
	P_tot(i) = const*gamma(1+m/k(i))*lambda(i)^m;
end

%for each month
for i =1:12
	%compute the ratio
	ratio(i) = tau_IS(100,i)/P_tot(i);
        
	%compute confidence interval
	LB_R(i) = LB_CV(100,i)/P_tot(i);
	UB_R(i) = UB_CV(100,i)/P_tot(i);
end

ratio
figure()
plot(1:1:12,ratio,'o','LineWidth',2.5)
title("Ratio by months")
xlabel("Months");
ylabel("Ratio");


%% 2g

max_power = 15*10^6;
for i =1:length(1:1:12)
    %compute the capacitity factor as the estimate
    %diivded by the max power
    capacity_factor(i) = tau_IS(100,i)/max_power;
    availability_factor(i) = probability(i,100);
end
capacity_factor_year = mean(capacity_factor)
availability_factor_year = mean(availability_factor)

%%

% Problem 3 - HA1
clear;clc;
load('powercurve_D236.mat')

maxN = 10000;
jump = 100;
alpha = 0.638;
p = 3;
q = 1.5;
k = 1.95;
lambda = 10.05;

fweib = @(x) wblpdf(x,lambda,k);
Fweib = @(x) wblcdf(x,lambda,k);

fcomb = @(X) fweib(X(:,1)).*fweib(X(:,2)).* ...
    (1 + alpha.*((1-Fweib(X(:,1)).^p).^(q-1)).*((1-Fweib(X(:,2)).^p).^(q-1)).* ...
    ((Fweib(X(:,1)).^p).*(1+p.*q)-1).*((Fweib(X(:,2)).^p).*(1+p.*q)-1))

%% Tuning the Weight Function -- Part A

X = 0:0.01:35;

mu = [1 1]*11.5;
sigma = eye(2)*21.5;

Pwr = P(X)*1e-9;
fcomb_test = fcomb([X' X']);
g = mvnpdf([X' X'],mu,sigma);
omega = fcomb_test./g;
Pwr_omega = Pwr.*omega;
Pwr_fcomb = Pwr.*fcomb([X',X'])*100;

figure(1)
hold on
plot(X,Pwr)
plot(X,fcomb_test)
plot(X,g)
% plot(X,omega)
% plot(X,Pwr_omega)
plot(X,Pwr_fcomb)
hold off
legend("Pwr","fcomb","g","Pwr*fcomb")

%% Monte Carlo - IS -- Part A
iters = jump:jump:maxN;
tau = zeros(1,length(iters));
for N = iters
    i = N/jump;
    
    X = mvnrnd(mu,sigma,N);
    g = mvnpdf(X,mu,sigma); 

    P1 = P(X(:,1));
    P2 = P(X(:,2));

    tau(i) = mean((P1+P2).*(fcomb(X)./g));
end

% figure(2)
% plot(iters,tau)
% title("Estimated Combined Power Generation vs. N Samples")
% xlabel("Number of Samples (N)")
% ylabel("Power Generated (1e7 W)")
% legend("IS Monte Carlo")

clc
fprintf('E(P(V1)+P(V2)): %.4f\n',tau(end)/10^6)

%% Covariance -- Part B
N=1000000;

X = mvnrnd(mu,sigma,N);
g = mvnpdf(X,mu,sigma); 

P1 = P(X(:,1));
P2 = P(X(:,2));

E_P1 = mean( P1.*(fcomb(X)./g) );
E_P2 = mean( P2.*(fcomb(X)./g) );
E_P1P2 = mean((P1.*P2).*(fcomb(X)./g));

COV = E_P1P2 - E_P1*E_P2;

% Variance & Std. Dev. -- Part C
P1weib = P1.*(fcomb(X)./g);
P2weib = P2.*(fcomb(X)./g);

VAR = sum(((P1weib+P2weib)-tau(end)).^2)/(N-1);
STD = sqrt(VAR);

clc
fprintf('Cov: %.2fe13\nVar: %.2fe13\nStd: %.2fe6\n',COV/10^13,VAR/10^13,STD/10^6)

%% Cutoff Energy Production -- Part D

% Finding new instrumental distribution to solve unbounded weights when X > 15MW.

X = 0:0.01:35;

mu = [1 1]*11.5;
sigma = eye(2)*21.5; 
sigma2 = eye(2)*50.5; % Increasing sigma seems to restrain the weights

Pwr = P(X)*1e-9;
fcomb_test = fcomb([X' X']);
g = mvnpdf([X' X'],mu,sigma2);
omega = fcomb_test./g;
Pwr_omega = Pwr.*omega;
Pwr_fcomb = Pwr.*fcomb([X',X'])*100;

figure(1)
hold on
plot(X,Pwr)
plot(X,fcomb_test)
plot(X,g)
% plot(X,omega)
plot(X,Pwr_omega)
plot(X,Pwr_fcomb)
hold off
legend("Pwr","fcomb","g","Pwr*omega","Pwr*fcomb")

%% Calculations -- Part D
N=1000000;

X = mvnrnd(mu,sigma,N);
g = mvnpdf(X,mu,sigma);
X2 = mvnrnd(mu,sigma2,N);
g2 = mvnpdf(X2,mu,sigma2);

P1 = P(X(:,1));
P2 = P(X(:,2));

P1_2 = P(X2(:,1));
P2_2 = P(X2(:,2));

tau_below = mean(((P1+P2)<15e6).*(fcomb(X)./g));
tau_above = mean(((P1_2+P2_2)>15e6).*(fcomb(X2)./g2));

clc
above_U = tau_above + abs(norminv(0.975)*std((P1_2+P2_2)>15e6)/sqrt(N))
above_L = tau_above - abs(norminv(0.975)*std((P1_2+P2_2)>15e6)/sqrt(N))
below_U = tau_below + abs(norminv(0.975)*std((P1+P2)<15e6)/sqrt(N))
below_L = tau_below - abs(norminv(0.975)*std((P1+P2)<15e6)/sqrt(N))
