clear;
close all;
clc;

load coal_mine_disasters.mat;

n1 = 751;

figure()
plot(T, 1:n1,'LineWidth',2.5)
title("Evolution of disasters over the years")
xlabel("Years");
ylabel("Disasters");

%% 1b and 1c
N = 10000; %number of samples

%initial and final brekpoints
t_start = 1658;
t_end = 1980;

d=6; %no of breakpoints

%rho = 0.03;
rho = 0.01*ones(d,1);
psi = 20;
tau = T;
burn_in_size = 1000;

for i = 2:d
    d = i;
    %initial breakpoints
    t_aux_final = 1980;
    t_aux_init = 1670;
    step = (t_aux_final-t_aux_init)/d;
    t_middle = t_aux_init:step:t_aux_final;
    t =[t_start t_middle(2:end-1) t_end]; %breakpoints vector
    
    updated_breakpoints = zeros(N, length(t));
    theta = gamrnd(2, 1/psi); %Gamma(2;psi)-hyperprior on theta
    lambda = gamrnd(2, 1/theta, 1, d); %Gamma(2,theta) prior on the intesities
    
    %do a burn-in first to get a stationary behavior
    for j = 1:burn_in_size
        %Gibbs sampling for theta
        d_new = length(lambda);
        %draw particle theta from the posterior f f(theta|lambda,t,tau)
        theta = gamrnd(2*d_new + 2, 1./(sum(lambda)+psi));
        %draw particle lambda from the posterior f f(lambda|tau,t,theta)
        lambda = lambda_posterior(theta, t, tau);
        %use Metropolis-Hastings algorithm to update the breakpoints
        [~, t] = metropolisHastings(lambda, t, tau, rho);
    end

    for j = 1:N
        %Gibbs sampling for theta
        d_new = length(lambda);
        %draw particle theta from the posterior f f(theta|lambda,t,tau)
        theta = gamrnd(2*d_new + 2, 1./(sum(lambda)+psi));
        %draw particle lambda from the posterior f f(lambda|tau,t,theta)
        lambda = lambda_posterior(theta, t, tau);
        %use Metropolis-Hastings algorithm to update the breakpoints
        [acc, t] = metropolisHastings(lambda, t, tau, rho);
        updated_breakpoints(j, :) = t;
    end
    
    
    figure
    plot(updated_breakpoints,'LineWidth',2.5)
    title(['Chain for d = ' num2str(d-1) ' breakpoints'])
    xlabel('Number of samples')
    ylabel('Years (coordinates) & Breakpoints (lines)')
    
    figure
    for j = 2:d
        histogram(updated_breakpoints(:,j), 20)
        hold on
    end
    title(['Histogram of the chain for d = ' num2str(d-1) ' breakpoints'])
    xlabel('Years (coordinates) & Breakpoints (histogram)')
    ylabel('Number of disasters')
    
    figure
    plot(T, 1:length(T),'LineWidth',2.5)
    hold on
    for j = 2:d
        line([mean(updated_breakpoints(:, j)) mean(updated_breakpoints(:, j))], [0 length(T)], 'Color', [rand rand rand],'LineWidth',2.5)
    end
    
    title('Number of disasters during 1658-1980')
    xlabel('Years (coordinates) & Breakpoints (lines)')
    ylabel('Number of disasters')
    axis([t_start T(end) 0 length(T)])
end

%% 1d - theta
clear
close all

load('coal_mine_disasters.mat')

N = 10000; %number of samples

%initial and final brekpoints
t_start = 1658;
t_end = 1980;
t_aux_final = 1980;
t_aux_init = 1670;

d=5; %no of breakpoints
tau = T;

step = (t_aux_final - t_aux_init)/d;
%t =t_start:step:t_end;
t_middle =t_aux_init:step:t_aux_final;
t =[t_start t_middle(2:end-1) t_end];

psi = 50;
rho = 0.01*ones(d,1);

mean_theta = zeros(psi, 1);
var_theta = zeros(psi, 1);
burn_in_size =1000;

for psi = 1:psi
    theta = gamrnd(2, 1/psi);
    lambda = gamrnd(2, 1/theta, 1, d);
    aux_theta = zeros(N, 1);
    
    %do a burn-in first to get a stationary behavior
    for j = 1:burn_in_size
        %Gibbs sampling for theta
        d_new = length(lambda);
        %draw particle theta from the posterior f f(theta|lambda,t,tau)
        theta = gamrnd(2*d_new + 2, 1./(sum(lambda)+psi));
        %draw particle lambda from the posterior f f(lambda|tau,t,theta)
        lambda = lambda_posterior(theta, t, tau);
        %use Metropolis-Hastings algorithm to update the breakpoints
        [~, t] = metropolisHastings(lambda, t, tau, rho);
    end
    
    for i = 1:N
        %Gibbs sampling for theta
        d_new = length(lambda);
        %draw particle theta from the posterior f f(theta|lambda,t,tau)
        theta = gamrnd(2*d_new + 2, 1./(sum(lambda)+psi));
        %draw particle lambda from the posterior f f(lambda|tau,t,theta)
        lambda = lambda_posterior(theta, t, tau);
        %use Metropolis-Hastings algorithm to update the breakpoints
        [~, t] = metropolisHastings(lambda, t, tau, rho);
        %updated_breakpoints(j, :) = t;
        aux_theta(i) = theta;
    end
    
    mean_theta(psi) = mean(aux_theta);
    var_theta(psi) = var(aux_theta);
end

figure
plot(mean_theta, '*')
title('Mean of \theta posterior when changing \psi')
xlabel('\psi')
ylabel('\theta')


figure
plot(var_theta, '*')
title('Variance of \theta posterior when changing \psi')
xlabel('\psi')
ylabel('\theta')


%% 1d - lambda
clear
close all

load('coal_mine_disasters.mat')
N = 10000; %number of samples

%initial and final brekpoints
t_start = 1658;
t_end = 1980;
t_aux_final = 1980;
t_aux_init = 1670;

d=5; %no of breakpoints
tau = T;

step = (t_aux_final - t_aux_init)/d;
%t =t_start:step:t_end;
t_middle =t_aux_init:step:t_aux_final;
t =[t_start t_middle(2:end-1) t_end];

psi = 50;
rho  = 0.01*ones(d,1);

mean_lambda = zeros(psi, d);
var_lambda = zeros(psi,d);
aux_lambda = zeros(N, d);

burn_in_size =1000;

for auxi = 1:psi
    theta = gamrnd(2, 1/psi);
    lambda = gamrnd(2, 1/theta, 1, d);
    aux_lambda = zeros(auxi, d);
    
    %do a burn-in first to get a stationary behavior
    for j = 1:burn_in_size
        %Gibbs sampling for theta
        d_new = length(lambda);
        %draw particle theta from the posterior f f(theta|lambda,t,tau)
        theta = gamrnd(2*d_new + 2, 1./(sum(lambda)+psi));
        %draw particle lambda from the posterior f f(lambda|tau,t,theta)
        lambda = lambda_posterior(theta, t, tau);
        %use Metropolis-Hastings algorithm to update the breakpoints
        [~, t] = metropolisHastings(lambda, t, tau, rho);
    end
    
    for i = 1:N
        %Gibbs sampling for theta
        d_new = length(lambda);
        %draw particle theta from the posterior f f(theta|lambda,t,tau)
        theta = gamrnd(2*d_new + 2, 1./(sum(lambda)+auxi));
        %draw particle lambda from the posterior f f(lambda|tau,t,theta)
        lambda = lambda_posterior(theta, t, tau);
        %use Metropolis-Hastings algorithm to update the breakpoints
        [~, t] = metropolisHastings(lambda, t, tau, rho);
        %updated_breakpoints(j, :) = t;
        
        aux_lambda(i, :) = lambda';
    end
    mean_lambda(auxi, :) = mean(aux_lambda);
    var_lambda(auxi,:) = var(aux_lambda);
end

psi1 = (1:psi)*1e-2;
psi_2 = zeros(d,psi);

for i = 1:d
   psi_2(i,:) = psi1; 
end
figure()
plot(psi_2(1,:),mean_lambda, '*')
title('Mean of \lambda posterior when changing \psi')
xlabel('\psi')
ylabel('\lambda')
legend('\lambda_1', '\lambda_2', '\lambda_3', '\lambda_4', '\lambda_5', 'Orientation', 'horizontal')

figure()
plot(psi_2(1,:),var_lambda, '*')
title('Variance of \lambda posterior when changing \psi')
xlabel('\psi')
ylabel('\lambda')
legend('\lambda_1', '\lambda_2', '\lambda_3', '\lambda_4', '\lambda_5', 'Orientation', 'horizontal')


%% 1d - t
clear
close all

load('coal_mine_disasters.mat')
N = 10000; %number of samples

%initial and final brekpoints
t_start = 1658;
t_end = 1980;
t_aux_final = 1980;
t_aux_init = 1670;

d=5; %no of breakpoints
tau = T;

step = (t_aux_final - t_aux_init)/d;
%t =t_start:step:t_end;
t_middle =t_aux_init:step:t_aux_final;
t =[t_start t_middle(2:end-1) t_end];

psi = 50;
%rho = 0.01;
rho = 0.01*ones(d,1);

t_aux = zeros(N, length(t));
mean_t = zeros(psi, length(t));
variance_t = zeros(psi, length(t));

burn_in_size = 1000;

for psi = 1:psi
    theta = gamrnd(2, 1/psi);
    lambda = gamrnd(2, 1/theta, 1, d);
    
    %do a burn-in first to get a stationary behavior
    for j = 1:burn_in_size
        %Gibbs sampling for theta
        d_new = length(lambda);
        %draw particle theta from the posterior f f(theta|lambda,t,tau)
        theta = gamrnd(2*d_new + 2, 1./(sum(lambda)+psi));
        %draw particle lambda from the posterior f f(lambda|tau,t,theta)
        lambda = lambda_posterior(theta, t, tau);
        %use Metropolis-Hastings algorithm to update the breakpoints
        [~, t] = metropolisHastings(lambda, t, tau, rho);
    end
    
    for i = 1:N
        %Gibbs sampling for theta
        d_new = length(lambda);
        %draw particle theta from the posterior f f(theta|lambda,t,tau)
        theta = gamrnd(2*d_new + 2, 1./(sum(lambda)+psi));
        %draw particle lambda from the posterior f f(lambda|tau,t,theta)
        lambda = lambda_posterior(theta, t, tau);
        %use Metropolis-Hastings algorithm to update the breakpoints
        [~, t] = metropolisHastings(lambda, t, tau, rho);
        %updated_breakpoints(j, :) = t;
        t_aux(i, :) = t;
    end
    mean_t(psi, :) = mean(t_aux);
    variance_t(psi, :) = var(t_aux);
end

figure
plot(mean_t(:, 2:end-1), '*')
title('Mean of the breakpoints t when changing \psi')
xlabel('\psi')
ylabel('Breakpoints')
legend('t_1', 't_2', 't_3', 't_4', 'Orientation', 'horizontal')

figure
plot(variance_t(:, 2:end-1), '*')
title('Variance of the breakpoints t when changing \psi')
xlabel('\psi')
ylabel('Variance of breakpoints')
legend('t_1', 't_2', 't_3', 't_4', 'Orientation', 'horizontal')

%% 1 e - lambda and theta
clear
close all

load('coal_mine_disasters.mat')
N = 10000; %number of samples

%initial and final brekpoints
t_start = 1658;
t_end = 1980;
t_aux_final = 1980;
t_aux_init = 1670;

d=5; %no of breakpoints
tau = T;

step = (t_aux_final - t_aux_init)/d;
%t =t_start:step:t_end;
t_middle =t_aux_init:step:t_aux_final;
t =[t_start t_middle(2:end-1) t_end];

auxi = 30;
psi = 20;
rho_2 = (1:auxi)*0.01;
rho = zeros(d,auxi);
for i = 1:d
   rho(i,:) = rho_2; 
end

mean_theta_2 = zeros(auxi, 1);
mean_lambda_2 = zeros(auxi, d);

aux_theta_2 = zeros(N, 1);
aux_lambda_2 = zeros(N, d);

burn_in_size = 1000;

for auxi = 1:auxi
    theta = gamrnd(2, 1/psi);
    lambda = gamrnd(2, 1/theta, 1, d);
    aux_theta_2 = zeros(auxi, 1);
    aux_lambda_2 = zeros(auxi, d);
    
    %do a burn-in first to get a stationary behavior
    for j = 1:burn_in_size
        %Gibbs sampling for theta
        d_new = length(lambda);
        %draw particle theta from the posterior f f(theta|lambda,t,tau)
        theta = gamrnd(2*d_new + 2, 1./(sum(lambda)+psi));
        %draw particle lambda from the posterior f f(lambda|tau,t,theta)
        lambda = lambda_posterior(theta, t, tau);
        %use Metropolis-Hastings algorithm to update the breakpoints
        [~, t] = metropolisHastings(lambda, t, tau, rho);
    end
    
    for i = 1:N
        %Gibbs sampling for theta
        d_new = length(lambda);
        %draw particle theta from the posterior f f(theta|lambda,t,tau)
        theta = gamrnd(2*d_new + 2, 1./(sum(lambda)+psi));
        %draw particle lambda from the posterior f f(lambda|tau,t,theta)
        lambda = lambda_posterior(theta, t, tau);
        %use Metropolis-Hastings algorithm to update the breakpoints
        [~, t] = metropolisHastings(lambda, t, tau, rho(:, auxi));
        
        aux_theta_2(i) = theta;
        aux_lambda_2(i, :) = lambda';
    end
    mean_theta_2(auxi) = mean(aux_theta_2);
    mean_lambda_2(auxi, :) = mean(aux_lambda_2);
end

figure()
plot(rho(1,:),mean_lambda_2, '*')
title('The mean of the posterior intensities \lambda when changing \rho')
xlabel('\rho')
ylabel('\lambda')

figure()
plot(rho(1,:),mean_theta_2, '*')
title('The mean of the posteriror \theta when changing \rho')
xlabel('\rho')
ylabel('\theta')

%% 1e - t
clear
close all

load('coal_mine_disasters.mat')
N = 10000; %number of samples

%initial and final brekpoints
t_start = 1658;
t_end = 1980;
t_aux_final = 1980;
t_aux_init = 1670;

d=5; %no of breakpoints
tau = T;

step = (t_aux_final - t_aux_init)/d;
%t =t_start:step:t_end;
t_middle =t_aux_init:step:t_aux_final;
t =[t_start t_middle(2:end-1) t_end];

psi = 20;
rho_new = [0.01 0.02 0.03 0.04];
auxi = length(rho_new);
rho = zeros(d,auxi);
for i = 1:d
   rho(i,:) = rho_new; 
end

aux_t = zeros(N, length(t));
burn_in_size = 1000;

for auxi = 1:auxi
    aux_t = zeros(N, length(t));
    theta = gamrnd(2, 1/psi);
    lambda = gamrnd(2, 1/theta, 1, d);
    
    %do a burn-in first to get a stationary behavior
    for j = 1:burn_in_size
        %Gibbs sampling for theta
        d_new = length(lambda);
        %draw particle theta from the posterior f f(theta|lambda,t,tau)
        theta = gamrnd(2*d_new + 2, 1./(sum(lambda)+psi));
        %draw particle lambda from the posterior f f(lambda|tau,t,theta)
        lambda = lambda_posterior(theta, t, tau);
        %use Metropolis-Hastings algorithm to update the breakpoints
        [~, t] = metropolisHastings(lambda, t, tau, rho);
    end
    
    for i = 1:N
        %Gibbs sampling for theta
        d_new = length(lambda);
        %draw particle theta from the posterior f f(theta|lambda,t,tau)
        theta = gamrnd(2*d_new + 2, 1./(sum(lambda)+psi));
        %draw particle lambda from the posterior f f(lambda|tau,t,theta)
        lambda = lambda_posterior(theta, t, tau);
        %use Metropolis-Hastings algorithm to update the breakpoints
        [~, t] = metropolisHastings(lambda, t, tau, rho(:, auxi));
        aux_t(i, :) = t;
    end
    
    figure
    subplot(2,2,1)
    autocorr(aux_t(:, 2), 500)
    aux = rho(1, auxi);
    title(['Correlation for t_1 with \rho = ' num2str(aux)])
    xlabel('Time lag')
    ylabel('Dependency')
    
    subplot(2,2,2)
    autocorr(aux_t(:, 3), 500)
    aux = rho(1, auxi);
    title(['Correlation for t_2 with \rho = ' num2str(aux)])
    xlabel('Time lag')
    ylabel('Dependency')
    
    subplot(2,2,3)
    autocorr(aux_t(:, 4), 500)
    aux = rho(1, auxi);
    title(['Correlation for t_3 with \rho = ' num2str(aux)])
    xlabel('Time lag')
    ylabel('Dependency')
    
    subplot(2,2,4)
    autocorr(aux_t(:, 4), 500)
    aux = rho(1, auxi);
    title(['Correlation for t_4 with \rho = ' num2str(aux)])
    xlabel('Time lag')
    ylabel('Dependency')
    
    figure
    plot(aux_t,'LineWidth',2.5)
    aux = rho(1, auxi);
    title(['Chain for d = ' num2str(d-1) ' breakpoints and \rho =' num2str(aux)])
    xlabel('Number of samples')
    ylabel('Years (coordinates) & Breakpoints (lines)')
end

%% 1e - acceptance rate
close all

rho = linspace(0.001,0.1);
accepted = zeros(length(rho),d-1);

N2 = 100;
for j = 1:length(rho)
    
   for i = 1:N2
       %Gibbs sampling for theta
       d_new = length(lambda);
       %draw particle theta from the posterior f f(theta|lambda,t,tau)
       theta = gamrnd(2*d_new + 2, 1./(sum(lambda)+psi));
       %draw particle lambda from the posterior f f(lambda|tau,t,theta)
       lambda = lambda_posterior(theta, t, tau);
       %use Metropolis-Hastings algorithm to update the breakpoints
       [arr_acc, t] = metropolisHastings(lambda, t, tau, rho(j));
       accepted(j,:) = accepted(j,:) + arr_acc;
   end
   
end

ratio = sum(accepted,2)/(N2*(d-1));
plot(rho,ratio,'*');
hold on
plot([0,0.1],[0.3,0.3],'LineWidth',2.5);
title('Acceptance rate with respect to \rho');
xlabel('\rho');
ylabel('Acceptance percentage');

%% 2 b
clear
close all
load atlantic.txt

rng(10)
figure()
plot(1:length(atlantic), atlantic,'LineWidth',2.5)
title("Gumbel distribution of data")

%Parameter estimates for Gumbel data
[beta,mu] = est_gumbel(atlantic);
%expression for F inv
Finv = @(u, mu, beta) mu - beta.*log(-log(u));
N2 = 1000; %number of samples

boot_beta = zeros(N2,1);
boot_mu = zeros(N2,1);

%lecture 13, slide 33
for i = 1:N2 % bootstrap
    U = rand(length(atlantic),1);
    %compute F inv
    f_inv = Finv(U, mu, beta);
    %estimate parameters for F inv
    [beta_new, mu_new]  = est_gumbel(f_inv);
    %save the new parameters
    boot_beta(i) = beta_new;
    boot_mu(i) = mu_new;
end

delta_beta = sort(beta - boot_beta); % sorting to obtain quantiles
delta_mu = sort(mu - boot_mu); % sorting to obtain quantiles
%compute 95% confidence interval
%beta
alpha_2 = 0.05; % CB level
L_beta = beta - delta_beta(ceil((1-alpha_2/2)*N2)); 
U_beta = beta - delta_beta(ceil(alpha_2*N2/2));
%mu
L_mu = mu - delta_mu(ceil((1-alpha_2/2)*N2)); 
U_mu = mu - delta_mu(ceil(alpha_2*N2/2));

%% 2c
T1 = 3*14*100;

wave = Finv(1-1/T1, mu, beta);
wave_boot = zeros(N2,1);

for i = 1:N2 % bootstrap
    wave_boot(i) = Finv(1-1/T1, boot_mu(i), boot_beta(i));
end

deltaWave = sort(wave-wave_boot); % sorting to obtain quantiles

alpha5 = 0.95*N2;
%compute 95% one-sided confidence interval
CI_wave =[0, wave + deltaWave(ceil(alpha5))]; 

x = linspace(0,1,length(atlantic));
figure()
plot(x,f_inv,'LineWidth',2.5);
title("Inverse of the Gumbel distribution")