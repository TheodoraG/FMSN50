clear all;
close all;
clc;

gamma_1 = 1;
gamma_2 = 43/32;


%% 3
N = 1000; %sample size, iterations
no_steps = 10; %observations - number of steps
roads = cell(N,2,no_steps); %all possible roads (walks)
cn = zeros(1,no_steps);

%for different no of steps
for ii = 1:no_steps
    for jj = 1:N % main loop
        %a path is a sequence of no_steps points
        directions = ceil(4*rand(1,ii)); %the four possible directions for 2d
        
        roads(jj,1,ii) = {[0 0]}; %start in the center - 0,0
        for kk = 1:ii
            %compute neighbors
            switch directions(kk)
                case 1 %left neighbor
                    roads(jj,1,ii) = {[roads{jj,1,ii}; [(roads{jj,1,ii}(kk,1)-1) roads{jj,1,ii}(kk,2)]]};
                case 2 %up
                    roads(jj,1,ii) = {[roads{jj,1,ii}; [roads{jj,1,ii}(kk,1) (roads{jj,1,ii}(kk,2)+1)]]};
                case 3 %right
                    roads(jj,1,ii) = {[roads{jj,1,ii}; [(roads{jj,1,ii}(kk,1)+1) roads{jj,1,ii}(kk,2)]]};
                case 4 %down
                    roads(jj,1,ii) = {[roads{jj,1,ii}; [roads{jj,1,ii}(kk,1) (roads{jj,1,ii}(kk,2)-1)]]};
            end
        end
        %check if the road is self avoiding - it does not have repeated steps
        roads(jj,2,ii) = {length(roads{jj,1,ii}) == length(unique(roads{jj,1,ii},'rows'))};
    end
    selfAvoidingPaths = sum([roads{:,2,ii}] == 1);
    cn(ii) = round((selfAvoidingPaths/N) * 4^ii); %estimate
end

%% 4
N = 1000; %sample size, iterations
no_steps = 10; %observations - number of steps
roads = cell(N,3,no_steps); %all possible roads
cn_4 = zeros(1,no_steps);

%for different no steps
for ii = 1:no_steps
    for jj = 1:N % main loop
        lattice = zeros(2*ii+1, 2*ii+1);  %define the lattice\the field
        weight = 1; %initialize the weights
        roads(jj,1,ii) = {[0 0]}; %start in the center - 0,0
        
        %initial coordinates
        xx = ii + 1 + roads{jj,1,ii}(1,1);
        yy = ii + 1 + roads{jj,1,ii}(1,2);
        lattice(xx,yy) = 1;
        
        for kk = 1:ii
            %define the possible directions
            up    = lattice(xx,yy+1);
            down  = lattice(xx,yy-1);
            left  = lattice(xx-1,yy);
            right = lattice(xx+1,yy);
            
            %compute available neighbours
            neighbors = [1 1 1 1] - [left, up, right, down];
            
            %avoid self-loops - eliminate the path
            if (sum(neighbors)==0)
                roads(jj,2,ii)={0}; %mark as not being a self-avoiding walk
                roads(jj,3,ii) = {weight}; %store the weight. 
                break;
            end
            
            weight = weight * sum(neighbors); %compute weights
            
            %choose a direction for moving
            directions = find(rand < (cumsum(neighbors)/sum(neighbors)),1,'first'); %generate the whole path
            
            switch directions
                %compute neighbors
                case 1 %left neighbor
                    roads(jj,1,ii) = {[roads{jj,1,ii}; [(roads{jj,1,ii}(kk,1)-1) roads{jj,1,ii}(kk,2)]]};
                case 2 %up
                    roads(jj,1,ii) = {[roads{jj,1,ii}; [roads{jj,1,ii}(kk,1) (roads{jj,1,ii}(kk,2)+1)]]};
                case 3 %right
                    roads(jj,1,ii) = {[roads{jj,1,ii}; [(roads{jj,1,ii}(kk,1)+1) roads{jj,1,ii}(kk,2)]]};
                case 4 %down
                    roads(jj,1,ii) = {[roads{jj,1,ii}; [roads{jj,1,ii}(kk,1) (roads{jj,1,ii}(kk,2)-1)]]};
            end
            
            %update coordinates
            xx = ii + 1 + roads{jj,1,ii}(kk+1,1);
            yy = ii + 1 + roads{jj,1,ii}(kk+1,2);
            lattice(xx,yy) = 1; %update lattice
        end
        if length(roads{jj,1,ii})==(ii+1) %if the walk is self avoiding
            roads(jj,2,ii) = {1}; %mark as self avoiding
            roads(jj,3,ii) = {weight}; %store the weights
        end
    end
    weight = [roads{:,3,ii}];
    cn_4(ii) = round(mean(weight));
end


%% 5
N = 1000; %sample size, iterations
no_steps = 10; %observations - number of steps
roads = cell(N,3); %all possible roads\ walks
cn_5 = zeros(1,no_steps);
xx = zeros(1,N); %coordinates
yy = zeros(1,N);
weights = zeros(1,N); %weights

%initialization - simultate the first component
for i = 1:N
    roads(i,1) = {[0 0]}; %start in the center
    roads(i,2) = {[]}; %initial weights
    roads(i,3) = {zeros(2*no_steps+1, 2*no_steps+1)}; %store the lattice
    xx(i) = no_steps + 1 + roads{i,1}(1,1); %initial coordinates
    yy(i) = no_steps + 1 + roads{i,1}(1,2);
    roads{i,3}(xx(i),yy(i)) = 1; %initial weights
end
%for different no of steps
for ii = 1:no_steps
    for jj = 1:N % main loop
        %start in the center - 0,0
        xx(jj) = no_steps + 1 + roads{jj,1}(ii,1);
        yy(jj) = no_steps + 1 + roads{jj,1}(ii,2);
        
        %define possible moving directions
        up    = roads{jj,3}(xx(jj),yy(jj)+1);
        down  = roads{jj,3}(xx(jj),yy(jj)-1);
        left  = roads{jj,3}(xx(jj)-1,yy(jj));
        right = roads{jj,3}(xx(jj)+1,yy(jj));
        
        %compute available neighbours
        neighbors = [1 1 1 1] - [left, up, right, down];
        
        weights(jj) = sum(neighbors); %update the weights
        
        %avoid self-loops
        if (weights(jj)==0)
            roads(jj,2) = {[roads{jj,2}; weights(jj)]}; 
            continue; %continue with the next step
        end
        
        roads(jj,2) = {[roads{jj,2}; weights(jj)]}; %store weights
        
        %sample a move direction
        directions = find(rand < (cumsum(neighbors)/sum(neighbors)),1,'first'); %generate the whole path
        
        switch directions
            %compute neighbors
            case 1 %left
                roads(jj,1) = {[roads{jj,1}; [(roads{jj,1}(ii,1)-1) roads{jj,1}(ii,2)]]};
            case 2 %up
                roads(jj,1) = {[roads{jj,1}; [roads{jj,1}(ii,1) (roads{jj,1}(ii,2)+1)]]};
            case 3 %right
                roads(jj,1) = {[roads{jj,1}; [(roads{jj,1}(ii,1)+1) roads{jj,1}(ii,2)]]};
            case 4 %down
                roads(jj,1) = {[roads{jj,1}; [roads{jj,1}(ii,1) (roads{jj,1}(ii,2)-1)]]};
        end
        %update coordinates
        xx(jj) = no_steps + 1 + roads{jj,1}(ii+1,1);
        yy(jj) = no_steps + 1 + roads{jj,1}(ii+1,2);
        roads{jj,3}(xx(jj),yy(jj)) = 1; %update lattice
    end
    cn_5(ii) = round(prod((1/N).*sum([roads{:,2}],2)));
    %resample with replacement
    ind  = randsample(1:N,N,true, weights./sum(weights));
    roads(:,1)=roads(ind,1); %update the walk
    roads(:,3)=roads(ind,3); %update the lattice
    
end

%% 6
steps = 1:no_steps;
v = steps';
v(:,2) = log(steps'); 
v(:,3) = ones(length(steps),1);
slope = v\log(cn_5)';
mu_2 = exp(slope(1));
gamma_2 = slope(2)+1;
A_2 = exp(slope(3));


%% 9
d = 3; %dimension
N = 1000; %sample size, iterations
no_steps = 10; %observations - number of steps
moves = [eye(d); -eye(d)];
roads = cell(N,2); %all possible roads
cn_9 = zeros(1,no_steps);
weights = zeros(1,N);

%initialization
for i = 1:N
    roads(i,1) = {zeros(1,d)}; %start in the center
    roads(i,2) = {[]}; %initial weights
end
%for different no of steps
for ii = 1:no_steps
    for jj = 1:N % main loop
        %compute available neighbours
        neighbors = zeros(1,2*d);
        for kk = 1:2*d %compute the possible moving directions
            aux = moves(kk,:)+ roads{jj,1}(ii,:); %all the directions
            %free neighbors
            neighbors(kk) = sum(sum(aux==roads{jj,1},2)==d)~=1;
        end
        
        weights(jj) = sum(neighbors); %update weights
        
        %if there are no available neighbors
        if (weights(jj)==0)
            roads(jj,2) = {[roads{jj,2}; weights(jj)]};
            continue; %continue with the next step
        end
        
        %compute importance weights:
        roads(jj,2) = {[roads{jj,2}; weights(jj)]};
        
        %sample a move direction
        directions = find(rand < (cumsum(neighbors)/sum(neighbors)),1,'first'); %generate the whole path
        
        %update the movements / walk
        roads(jj,1) = {[roads{jj,1}; moves(directions,:)+roads{jj,1}(ii,:)]};
    end
    cn_9(ii) = round(prod((1/N).*sum([roads{:,2}],2)));
    %resample with replacement
    ind  = randsample(1:N,N,true, weights./sum(weights));
    roads(:,1)=roads(ind,1); %update the walk   
end

mu_calculation = 2*d - 1 - 1/(2*d) - 3/(2*d)^2 - 16/(2*d)^3 - 102/(2*d)^4;

if d~=4
    steps = 1:no_steps;
    v = steps';
    v(:,2) = log(steps');
    v(:,3) = ones(length(steps),1);
    slope_9 = v\log(cn_9)';
    mu_9 = exp(slope_9(1));
    gamma_9 = slope_9(2)+1;
    A_9 = exp(slope_9(3));
elseif d == 4
    steps = 1:no_steps;
    v = steps';
    %v(:,2) = log(steps');
    v(:,2) = ones(length(steps),1);
    slope_9 = v\log(cn_9)';
    mu_9 = exp(slope_9(1));
    %gamma_9 = slope_9(2)+1
    A_9 = exp(slope_9(2));
end


%% 10
load('population.mat')

N = 1000;
n = 50; %number of observation
tau = zeros(1,n+1); %filter expectation
w = zeros(N,1); %weights 
p = @(x,y) unifpdf(y, x/2, x); %observation density, for weights
part = unifrnd(1/5, 3/5, [N 1]); %initialization
w = p(part,Y(1)); %weighting
tau(1) = sum(part.*w)/sum(w); %estimation

%conf interval
[x_sorted,index]=sort(part); % sort data
cum_weight=cumsum(w(index))/sum(w); %cumulative normalized sum of the weights
% for sorted data
indexLower=find(cum_weight>=0.025,1); % index for lower 2.5% quantile
indexUpper=find(cum_weight>=0.975,1); % index upper 2.5% quantile
tauLower(1)=x_sorted(indexLower); % lower 2.5% quantile
tauUpper(1)=x_sorted(indexUpper); % upper 2.5% quantile

ind = randsample(N,N,true,w); %selection
part = part(ind);
B = unifrnd(1/2, 3, [N,n+1]);  %stochastic repoduction rate
for k = 1:n % main loop
    part = B(:,k+1).* part.*(1-part); %mutation
    w = p(part,Y(k + 1)); %weighting
    tau(k + 1) = sum(part.*w)/sum(w); %estimation
    ind = randsample(N,N,true,w); %selection
    part = part(ind);
    % conf interval
    [x_sorted,index]=sort(part); % sort data
    cum_weight=cumsum(w(index))/sum(w); %cumulative normalized sum of the weights
    % for sorted data
    Ilower=find(cum_weight>=0.025,1); % index for lower 2.5% quantile
    Iupper=find(cum_weight>=0.975,1); % index upper 2.5% quantile
    tauLower(k+1)=x_sorted(indexLower); % lower 2.5% quantile
    tauUpper(k+1)=x_sorted(indexUpper); % upper 2.5% quantile
end

figure()
hold on
plot(1:n+1, tau, 'LineWidth', 2.5);
plot(1:n+1, X, '--d','LineWidth', 2.5,'color','#77AC30');
hold off
xlabel('Simulations')
ylabel('Expectation')
legend('estimation','X')
title('Filter expectation estimate')

figure()
hold on
plot(1:n+1, tau, 'LineWidth', 2.5);
plot(1:n+1, tauLower, '--','LineWidth', 2.5,'color','#7E2F8E');
plot(1:n+1, tauUpper, '--','LineWidth', 2.5,'color','#EDB120');
plot(1:n+1, X, '--d','LineWidth', 2.5,'color','#77AC30');
hold off
xlabel('Simulations')
ylabel('Expectation')
legend('estimation','estimation - lower bound','estimation - upper bound','X')
title('Confidence interval')
