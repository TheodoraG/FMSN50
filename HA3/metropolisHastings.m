function [accepted_no, t] = metropolisHastings(lambda, t, tau, rho)
    d = length(t) - 1;
    accepted_no = zeros(1,d-1);
    
    %Metropolis-Hasting algorithm
    for i = 2:d        
        %Random walk proposal one at a time
        R = rho(1)*(t(i+1)-t(i-1));
        %create a canditate t star
        epsilon = -R + (R+R)*rand(1,1); %choose something from uniform distribution
        t_star = t(i) + epsilon;
        
        %draw until the order of the breakpoints is not changed
        while(t_star < t(i-1) || t_star >  t(i+1))
            epsilon = -R + (R+R)*rand(1,1); %choose something from uniform distribution
            t_star = t(i) + epsilon;
        end
    
        num_alpha = t_posterior(lambda, [t(1:i-1) t_star t(i+1:end)], tau);
        den_alpha = t_posterior(lambda, t, tau);
        %the probability of accepting the new state t_star
        alpha = min(1, num_alpha/den_alpha);
        
        %draw a radom u
        U = rand(1);
        %update the breakpoints
        if U <= alpha
            t(i) = t_star;
            accepted_no(i-1) = accepted_no(i-1)+1;
        end
    end
end

