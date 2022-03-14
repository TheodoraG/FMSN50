function ft = t_posterior(lambda, t, tau)
    d = length(t) - 1;
    t_diff = zeros(d, 1);

    disasters = zeros(d, 1);
    %compute the number of disasters in the given interval
    for i = 1:d
        disasters(i) = sum((t(i) <= tau) & (tau < t(i+1)));
    end

    for j = 1:d
        t_diff(j) = t(j+1) - t(j);
    end    
    %posterior f(t|lambda,tau,theta)
    ft = exp(sum(log(lambda).*disasters + log(t_diff) - lambda.*t_diff));
end

