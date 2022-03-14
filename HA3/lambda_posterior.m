function lambda_post = lambda_posterior(theta, t, tau)
    t_diff = t(2:end) - t(1:end-1);
    d = length(t) - 1;
    disasters = zeros(1, d);
    %compute the number of disasters in the given interval
    for i = 1:d
        disasters(i) = sum((t(i) <= tau) & (tau < t(i+1)));
    end
    %posterior f(lambda|tau,t,theta)
    lambda_post = gamrnd(2 + disasters', 1./(t_diff' + theta));
end

