% mmse - Calculates the Mean Squared Error (MSE) of an MLE Estimator
% for a given deterministic parameter theta, for the following signal model:
%   x[n] = h * theta + w[n]
%
% Other parameters:
%   N_max: Maximum number of samples to be used for estimation
%   Exp_per_N: Number of experiments to be performed for each N
%
% Output:
%   - estimator_mse: Array of MSE values for different sample sizes (N)
function estimator_mse = mle(theta, h, m_w, sigma_w, N_max, Exp_per_N)
    estimator_mse = zeros(1, N_max);

    for i = 1:N_max
        theta_est = zeros(1, Exp_per_N);
        for j = 1:Exp_per_N
            N = i;
            % Generate random samples
            w = m_w + sigma_w * randn(N, 1);
            x = h * theta + w;
            theta_est(j) = 1 / (N * h^2) * sum(h * x);
        end
        estimator_mse(i) = mean((theta_est - theta).^2);
    end
end