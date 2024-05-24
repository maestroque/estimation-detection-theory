addpath ..\src\
close all;

N_max = 50; % Number of samples
Exp_per_N = 10; % Number of experiments per N

%%% x[n] = h * theta + w[n]
h = 0.5;
m_w = 0;
sigma_w = 2;

% Theta is chosen from a normal distribution N(m_theta, sigma_theta^2)
m_theta = 4;
sigma_theta = 1;
theta_norm = m_theta + sigma_theta * randn(1, 1);

theta_est_mse_mle = mle_implementation(theta_norm, h, m_w, sigma_w, N_max, Exp_per_N);
theta_est_mse_mmse = mmse_implementation(theta_norm, m_theta, sigma_theta, h, m_w, sigma_w, N_max, Exp_per_N);

theta_not_norm = 8;

theta_est_mse_mle_not_norm = mle_implementation(theta_not_norm, h, m_w, sigma_w, N_max, Exp_per_N);
theta_est_mse_mmse_not_norm = mmse_implementation(theta_not_norm, m_theta, sigma_theta, h, m_w, sigma_w, N_max, Exp_per_N);


exp_dist = makedist('Exponential', 4);
uniform_dist = makedist('Uniform', 0, 8);

theta_exp = random(exp_dist, 1);
theta_uniform = random(uniform_dist, 1);

theta_est_mse_exp = mmse_implementation(theta_exp, m_theta, sigma_theta, h, m_w, sigma_w, N_max, Exp_per_N);
theta_est_mse_uniform = mmse_implementation(theta_uniform, m_theta, sigma_theta, h, m_w, sigma_w, N_max, Exp_per_N);

f2_3_1 = figure;
plot(1:N_max, theta_est_mse_mle);
hold on;
plot(1:N_max, theta_est_mse_mmse);
xlabel('Number of Samples (N)');
ylabel('MSE of \theta Estimate');
title('MSE of \theta (from N(4,1)) Estimate vs Number of Samples');
grid on;
grid minor;
legend('MLE Estimator', 'MMSE Estimator');
saveas(f2_3_1, 'media/ex2_mle_mmse_comparison_norm.jpg');

f2_3_2 = figure;
plot(1:N_max, theta_est_mse_mle_not_norm);
hold on;
plot(1:N_max, theta_est_mse_mmse_not_norm);
xlabel('Number of Samples (N)');
ylabel('MSE of \theta Estimate');
title('MSE of \theta (not from N(4, 1)) Estimate vs Number of Samples');
grid on;
grid minor;
legend('MLE Estimator', 'MMSE Estimator');
saveas(f2_3_2, 'media/ex2_mle_mmse_comparison_not_norm.jpg');

f2_4 = figure;
plot(1:N_max, theta_est_mse_mmse);
hold on;
plot(1:N_max, theta_est_mse_exp);
hold on;
plot(1:N_max, theta_est_mse_uniform);
xlabel('Number of Samples (N)');
ylabel('MSE of \theta Estimate');
title('MSE of \theta Estimate vs Number of Samples', 'For Different Priors and MMSE Estimator for N(4, 1)');
grid on;
grid minor;
legend('Gaussian Prior', 'Exponential Prior', 'Uniform Prior');
saveas(f2_4, 'media/ex2_different_priors.jpg');