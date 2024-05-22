addpath ..\src\
close all;

N_max = 50; % Number of samples
Exp_per_N = 100; % Number of experiments per N

%%% x[n] = h * theta + w[n]
h = 5;
m_w = 0;
sigma_w = 2;

% Theta is chosen from a normal distribution N(m_theta, sigma_theta^2)
m_theta = 4;
sigma_theta = 1;
theta = m_theta + sigma_theta * randn(1, 1);

theta_est_mse1 = mle(theta, h, m_w, sigma_w, N_max, Exp_per_N);
theta_est_mse2 = mmse(theta, m_theta, sigma_theta, h, m_w, sigma_w, N_max, Exp_per_N);

figure;
plot(1:N_max, theta_est_mse1);
hold on;
plot(1:N_max, theta_est_mse2);
xlabel('Number of Samples (N)');
ylabel('MSE of \theta Estimate');
title('MSE of \theta Estimate vs Number of Samples');
grid on;
grid minor;

legend('MMSE Estimator', 'MLE Estimator');