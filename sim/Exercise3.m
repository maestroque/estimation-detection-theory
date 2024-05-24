% Load the data from the given file. They will be stored in a variable
% called data.
load('../assets/data.mat');

data_mean = mean(data);
sum(data)
n = size(data, 2);
data_variance = var(data);

% Create a figure featuring the histogram of the data.
f3_1 = figure;
data_hist = histogram(data, 'Normalization','pdf');
hold on;
xline(data_mean, 'Color', 'r', 'LineStyle', '--');

% === Comparison Strategy ===
%
% Given that the log-likelihood function expresses how well the parameter
% explains the observed data, we will use it to determine which
% distribution was used to produce the given data. Specifically, the
% distribution with the larger log-likelihood value is more likely to 
% correspond to the given data.

% === Exponential MLE ===
%
% The function corresponding to the vectorised data form is the following:
%
%               N-1                                     N-1
%   pe(x;b) = PRODUCT b * exp(-b*x[i]) = b^N * exp(-b * SUM x[i])
%               i=0                                     i=0
%
% Therefore the log-likelihood function is of the following form:
%
%                                 N-1
%   ln(pe(x;b)) = N * ln(b) - b * SUM x[i]
%                                 i=0
%
% To find the MLE estimator, we need calculate the derivative of the log-
% likelihood function:
%
%                              N-1
%   D_b(ln(pe(x;b))) = N / b - SUM x[i]
%                              i=0
%                                   N-1                 N-1
%   D_b(ln(pe(x;b))) = 0 => N / b = SUM x[i] => b = N / SUM x[i]
%                                   i=0                 i=0
%   => b = 1 / E[x]
%
% Or in code:
exp_estimator = 1 / data_mean;
exp_log_likelihood = @(lambda, x) size(x, 2) * log(lambda) - lambda * sum(x);
% === Rayleigh MLE ===
%
% The function corresponding to the vectorised data form is the following:
%
%                         N-1                               N-1
%   pr(x;b) = b^(-2N) * PRODUCT(x[i]) * exp(-0.5 * b^(-2) * SUM(x[i]^2))
%                         i=0                               i=0
%
% Therefore the log-likelihood function is of the following form:
%
%                 N-1                                         N-1
%   ln(pr(x;b)) = SUM(ln(x)) - 2 * N * ln(b) - 0.5 * b^(-2) * SUM(x[i]^2)
%                 i=0                                         i=0
%
% To find the MLE estimator, we need calculate the derivative of the log-
% likelihood function:
%
%                                           N-1
%   D_b(ln(pe(x;b))) = 2 * N / b + b^(-3) * SUM(x[i]^2)
%                                           i=0
%                                                N-1
%   D_b(ln(pe(x;b))) = 0 => 2 * N / b = b^(-3) * SUM(x[i]^2)
%                                                i=0
%                               N-1
%   => b = sqrt((1 / (2 * N)) * SUM(x[i]^2))
%                               i=0
% Or in code:
rayleigh_estimator = sqrt(0.5 * sum(data.^2) / n);
rayleigh_log_likelihood = @(lambda, x) sum(log(x)) - 2 * size(x, 2) * log(lambda) - sum(data.^2) / (2 * lambda^2);

exp_log_likelihood_est = exp_log_likelihood(exp_estimator, data);
rayleigh_log_likelihood_est = rayleigh_log_likelihood(rayleigh_estimator, data);

% Print the results
fprintf("Exponential log likelihood after estimation: %f\n", exp_log_likelihood_est);
fprintf("Rayleigh log likelihood after estimation: %f\n", rayleigh_log_likelihood_est);

% Plot the pdfs of the two distributions to show which one better fits the
% dataset.
horizontal_axis = 0:0.01:max(data);
exp_pdf_curve = exppdf(horizontal_axis, 1 / exp_estimator);
rayleigh_pdf_curve = raylpdf(horizontal_axis, rayleigh_estimator);
plot(horizontal_axis, exp_pdf_curve, 'Color', 'g');
hold on;
plot(horizontal_axis, rayleigh_pdf_curve, 'Color', 'm');
xlabel('Data value');
ylabel('Value of pdf');
title('Histogram of the data');
grid on;
grid minor;
legend('Data', 'Mean', 'Exponential Distribution', 'Rayleigh Distribution');
saveas(f3_1, 'media/ex3_data_histogram.jpg');
