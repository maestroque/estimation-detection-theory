% Detection and Estimation Theory Project - Exercise 1

addpath ..\src\
close all;

% 1. Find CRLB and for a set number of the parameter, run 5000 experiments.

% For an exponential distribution, the CRLB can be calculated as follows:
%               N-1                                         N-1
%   pe(x; b) = PRODUCT(b * exp(-b * x[i])) = b^N * exp(-b * SUM(x[i]))
%               i=0                                         i=0
%                                  N-1
%   ln(pe(x; b)) = N * ln(b) - b * SUM(x[i])
%                                  i=0
%                               N-1
%   D_b(ln(pe(x; b))) = N / b - SUM(x[i])
%                               i=0
%
%   D2_b(ln(pe(x; b))) = - N / b^2
%
% The fisher information is equal to:
%
%   I(b) = -E[D2_b(ln(pe(x; b)))] = -E[-N / b^2] = N / b^2
%
% Therefore: var(b) >= 1 / I(b) => var(b) >= b^2 / N
%
% Or, in code:
crlb_for_exponential_dist = @(lambda, n) lambda.^2 / n;

% Set the parameter, the number of experiments and the samples
true_lambda = 2;
number_of_experiments = 5000;
number_of_samples = 10;

% An ideal scenario will be simulated for this question, using a Normal
% distribution. The estimator in this case, is the mean of the data (shown 
% in in-class theory.).
sigma_sq = crlb_for_exponential_dist(true_lambda, number_of_samples);
ideal_dist = makedist('Normal', true_lambda, sqrt(number_of_samples * sigma_sq));
ideal_pdf = @(n) random(ideal_dist, [n, 1]);
ideal_estimator = @(samples) mean(samples);

ideal_estimation = experiment(ideal_pdf, ideal_estimator, number_of_samples, number_of_experiments);
fprintf("Ideal estimation: mean: %f    variance: %f\n", mean(ideal_estimation), var(ideal_estimation));
f1_1 = figure;
ideal_histogram = histogram(ideal_estimation, 'Normalization', 'pdf');
hold on;
xline(true_lambda, 'Color', 'r', 'LineStyle', '--');
hold on;
xline(mean(ideal_estimation), 'Color', 'y', 'LineStyle', '--');
xlabel('Estimate values');
ylabel('Estimate pdf value');
title('Ex.1.1 - Ideal unbiased estimation histogram - 10 samples per experiment');
grid on;
grid minor;
legend('Estimation results', 'True mean', 'Estimation mean');
saveas(f1_1, 'media/ex1_ideal_10_samples.jpg');

% 2. Find the MLE and run the experiment again.

% Using the formulas described in the first part of the exercise, the MLE
% estimator q can be found:
%                                    N-1                  N-1
%   D_b(ln(pe(x; q))) = 0 => N / q = SUM(x[i]) => q = N / SUM(x[i])
%                                    i=0                  i=0
%   => q = 1 / E[x]
exp_dist = makedist('Exponential', 1 / true_lambda);
exp_pdf = @(n) random(exp_dist, [n, 1]);
mle_estimator = @(samples) 1 / mean(samples);
mle_estimation = experiment(exp_pdf, mle_estimator, number_of_samples, number_of_experiments);

fprintf("MLE Estimation: mean: %f    variance: %f\n", mean(mle_estimation), var(mle_estimation));
f1_2 = figure;
mle_histogram = histogram(mle_estimation, 'Normalization', 'pdf');
hold on;
xline(true_lambda, 'Color', 'r', 'LineStyle', '--');
hold on;
xline(mean(mle_estimation), 'Color', 'y', 'LineStyle', '--');
xlabel('MLE estimate values');
ylabel('MLE estimate pdf value');
title('Ex.1.2 - MLE estimation histogram - 10 samples per experiment');
grid on;
grid minor;
legend('Estimation results', 'True mean', 'Estimation mean');
saveas(f1_2, 'media/ex1_mle_10_samples.jpg');

% 3. Repeat the previous procedures for 20, 50 and 100 samples per
% experiment

samples_n = [20 50 100];
squared_sigma = [ 0 0 0 ];
for i=1:3
    squared_sigma(i) = crlb_for_exponential_dist(true_lambda, samples_n(i));
end

% Ideal estimator - 20 samples per experiment
ideal_dist_20 = makedist('Normal', true_lambda, sqrt(samples_n(1) * squared_sigma(1)));
ideal_pdf_20 = @(n) random(ideal_dist_20, [n, 1]);
ideal_estimation_20 = experiment(ideal_pdf_20, ideal_estimator, samples_n(1), number_of_experiments);
fprintf("Ideal estimation (20 samples): mean: %f    var: %f\n", mean(ideal_estimation_20), var(ideal_estimation_20));
f1_3_1 = figure;
ideal_histogram_20 = histogram(ideal_estimation_20, 'Normalization', 'pdf');
hold on;
xline(true_lambda, 'Color', 'r', 'LineStyle', '--');
hold on;
xline(mean(ideal_estimation_20), 'Color', 'y', 'LineStyle', '--');
xlabel('Estimate values');
ylabel('Estimate pdf value');
title('Ex.1.3.1 - Ideal unbiased estimation histogram - 20 samples per experiment');
grid on;
grid minor;
legend('Estimation results', 'True mean', 'Estimation mean');
saveas(f1_3_1, 'media/ex1_ideal_20_samples.jpg');

% Ideal estimator - 50 samples per experiment
ideal_dist_50 = makedist('Normal', true_lambda, sqrt(samples_n(2) * squared_sigma(2)));
ideal_pdf_50 = @(n) random(ideal_dist_50, [n, 1]);
ideal_estimation_50 = experiment(ideal_pdf_50, ideal_estimator, samples_n(2), number_of_experiments);
fprintf("Ideal estimation (50 samples): mean: %f    var: %f\n", mean(ideal_estimation_50), var(ideal_estimation_50));
f1_3_2 = figure;
ideal_histogram_50 = histogram(ideal_estimation_50, 'Normalization', 'pdf');
hold on;
xline(true_lambda, 'Color', 'r', 'LineStyle', '--');
hold on;
xline(mean(ideal_estimation_50), 'Color', 'y', 'LineStyle', '--');
xlabel('Estimate values');
ylabel('Estimate pdf value');
title('Ex.1.3.2 - Ideal unbiased estimation histogram - 50 samples per experiment');
grid on;
grid minor;
legend('Estimation results', 'True mean', 'Estimation mean');
saveas(f1_3_2, 'media/ex1_ideal_50_samples.jpg');

% Ideal estimator - 100 samples per experiment
ideal_dist_100 = makedist('Normal', true_lambda, sqrt(samples_n(3) * squared_sigma(3)));
ideal_pdf_100 = @(n) random(ideal_dist_100, [n, 1]);
ideal_estimation_100 = experiment(ideal_pdf_100, ideal_estimator, samples_n(3), number_of_experiments);
fprintf("Ideal estimation (100 samples): mean: %f    var: %f\n", mean(ideal_estimation_100), var(ideal_estimation_100));
f1_3_3 = figure;
ideal_histogram_100 = histogram(ideal_estimation_100, 'Normalization', 'pdf');
hold on;
xline(true_lambda, 'Color', 'r', 'LineStyle', '--');
hold on;
xline(mean(ideal_estimation_100), 'Color', 'y', 'LineStyle', '--');
xlabel('Estimate values');
ylabel('Estimate pdf value');
title('Ex.1.3.3 - Ideal unbiased estimation histogram - 100 samples per experiment');
grid on;
grid minor;
legend('Estimation results', 'True mean', 'Estimation mean');
saveas(f1_3_3, 'media/ex1_ideal_100_samples.jpg');

% MLE estimator - 20 samples per experiment
mle_estimation_20 = experiment(exp_pdf, mle_estimator, samples_n(1), number_of_experiments);
fprintf("MLE estimation (20 samples): mean: %f    var: %f\n", mean(mle_estimation_20), var(mle_estimation_20));
f1_3_4 = figure;
mle_histogram_20 = histogram(mle_estimation_20, 'Normalization', 'pdf');
hold on;
xline(true_lambda, 'Color', 'r', 'LineStyle', '--');
hold on;
xline(mean(mle_estimation_20), 'Color', 'y', 'LineStyle', '--');
xlabel('MLE estimate values');
ylabel('MLE estimate pdf value');
title('Ex.1.3.4 - MLE estimation histogram - 20 samples per experiment');
grid on;
grid minor;
legend('Estimation results', 'True mean', 'Estimation mean');
saveas(f1_3_4, 'media/ex1_mle_20_samples.jpg');

% MLE estimator - 50 samples per experiment
mle_estimation_50 = experiment(exp_pdf, mle_estimator, samples_n(2), number_of_experiments);
fprintf("MLE estimation (50 samples): mean: %f    var: %f\n", mean(mle_estimation_50), var(mle_estimation_50));
f1_3_5 = figure;
mle_histogram_50 = histogram(mle_estimation_50, 'Normalization', 'pdf');
hold on;
xline(true_lambda, 'Color', 'r', 'LineStyle', '--');
hold on;
xline(mean(mle_estimation_50), 'Color', 'y', 'LineStyle', '--');
xlabel('MLE estimate values');
ylabel('MLE estimate pdf value');
title('Ex.1.3.5 - MLE estimation histogram - 50 samples per experiment');
grid on;
grid minor;
legend('Estimation results', 'True mean', 'Estimation mean');
saveas(f1_3_5, 'media/ex1_mle_50_samples.jpg');

% MLE estimator - 100 samples per experiment
mle_estimation_100 = experiment(exp_pdf, mle_estimator, samples_n(3), number_of_experiments);
fprintf("MLE estimation (100 samples): mean: %f    var: %f\n", mean(mle_estimation_100), var(mle_estimation_100));
f1_3_6 = figure;
mle_histogram_100 = histogram(mle_estimation_100, 'Normalization', 'pdf');
hold on;
xline(true_lambda, 'Color', 'r', 'LineStyle', '--');
hold on;
xline(mean(mle_estimation_100), 'Color', 'y', 'LineStyle', '--');
xlabel('MLE estimate values');
ylabel('MLE estimate pdf value');
title('Ex.1.3.6 - MLE estimation histogram - 100 samples per experiment');
grid on;
grid minor;
legend('Estimation results', 'True mean', 'Estimation mean');
saveas(f1_3_6, 'media/ex1_mle_100_samples.jpg');

