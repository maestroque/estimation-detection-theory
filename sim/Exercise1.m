addpath ..\src\
close all;

crlb_for_exponential_dist = @(lambda, n) lambda.^2 / n;

true_lambda = 2;
number_of_experiments = 5000;
number_of_samples = 10;

sd = crlb_for_exponential_dist(true_lambda, number_of_samples);
virtual_dist = makedist('Normal', true_lambda, sqrt(sd));
virtual_pdf = @(n) random(virtual_dist, [number_of_samples, 1]);
virtual_estimator = @(samples) mean(samples);

virtual_estimation = experiment(virtual_pdf, virtual_estimator, number_of_samples, number_of_experiments);
h_1_1 = figure;
virtual_histogram = histogram(virtual_estimation);
hold on;
xline(true_lambda, 'Color', 'r');

exp_dist = makedist('Exponential', 1 / true_lambda);
exp_pdf = @(n) random(exp_dist, [number_of_samples, 1]);
mle_estimator = @(samples) 1 / mean(samples);
mle_estimation = experiment(exp_pdf, mle_estimator, number_of_samples, number_of_experiments);

h_1_2 = figure;
mle_histogram = histogram(mle_estimation);
hold on;
xline(true_lambda, 'Color', 'r')
