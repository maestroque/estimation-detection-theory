function estimation = experiment(data_pdf, estimator, n_samples, n_experiments)
% EXPERIMENT Runs an estimation experiment on generated samples.
%   ESTIMATION = EXPERIMENT(pdf, estimator, n_samples, n_experiments)
%       generates n_samples samples using pdf and calculates an estimation
%       using estimator. This is repeated n_experiment times. In the end,
%       the estimation vector contains the estimations of each experiment.
    estimation = zeros(n_experiments, 1);
    for i = 1:n_experiments
        samples = feval(data_pdf, n_samples);
        estimation(i) = feval(estimator, samples);
    end
end

