function estimation = experiment(data_pdf, estimator, n_samples, n_experiments)
    estimation = zeros(n_experiments, 1);
    for i = 1:n_experiments
        samples = feval(data_pdf, n_samples);
        estimation(i) = feval(estimator, samples);
    end
end

