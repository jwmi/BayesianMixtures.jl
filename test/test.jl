

using BayesianMixtures
using Statistics
using PyPlot

n_total = 1000  # total number of MCMC sweeps to run


# _____________________________________________________________________________________
# Check samplers on univariate data

n = 50
x = randn(n)
for mode_ in ["Normal","NormalNonoptimized"]
	for model_type_ in ["MFM","DPM"]
		# Run MCMC sampler
		options = BayesianMixtures.options(mode_,model_type_,x,n_total)
		result = BayesianMixtures.run_sampler(options)
	end
end


# _____________________________________________________________________________________
# Check samplers on multivariate data

d = 3
n = 50
x = randn(n,d)
x = (x .- mean(x,dims=1))./std(x,dims=1; corrected=false)
x = [x[i,:] for i=1:n]
for mode_ in ["MVN","MVNaaC","MVNaaN","MVNaaRJ"]
	for model_type_ in ["MFM","DPM"]
		if model_type_=="DPM" && mode_=="MVNaaRJ"; continue; end
		# Run MCMC sampler
		options = BayesianMixtures.options(mode_,model_type_,x,n_total)
		result = BayesianMixtures.run_sampler(options)
	end
end


# _____________________________________________________________________________________
# Check other functions using univariate data


# Settings
mode = "Normal"
model_type = "MFM"
n = 100
B = BayesianMixtures

# Generate data and run sampler
x = randn(n)
options = B.options(mode,model_type,x,n_total)
result = B.run_sampler(options)

# Save and load
B.save_result("temp.jld",result)
result2 = B.load_result("temp.jld")

B.t_running(result)
B.t_posterior(result)
B.k_posterior(result)

B.plot_density_estimate(result)


k,a,w,theta = B.sample_mixture_parameters(result,50)
loglik = B.log_likelihoods(x,w,theta,mode)


B.open_figure(1)
B.traceplot(result.t); pause(0.5); clf()
B.traceplot_timewise(result,1); pause(0.5); clf()
B.plot_t_running(result); pause(0.5); clf()
B.plot_autocorrelation(result,20); pause(0.5); clf()
B.plot_t_posterior(result); pause(0.5); clf()
B.plot_t_posterior_average([result]); pause(0.5); clf()
B.plot_k_posterior(result); pause(0.5); clf()
B.plot_k_posterior_average([result]); pause(0.5); clf()
B.plot_similarity_matrix(result); pause(0.5); clf()
B.plot_clusters(x,result.z[:,end]); pause(0.5); clf()
B.rug_plot(x); pause(0.5); clf()
B.plot_density_estimate(result); pause(0.5); clf()



# _____________________________________________________________________________________
# Check other functions using multivariate data


# Settings
mode = "MVN"
model_type = "MFM"
d = 2
n = 100
B = BayesianMixtures

# Generate data and run sampler
x = randn(n,d)
x = (x .- mean(x,dims=1))./std(x,dims=1; corrected=false)
x = [x[i,:] for i=1:n]
options = B.options(mode,model_type,x,n_total)
result = B.run_sampler(options)

# Save and load
B.save_result("temp.jld",result)
result2 = B.load_result("temp.jld")

B.t_running(result)
B.t_posterior(result)
B.k_posterior(result)

B.plot_density_estimate(result)


k,a,w,theta = B.sample_mixture_parameters(result,50)
loglik = B.log_likelihoods(x,w,theta,mode)


B.open_figure(1)
B.traceplot(result.t); pause(0.5); clf()
B.traceplot_timewise(result,1); pause(0.5); clf()
B.plot_t_running(result); pause(0.5); clf()
B.plot_autocorrelation(result,20); pause(0.5); clf()
B.plot_t_posterior(result); pause(0.5); clf()
B.plot_t_posterior_average([result]); pause(0.5); clf()
B.plot_k_posterior(result); pause(0.5); clf()
B.plot_k_posterior_average([result]); pause(0.5); clf()
B.plot_similarity_matrix(result); pause(0.5); clf()
B.plot_clusters(x,result.z[:,end]); pause(0.5); clf()
B.plot_density_estimate(result); pause(0.5); clf()








