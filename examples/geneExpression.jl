# Run collapsed Jain-Neal sampler on a gene expression data set.
module GeneExpression

using BayesianMixtures
B = BayesianMixtures
can_plot = (Pkg.installed("PyPlot")!=nothing)
can_plot ? using PyPlot : warn("Skipping plots since PyPlot is not installed.")

# Load data
collection = "deSouto"
data_ID = "armstrong-2002-v1"
# collection = "Siri"
# data_ID = "Blood1"
if collection=="deSouto"
    dataset = readdlm("datasets/gene-deSouto/$(data_ID)_database.txt",'\t')
    x = convert(Array{Float64,2}, dataset[3:end,2:end])
    x = log2(x)
    x = [x[:,i]::Array{Float64,1} for i = 1:size(x,2)]
elseif collection=="Siri"
    dataset = readdlm("datasets/gene-Siri/$(data_ID).csv",',')
    x = convert(Array{Float64,2}, dataset[2:end,2:end])
    x = [x[:,i]::Array{Float64,1} for i = 1:size(x,2)]
else
    error("Unknown collection: $collection")
end
# Normalize to zero mean, unit variance
mu = mean(x)
v = mean([xi.*xi for xi in x]) - mu.*mu  # sample variance
x = [((xi-mu)./sqrt(v))::Array{Float64,1} for xi in x]
    

# Run sampler
n_total = 1000  # number of MCMC sweeps to run the algorithm
options = B.options("MVNaaC","MFM",x,n_total; n_keep=min(1000,n_total),t_max=20)
result = B.run_sampler(options)


# Posteriors on t and k
B.open_figure(1)
B.plot_t_posterior(result; color="b", marker="s", label="p(t|data)")
B.plot_k_posterior(result; color="g", marker="o", label="p(k|data)")
B.labels("Posteriors on t and k","t (clusters) or k (components)","")
legend(loc="upper right",numpoints=1)
xlim(0,15)

# Posterior similarity matrix (probability that i and j are in same cluster)
B.open_figure(2; figure_size=(5,4))
B.plot_similarity_matrix(result)


end

