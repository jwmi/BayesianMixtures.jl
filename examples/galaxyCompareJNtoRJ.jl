# Compare Jain-Neal sampler to Reversible Jump on the galaxy dataset.
module GalaxyCompareJNtoRJ

using PyPlot
using DelimitedFiles
using BayesianMixtures
B = BayesianMixtures

using Pkg
packages_installed = [pkg.name for pkg in collect(values(Pkg.dependencies()))]
can_plot = ("PyPlot" in packages_installed)
can_plot ? using PyPlot : @warn("Skipping plots since PyPlot is not installed.")

# ==================== Run Jain-Neal sampler ====================

# Read the data from file
x = readdlm("datasets/galaxy.dat",' ',Float64)[:]

# Specify model, data, and MCMC options
n_total = 1000000  # total number of MCMC sweeps to run
log_pk = "k -> log(k in (1:30) ? 1/30 : 0)"   # log prior on k is Uniform{1,...,30}
options = B.options("Normal","MFM",x,n_total; n_keep=1000,log_pk=log_pk,t_max=30)

# Run MCMC sampler
result = B.run_sampler(options)

# ==================== Run RJMCMC sampler using Nmix ====================

# Remove old output files
if (!isdir("galx")); mkdir("galx"); end
outfiles = ["1.bdlog","1.bk","1.ent","1.k","1.log","1.out","1.pe"]
for outfile in outfiles; try rm(joinpath("galx",outfile)); catch; end; end

# Run Peter Green's Nmix program
if Sys.isapple() # if running Mac OSX
    run(`chmod +x Nmix`)
    elapsed_time_RJ = (@elapsed run(`./Nmix -n$n_total -nb0 -ns$n_total galx`))
elseif Sys.iswindows() # if running Windows
    elapsed_time_RJ = (@elapsed run(`Nmix.exe -n$n_total -nb0 -ns$n_total galx`))
else
    error("Only Mac and Windows are currently supported for running Green's Nmix program.")
end

# Read and delete the output
output = readdlm(joinpath("galx","1.out"))
for outfile in outfiles; rm(joinpath("galx",outfile)); end

# Parse the output
k_RJ = zeros(Int64,n_total)
t_RJ = zeros(Int64,n_total)
i = 1
for j = 1:n_total
	global i
    k,t = output[i,findall(output[i,:].!="")]
    k_RJ[j] = k
    t_RJ[j] = t
    i += k+1
end

# ==================== Compare posterior on k ====================

# Compute the posterior on k from Jain-Neal sampler
k_posterior_JN = B.k_posterior(result)

# Compute posterior on k from RJMCMC sampler
use = options.n_burn+1:n_total
k_posterior_RJ = B.histogram(k_RJ[use], 0:options.t_max)[1]/length(use)

# Print posteriors on k
println("Posteriors on k:")
println("  Jain-Neal  RJMCMC")
using Printf
for k=1:15
    @printf "k=%d: %.4f %.4f\n" k k_posterior_JN[k] k_posterior_RJ[k]
end

if can_plot
    
    # ==================== Compare traceplots of t ====================
    # Traceplot of t for Jain-Neal
    B.open_figure(1)
    B.traceplot(result.t[1:min(5000,n_total)])
    B.labels("Traceplot for t (Jain-Neal)","sweep number","t")

    # Traceplot of t for RJMCMC
    B.open_figure(2)
    B.traceplot(t_RJ[1:min(5000,n_total)])
    B.labels("Traceplot for t (RJMCMC)","sweep number","t")

    # ==================== Compare autocorrelation of t ====================
    # Plot autocorrelation for Jain-Neal
    B.open_figure(3)
    B.plot_autocorrelation(result,0.01)

    # Plot autocorrelation for RJMCMC
    B.plot_autocorrelation(t_RJ[use],0.01,elapsed_time_RJ/n_total; color="r")
    legend(["Jain-Neal","RJMCMC"],loc="upper right",fontsize=12)
   
    # ==================== Plot density estimate ====================
    # Plot density estimate using Jain-Neal output
    B.open_figure(4)
    hist(x; color="w",edgecolor="g",linewidth=0.5)
    B.rug_plot(x)
    B.plot_density_estimate(result; resolution=500)
end

end # module GalaxyCompareJNtoRJ



