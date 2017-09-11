# Main module for BayesianMixtures package
module BayesianMixtures

include("MFM.jl")
include("Random.jl")

include("Normal.jl")
include("MVN.jl")
include("MVNaaC.jl")
include("MVNaaN.jl")
include("MVNaaRJ.jl")
include("NormalNonoptimized.jl")

# ===================================================================
# ===================================================================
# ================== Functions to generate results ==================
# ===================================================================
# ===================================================================

# Create an options object to specify model, data, and MCMC parameters.
function options(
        mode, # "Normal", "MVN", "MVNaaC", "MVNaaN", or "MVNaaRJ" 
        model_type, # "MFM" or "DPM"
        x, # data
        n_total; # total number of MCMC sweeps to run the sampler
        n_keep=n_total, # number of MCMC sweeps to keep after thinning
        n_burn=round(Int,n_total/10), # number of MCMC sweeps (out of n_total) to discard as burn-in
        verbose=true, # display information or not
        use_hyperprior=true, # put hyperprior on the base distn parameters or not
        t_max=40, # a guess at an upper bound on # of clusters that will be encountered during MCMC

        # MFM options:
        gamma=1.0, # Dirichlet_k(gamma,...,gamma)
        log_pk="k -> log(0.1)+(k-1)*log(0.9)", # string representation of log(p(k))
            # (where p(k) is the log of the prior on # of components, K)

        # DPM options:
        alpha_random=true, # put prior on alpha (DPM concentration parameter) or not
        alpha=1.0, # value of alpha (initial value if alpha_random=true)

        # Jain-Neal split-merge options:
        use_splitmerge=true, # use split-merge or not
        n_split=5, # number of intermediate sweeps for split launch state
        n_merge=5,  #                 "         "       merge    "     "  
        
        # RJMCMC options:
        k_max=t_max # a guess at an upper bound on # of components that will be encountered during MCMC
    )

    # Compute partition distribution values
    n = length(x)
    if model_type=="MFM"
        lpk = eval(parse(log_pk))
        log_pk_fn(k) = Base.invokelatest(lpk,k)
        log_v = MFM.coefficients(log_pk_fn,gamma,n,t_max+1)
        a = b = gamma
    elseif model_type=="DPM"
        log_v = float(1:t_max+1)*log(alpha) - lgamma(alpha+n) + lgamma(alpha)
        a,b = 1.,0.
    else
        error("Invalid model_type: $model_type.")
    end
    if mode=="MVNaaRJ"; @assert(model_type=="MFM", "RJMCMC is not implemented for DPMs."); end

    n_keep = min(n_keep,n_total)
    module_ = getfield(BayesianMixtures,Symbol(mode))
    return module_.Options(mode, model_type, x, n_total, n_keep, n_burn, verbose,
                           use_hyperprior, t_max, gamma, log_pk, alpha_random, alpha,
                           use_splitmerge, n_split, n_merge, k_max, a, b, log_v, n)
end


# Run the MCMC sampler with the specified options.
function run_sampler(options)
    o = options
    n,n_total,n_keep = o.n,o.n_total,o.n_keep
    module_ = getfield(BayesianMixtures,Symbol(o.mode))

    # Short run to precompile (calling precompile doesn't seem to work...)
    module_.sampler(o,1,1)

    if o.verbose
        println(o.mode, " ", o.model_type)
        println("n = $n, n_total = $n_total, n_keep = $n_keep")
        print("Running... ")
    end
    
    # Main run
    tic()
    t_r,N_r,z_r,theta_r,keepers = module_.sampler(o,n_total,n_keep)
    elapsed_time = toq()
    time_per_step = elapsed_time/(n_total*n)

    if o.verbose
        println("complete.")
        println("Elapsed time = $elapsed_time seconds")
        println("Time per step ~ $time_per_step seconds")
    end

    # # Run profiler
    # @profile sampler(x,n_total,n_keep,o)
    # Profile.print(format = :flat)
    # Profile.clear()
    
    return module_.Result(o,t_r,N_r,z_r,theta_r,keepers,elapsed_time,time_per_step)
end


# ===================================================================
# ===================================================================
# ================== Functions to analyze results ===================
# ===================================================================
# ===================================================================

can_save = (Pkg.installed("JLD")!=nothing)
if can_save; using JLD; end

# Save result to file.
function save_result(filename,result)
    @assert(can_save, "Results cannot be saved to file since JLD is not installed.")
    JLD.save(filename,"result",result)
end

# Load result from file.
function load_result(filename)
    @assert(can_save, "Results cannot be loaded from file since JLD is not installed.")
    return JLD.load(filename,"result")
end

# Compute running average of posterior CDF on t (# of clusters)
function t_running(result)
    o = result.options
    n_total = o.n_total
    T = zeros(n_total,o.t_max)
    for i = 1:n_total; T[i,result.t[i]] = 1; end
    cdfs = cumsum(cumsum(T,1),2)./(1:n_total)
    return cdfs
end

# Compute histogram with the specified bin edges,
# where x[i] is in bin j if edges[j] < x[i] <= edges[j+1].
function histogram(x, edges=[]; n_bins=50)
    if isempty(edges)
        mn,mx = minimum(x),maximum(x)
        r = mx-mn
        edges = linspace(mn-r/n_bins, mx+r/n_bins, n_bins+1)
    else
        n_bins = length(edges)-1
    end
    counts = zeros(Int,n_bins)
    for i=1:length(x)
        for j=1:n_bins
            if (edges[j] < x[i] <= edges[j+1])
                counts[j] += 1
                break
            end
        end
    end
    return counts,edges
end

# Compute posterior on t (# of clusters)
function t_posterior(result)
    o = result.options
    use = o.n_burn+1:o.n_total
    counts,edges = histogram(result.t[use], 0:o.t_max) # bins are: (0,1],(1,2],...,(t_max-1,t_max]
    return counts/length(use)
end

# Compute posterior on k (# of components)
function k_posterior(result; upto=result.options.t_max)
    o = result.options
    @assert(o.model_type=="MFM", "The posterior on k is not defined for the DPM.")
    lpk = eval(parse(o.log_pk))
    log_pk_fn(k) = Base.invokelatest(lpk,k)
    p_kt = MFM.p_kt(log_pk_fn,o.gamma,o.n,upto,o.t_max)
    return p_kt*t_posterior(result)
end

# Compute density estimate at a point x (using Green's pseudo-predictive density).
function density_estimate(x,result)
    o = result.options
    n,b = o.n,o.b
    t,N,theta,keepers = result.t, result.N, result.theta, result.keepers
    n_keep = length(keepers)
    loglik = getfield(BayesianMixtures,Symbol(o.mode)).log_likelihood
    r = 0.0
    n_used = 0
    for ik = 1:n_keep
        i = keepers[ik]
        if i>o.n_burn
            for c in find(N[:,i].>0)
                r += exp(loglik(x,theta[c,ik])) * (N[c,i]+b)/(n + b*t[i])
            end
            n_used += 1
        end
    end
    return r/n_used
end

rand_categorical(p) = (s=rand(); i=0; while s>0; i+=1; s-=p[i]; end; i)
rand_dirichlet(a) = (w=Random.gamma.(a,1); w/sum(w))
counts(a,k) = (N=zeros(Int,k); for ai in a; N[ai]+=1; end; N)
logsumexp(x) = (m = maximum(x); m == -Inf ? -Inf : log.(sum(exp.(x-m))) + m)

# Convert MFM result into posterior samples of k,w,a,theta,
# where k = number of components
#       w = mixture weights
#       a = assignments of data points to mixture components
#       theta = mixture component parameters
function sample_mixture_parameters(result,kmax)
    r,o = result,result.options
    @assert(o.model_type=="MFM", "This function is not defined for the DPM.")
    lpk = eval(parse(o.log_pk))
    log_pk_fn(k) = Base.invokelatest(lpk,k)
    p_kt = MFM.p_kt(log_pk_fn,o.gamma,o.n,kmax,o.t_max)

    module_ = getfield(BayesianMixtures,Symbol(o.mode))
    Theta = module_.Theta
    H = module_.construct_hyperparameters(o)

    k = zeros(Int16,o.n_keep)
    a = zeros(Int16,o.n,o.n_keep)
    w = Array{Float64,1}[]
    theta = Array{Theta,1}[]
    for i = 1:o.n_keep
        # sample k from p(k|t)
        ikeep = r.keepers[i]
        k[i] = rand_categorical(p_kt[:,r.t[ikeep]])

        # sample a from p(a|k,C) where C = clustering
        uz = sort(unique(r.z[:,i]))
        assign = zeros(Int,o.t_max+3)
        assign[uz] = randperm(k[i])[1:length(uz)]
        a[:,i] = assign[r.z[:,i]]

        # sample w from p(w|k,a)
        push!(w, rand_dirichlet(o.gamma + counts(a[:,i],k[i])))

        # sample from p(theta|k,a,phi)
        theta_i = [module_.new_theta(H)::Theta for j=1:k[i]]
        for j = 1:k[i]
            c_j = findfirst(assign,j)
            if c_j==0
                module_.prior_sample!(theta_i[j],H)
            else
                theta_i[j] = r.theta[c_j,i]
            end
        end
        push!(theta,theta_i)
    end
    return k,a,w,theta
end

# Compute deviances
function deviances(result,w,theta)
    r,o = result,result.options
    module_ = getfield(BayesianMixtures,Symbol(o.mode))
    log_f = module_.log_likelihood
    m = length(w)
    x,n = o.x,o.n
    dev = zeros(m)
    for iter = 1:m
        k = length(w[iter])
        logL = 0.0
        for i = 1:n
            l = Float64[log(w[iter][j]) + log_f(x[i],theta[iter][j]) for j = 1:k]
            logL += logsumexp(l)
        end
        dev[iter] = -2*logL
    end
    return dev
end

# Compute the posterior similarity matrix (probability that i and j are in same cluster).
function similarity_matrix(result)
    o = result.options
    n = o.n
    z = result.z
    C = zeros(Int,n,n)
    n_used = 0
    for index = 1:o.n_keep
        if result.keepers[index]>o.n_burn
            for i = 1:n, j = 1:n
                if z[i,index]==z[j,index]; C[i,j] += 1; end
            end
            n_used += 1
        end
    end
    return C/n_used
end

# Compute autocorrelation curve for the sequence y=[y_1,...,y_T]
function autocorrelation(y)
    T = length(y)
    u = ones(typeof(y[1]),T)
    d = [1:T; T-1:-1:1]
    autocov = xcorr(y,y)./d - (xcorr(y,u)./d).*(xcorr(u,y)./d)
    a = autocov[T:-1:1]
    return a/a[1]
end
#autocorrelation_time(auto,c) = (tau=1+2*cumsum(auto[2:end]); M=findfirst([1:length(tau)] .>= c*tau); tau[M])


# ===================================================================
# ===================================================================
# ==================== Functions to plot results ====================
# ===================================================================
# ===================================================================
    
can_plot = (Pkg.installed("PyPlot")!=nothing)
if can_plot; using PyPlot; end
checkplotting() = (if !can_plot; error("Plotting is disabled since PyPlot is not installed."); end)

# ---------- General plotting functions ----------

function draw_now(number=0)
    checkplotting()
    if number>0; figure(number); end
    pause(0.001) # this forces the figure to update (draw() is supposed to do this, but doesn't work for me)
    get_current_fig_manager()[:window][:raise_]() # bring figure window to the front
    #w=get_current_fig_manager()[:window][:attributes]; w("-topmost",1); w("-topmost", 0)
end

function open_figure(number;clear_figure=true,figure_size=(5,2.5))
    checkplotting()
    fig = figure(number, figsize=figure_size)
    subplots_adjust(top=0.85, bottom=0.2, left=0.15)
    clear_figure? clf() : nothing
    # get_current_fig_manager()[:window][:showMaximized]() # maximize figure window
    return fig
end
function labels(title_string,xlabel_string="",ylabel_string="")
    checkplotting()
    title(title_string,fontsize=14)
    xlabel(xlabel_string,fontsize=14)
    ylabel(ylabel_string,fontsize=14)
end
function density_grid(x0,x1,y0,y1,xres,yres)
    xs = linspace(x0,x1,xres)
    ys = linspace(y0,y1,yres)
    X = [x for x=xs, y=ys]
    Y = [y for x=xs, y=ys]
    return xs,ys,X,Y
end

# ---------- Plotting results ----------
color_list = "cgyrbmkw"^1000
marker_list = "ov^Ds<>p*+x1234hHd"

function traceplot(values)
    checkplotting()
    n_values = length(values)
    n_subset = min(2000,n_values)
    if (n_subset < n_values); warn("Traceplot shows only a subset of points."); end
    subset = round.(Int,linspace(1,n_values,n_subset))
    jitter = (rand(n_subset)-0.5)/2
    PyPlot.plot(subset,values[subset]+jitter, "k.", markersize=1.0)
    draw_now()
end

function traceplot_timewise(result,t_show)
    checkplotting()
    timestep = result.elapsed_time/result.options.n_total
    t_show = min(result.elapsed_time,t_show)
    n_show = round(Int,t_show/timestep)
    n_mark = min(n_show,10000)
    subset = round.(Int,linspace(round(Int,n_show/n_mark),n_show,n_mark))
    PyPlot.plot(subset*timestep,result.t[subset]+rand(n_mark)*0.5-0.25,"k.",markersize=1)
    xlim(0,t_show)
    yticks(1:maximum(result.t[subset]+1))
    ylim(0,ylim()[2])
    labels("Traceplot of t","time (seconds)","t")
    draw_now()
end

function plot_t_running(result)
    checkplotting()
    o = result.options
    cdfs = t_running(result)
    tcolors = "brgycmk"^1000
    n_subset = min(2000,o.n_total)
    subset = round.(Int,linspace(1,o.n_total,n_subset))
    for i = 1:size(cdfs,2)
        PyPlot.plot(subset,cdfs[subset,i],"$(tcolors[i])-")
    end
    yticks(0:0.2:1)
    labels("Running average of CDF(t|x)", "sweep number", "CDF(t|x)")
    draw_now()
end

function plot_autocorrelation(values,t_show,timestep;kwargs...)
    checkplotting()
    a = autocorrelation(values)
    n_plot = min(round(Int,t_show/timestep+1),length(a))
    PyPlot.plot(1000*(0:n_plot-1)*timestep,a[1:n_plot],linewidth=1.5;kwargs...)
    labels("Autocorrelation","time (milliseconds)","")
    xlim(0,1000*t_show)
    ylim(0,1)
    draw_now()
end

function plot_autocorrelation(result,t_show;kwargs...)
    checkplotting()
    o = result.options
    timestep = result.elapsed_time/o.n_total
    plot_autocorrelation(result.t[o.n_burn+1:o.n_total],t_show,timestep;kwargs...)
    labels("Autocorrelation of t","time (milliseconds)","")
    draw_now()
end

function plot_t_posterior(result; kwargs...)
    checkplotting()
    pt = t_posterior(result)
    PyPlot.plot(1:length(pt), pt; kwargs...)
    labels("Posterior on t","t (# of clusters)","p(t|data)")
    xticks(1:length(pt))
    draw_now()
end

function plot_t_posterior_average(results; kwargs...)
    checkplotting()
    t_posteriors = Array{Array{Float64,1},1}()
    for result in results
        push!(t_posteriors, t_posterior(result))
    end
    T = minimum(map(length,t_posteriors))
    avg = mean([p[1:T] for p in t_posteriors])[:]
    PyPlot.plot(1:T,avg; kwargs...)
    xticks(1:T)
    labels("Average posterior on t","t (# of clusters)","p(t|data)")
    draw_now()
end

function plot_k_posterior(result; kwargs...)
    checkplotting()
    pk = k_posterior(result)
    PyPlot.plot(1:length(pk), pk; kwargs...)
    labels("Posterior on k","k (# of components)","p(k|data)")
    xticks(1:length(pk))
    draw_now()
end

function plot_k_posterior_average(results; kwargs...)
    checkplotting()
    k_posteriors = Array{Array{Float64,1},1}()
    for result in results
        push!(k_posteriors, k_posterior(result))
    end
    T = minimum(map(length,k_posteriors))
    avg = mean([p[1:T] for p in k_posteriors])[:]
    PyPlot.plot(1:T,avg; kwargs...)
    xticks(1:T)
    labels("Average posterior on k","k (# of components)","p(k|data)")
    draw_now()
end

# Posterior similarity matrix (probability that i and j are in same cluster)
function plot_similarity_matrix(result; step=10)
    checkplotting()
    n = result.options.n
    title("Posterior similarity matrix")
    imshow(1-similarity_matrix(result),cmap="gray",interpolation="nearest",vmin=0,vmax=1,extent=(1,n,n,1),origin="upper")
    tick_params(axis="x",direction="out",top="off")
    tick_params(axis="y",direction="out",top="off")
    xticks(step:step:n)
    yticks(step:step:n)
    draw_now()
end

function plot_clusters(x,z; colors=color_list, markers=marker_list, markersize=2)
    checkplotting()
    n = length(x)
    d = length(x[1]) # dimension of data
    N = zeros(n); for i = 1:n; N[z[i]] += 1; end
    zcs = sortperm(N,rev=true)
    t = countnz(N)
    if d==1 # univariate data
        for c = 1:t
            xc = x[z.==zcs[c]]
            jitter = (rand(length(xc))-0.5)/4
            PyPlot.plot(xc,c+jitter,"bx",markersize=2.0)
        end
        yticks(1:t)
        ylim(0.5,t+0.5)
        grid(true,axis="y")
    elseif d==2 # bivariate data
        T = min(length(colors),length(markers))
        if t>T; warn("Only plotting the largest $T clusters for $label."); end
        for c = 1:min(t,T)
            xc = x[z.==zcs[c]]
            PyPlot.plot([xi[1] for xi in xc],[xi[2] for xi in xc],
                "$(colors[c])$(markers[c])",markersize=markersize,markeredgewidth=0.1)
        end
    else
        error("Plotting clusters is only enabled for univariate or bivariate data.")
    end
    title("Clusters")
    draw_now()
end

function rug_plot(data) # Rug plot of data
    checkplotting()
    @assert(all(map(length,data).==1),"Rug plot only enabled for univariate data.")
    PyPlot.plot(data,zeros(length(data)),"k+",markersize = 18)
    tick_params(axis="x",direction="out",top="off")
    title("Rug plot")
    draw_now()
end 

function plot_histogram(data; kwargs...) # Histogram of data
    checkplotting()
    @assert(all(map(length,data).==1),"Histogram only enabled for univariate data.")
    counts,edges = histogram(data; n_bins=50)
    PyPlot.bar(edges[1:end-1],counts./(length(data)*diff(edges)),diff(edges); kwargs...)
    title("Histogram")
    draw_now()
end

function plot_density_estimate(result; resolution=-1, kwargs...)
    checkplotting()
    bounds(y) = (mn=minimum(y); mx=maximum(y); s=(mx-mn)/8; (round(mn-s,0),round(mx+s,0)))
    o = result.options
    data,n = o.x,o.n
    d = length(data[1]) # dimension of data
    if d==1
        # Density estimate
        if resolution==-1; resolution = 500; end
        xmin,xmax = bounds(data)
        xs = linspace(xmin,xmax,resolution)
        fs = Float64[density_estimate(xi,result) for xi in xs]
        PyPlot.plot(xs,fs; kwargs...)
    elseif d==2
        if resolution==-1; resolution = 40; end
        n_contours = 10
        power = 0.25
        data_x = Float64[data[i][1] for i=1:n]
        data_y = Float64[data[i][2] for i=1:n]
        xmin,xmax = bounds(data_x)
        ymin,ymax = bounds(data_y)
        xs,ys,X,Y = density_grid(xmin,xmax,ymin,ymax,resolution,resolution)
        subplots_adjust(top=0.845)
        Z = Float64[density_estimate(Float64[x;y],result) for x=xs, y=ys]
        PyPlot.contour(X,Y,Z.^power,n_contours)
    else
        error("Plotting density estimates is only enabled for dimensions 1 and 2.")
    end
    title("Density estimate")
    draw_now()
end


end # module BayesianMixtures


