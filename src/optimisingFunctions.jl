#######################################################################################
## Functions involved in maximising the log-likelihood function
#######################################################################################

function sum_loglikelihood(time,params_combined,ICs,ode_system,data,error_type)
    # Sum over the loglikelihood values at each time point for a given set of parameters
    if error_type == "Normal" || error_type == "Lognormal"
        pars = params_combined[1:end-1]     # Extract the model parameters
    end
    y=odesolver(time,pars,ICs,ode_system);  # Find the model solution for the given parameter values
    e=0;                                    # Initialise sum of loglikelihood values

    if error_type == "Normal"
        Sd = params_combined[end]
        model_params = params_combined[1:end-1]
        data_dists=[Normal(mi,Sd) for mi in y];     # Generate Normal distributions for each time point
        e+=sum([loglikelihood(data_dists[i],data[i]) for i in 1:length(data_dists)]) 
    elseif error_type == "Lognormal"
        Sd = params_combined[end]
        model_params = params_combined[1:end-1]
        data_dists=[LogNormal(0,a[3]) for mi in y]; # Generate LogNormal distributions for each time point
        e+=sum([loglikelihood(data_dists[i],data[i]./y[i]) for i in 1:length(data_dists)]) 
    end

    return e
end


function optimise(fun,θ₀,lb,ub) 
    # Find the paramter values θ that maximise the loglikelihood of the data 
    tomax = (θ,∂θ) -> fun(θ)                     # Function to optimise
    opt = Opt(:LN_NELDERMEAD,length(θ₀))
    opt.max_objective = tomax                   # Optimisation method is maximisation
    opt.lower_bounds = lb                       # Lower bound
    opt.upper_bounds = ub                       # Upper bound
    opt.maxtime = 30.0;                         # maximum time in seconds
    res = optimize(opt,θ₀)                      # Results of optimisation include: 
                                                # the otpimising paramter values θ and 
                                                # the value of the loglikelihood sum given θ
    return res[[2,1]]
end

