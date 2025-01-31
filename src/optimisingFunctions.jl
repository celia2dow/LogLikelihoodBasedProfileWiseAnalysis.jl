#######################################################################################
## Functions involved in maximising the log-likelihood function
#######################################################################################

function sum_loglikelihood(times,params_combined,ICs,ode_system,data,error_type,num_X,restricted_msrmnts_per_tp)
    # Sum over the loglikelihood values at each time point for a given set of parameters
    if error_type == "Normal" || error_type == "Lognormal"
        pars = params_combined[1:end-1]             # Extract the model parameters
    end
    y=odesolver(times,pars,ICs,ode_system,num_X,restricted_msrmnts_per_tp);   # Find the model solution for the given parameter values
    e=0;                                            # Initialise sum of loglikelihood values
    # Loop through all state variables
    for j in 1:num_X
        if !restricted_msrmnts_per_tp               # If measurements are not restricted per time-point and data is stored in a dictionary
            data_1D = data[j]                       # Exctract into a 1D array
            y_1D = y[j]
        else                                        # Otherwise if measurements are restricted per time-point
            data_1D = data[j,:]                     # Extract the row
            y_1D = y[j,:]
        end
        # Calculate the log-likeilhood value for each data-point and add it to the total
        if error_type == "Normal"
            Sd = params_combined[end]
            data_dists=[Normal(mi,Sd) for mi in y_1D];     # Generate Normal distributions for each times point
            e+=sum([loglikelihood(data_dists[i],data_1D[i]) for i in eachindex(data_dists)]) 
        elseif error_type == "Lognormal"
            Sd = params_combined[end]
            data_dists=[LogNormal(0,Sd) for mi in y_1D]; # Generate LogNormal distributions for each times point
            e+=sum([loglikelihood(data_dists[i],data_1D[i]./y_1D[i]) for i in eachindex(data_dists)]) # RYAN TO IMPROVE
        end
    end

    return e
end


function optimise(fun,θ₀,lb,ub) 
    # Find the paramter values θ that maximise the loglikelihood of the data 
    tomax = (θ,∂θ) -> fun(θ)                    # Function to optimise
    opt = Opt(:LN_NELDERMEAD,length(θ₀))
    opt.max_objective = tomax                   # Optimisation method is maximisation
    opt.lower_bounds = lb                       # Lower bound
    opt.upper_bounds = ub                       # Upper bound
    opt.maxtime = 180.0;                        # maximum times in seconds 
    res = optimize(opt,θ₀)                      # Results of optimisation include: 
                                                # the otpimising paramter values θ and 
                                                # the value of the loglikelihood sum given θ
    # PRINT OUTPUT TO SEE IF OPTIMA IS REACHED
    # @time res = optimize(opt,θ₀)
    return res[[2,1]]
end

