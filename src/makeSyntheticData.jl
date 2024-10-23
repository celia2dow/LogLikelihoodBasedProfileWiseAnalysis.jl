function makeSyntheticData(args...)
    # Generates time series data using user-defined models and parameters.
    
    # Dowling (2024)
    # Adapted from Murphy et al. (2023) 

    # Assign arguments to appropriate labels
    num_args = length(args) # Check the number of arguemtns
    if num_args == 4        # The case with no noise
        time, ICs, ode_system, model_params = args
        error_type = "NoNoise"
    elseif num_args == 6    # The case with noise
        time, ICs, ode_system, model_params, error_type, noise_params = args
    end
    num_x = length(ICs);    # Find the number of state variables
    
    # Preallocate space
    if num_x == 1
        data=0.0*zeros(length(t));          # Column vector
    else
        data=0.0*zeros(num_x,length(t));    # One state variable per row
    end
    
    # Generate data
    if error_type == "Normal"
        Sd = noise_params
        data_dists=[Normal(mi,Sd) for mi in odesolver(time,model_params,ICs,ode_system)];   # Generate normal distributions for each time point
        data=[rand(data_dist) for data_dist in data_dists];                                 # Generate the noisy data from these distributions
        data[data .< 0] .= zero(eltype(data))                                               # Set any negative values to zero
    elseif error_type == "Lognormal"
        Sd = noise_params
        data = [mi*rand(LogNormal(0,Sd)) for mi in odesolver(time,model_params,ICs,ode_system)]; # Generate the noisy data
    elseif error_type == "NoNoise"
        data = odesolver(time,model_params,ICs,ode_system)  # Generate smooth data
    end

    if num_x == 1
        data = data' # Make column vector a row vector
    end
   
end