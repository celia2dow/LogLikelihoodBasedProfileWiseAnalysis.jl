function makeSyntheticData(args...)
    # Generates time series synthetic data using user-defined models and parameters.
    
    # Dowling (2024)
    # Adapted from Murphy et al. (2023) 

    # Assign arguments to appropriate labels
    num_args = length(args) # Check the number of arguemtns
    if num_args == 4        # The case with no noise
        time, ICs, ode_system, model_params = args
        error_type = "NoNoise"
    elseif num_args == 5    # The case with no noise inputted as a parameter
        time, ICs, ode_system, model_params, error_type = args
    elseif num_args == 6    # The case with noise
        time, ICs, ode_system, model_params, error_type, noise_params = args
    else
        println("makeSyntheticData: not enough inputs \n")
    end
    num_x = length(ICs);    # Find the number of state variables
    
    # Preallocate space
    if num_x == 1
        synthetic_data=0.0*zeros(length(time));          # Column vector
    else
        synthetic_data=0.0*zeros(num_x,length(time));    # One state variable per row
    end
    
    # Generate synthetic data
    if error_type == "Normal"
        Sd = noise_params
        data_dists=[Normal(mi,Sd) for mi in odesolver(time,model_params,ICs,ode_system)];   # Generate normal distributions for each time point
        synthetic_data=[rand(data_dist) for data_dist in data_dists];                       # Generate the noisy synthetic data from these distributions
        synthetic_data[synthetic_data .< 0] .= zero(eltype(synthetic_data))                 # Set any negative values to zero
    elseif error_type == "Lognormal"
        Sd = noise_params
        synthetic_data = [mi*rand(LogNormal(0,Sd)) for mi in odesolver(time,model_params,ICs,ode_system)]; # Generate the noisy synthetic data
    elseif error_type == "NoNoise"
        synthetic_data = odesolver(time,model_params,ICs,ode_system)  # Generate smooth synthetic data
    end

    if num_x == 1
        synthetic_data = synthetic_data' # Make column vector a row vector
    end
   
    return synthetic_data
end