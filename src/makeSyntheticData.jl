function makeSyntheticData(args...)
    # Generates time series synthetic data using user-defined models and parameters.
    
    # Dowling (2024)
    # Adapted from Murphy et al. (2023) 

    # Assign arguments to appropriate labels
    num_args = length(args) # Check the number of arguemtns
    if num_args == 5        # The case with no noise
        time, ICs, ode_system, model_params, noNegatives = args
        error_type = "NoNoise"
    elseif num_args == 6    # The case with no noise paramter inputs
        time, ICs, ode_system, model_params, error_type, noNegatives = args
    elseif num_args == 7    # The case with noise
        time, ICs, ode_system, model_params, error_type, noise_params, noNegatives = args
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
        sol1 = odesolver(time,model_params,ICs,ode_system)
        data_dists=[Normal(mu,Sd) for mu in sol1];   # Generate normal distributions for each time point
        synthetic_data=[rand(data_dist) for data_dist in data_dists];                       # Generate the noisy synthetic data from these distributions
        if noNegatives
            # Set any unrealistically negative values to zero
            synthetic_data[synthetic_data .< 0] .= zero(eltype(synthetic_data))                 # Set any negative values to zero
        end
    elseif error_type == "Lognormal"
        Sd = noise_params
        synthetic_data = [mu*rand(LogNormal(0,Sd)) for mu in odesolver(time,model_params,ICs,ode_system)]; # Generate the noisy synthetic data
    elseif error_type == "NoNoise"
        synthetic_data = odesolver(time,model_params,ICs,ode_system)  # Generate smooth synthetic data
    end

    return synthetic_data
end