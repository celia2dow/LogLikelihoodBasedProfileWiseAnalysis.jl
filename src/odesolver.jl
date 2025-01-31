function odesolver(times,pars,ICs,ode_system,num_X,restricted_msrmnts_per_tp)
    # Function for solving a system of ODEs

    # Dowling (2024)
    # Adapted from Murphy et al. (2023) 

    if isnothing(ICs)                               # If the ICs are unknown ...
        ode_ICs = pars[1:num_X]                     # ... they will be the first values in the array
        ode_pars = pars[num_X+1:end]                # Define the ODE parameters as the remaining values
    else
        ode_ICs = copy(ICs)
        ode_pars = copy(pars)
    end
    if !restricted_msrmnts_per_tp                   # If times is stored as a dictionary
        sol_reshaped = Dict{Int, Vector{Float64}}()
        num_X = length(times)
        for i in 1:num_X
            times_i = times[i]                      # The array of time-points for state variable i
            time_span_i=(0.0,maximum(times_i));
            prob_i=ODEProblem(ode_system,ode_ICs,time_span_i,ode_pars);
            sol_i=solve(prob_i,saveat=times_i);                 # generate solution from known parameters
            sol_reshaped_i= reshape(sol_i,num_X,length(times_i))  # reshape array into vector
            sol_reshaped[i] = sol_reshaped_i[i,:]
        end
    else                                            # If times is a vector
        time_span=(0.0,maximum(times));
        prob=ODEProblem(ode_system,ode_ICs,time_span,ode_pars);
        sol=solve(prob,saveat=times);                  # generate solution from known parameters
        sol_reshaped = reshape(sol,num_X,length(times)) # reshape array into vector
    end
    return sol_reshaped #sol_reshaped[1:num_X,:];
end
    