function odesolver(times,pars,ICs,ode_system)
    # Function for solving a system of ODEs

    # Dowling (2024)
    # Adapted from Murphy et al. (2023) 

    num_x = length(ICs);                            # Find the number of state variables
    time_span=(0.0,maximum(times));
    prob=ODEProblem(ode_system,ICs,time_span,pars);
    sol=solve(prob,saveat=times,);                  # generate solution from known parameters
    sol_reshaped = reshape(sol,num_x,length(times)) # reshape array into vector
    return sol_reshaped[1:num_x,:];
end
    