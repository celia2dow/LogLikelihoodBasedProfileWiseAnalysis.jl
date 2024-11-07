function odesolver(time,pars,ICs,ode_system)
    # Function for solving a system of ODEs

    # Dowling (2024)
    # Adapted from Murphy et al. (2023) 

    num_x = length(ICs);               # Find the number of state variables
    time_span=(0.0,maximum(time));
    prob=ODEProblem(ode_system,ICs,time_span,pars);
    sol=0.0*zeros(num_x,length(time)); # generate vector to solution 
    sol=solve(prob,saveat=time,);       # generate solution from known parameters

    return sol[1:num_x,:];
end
    