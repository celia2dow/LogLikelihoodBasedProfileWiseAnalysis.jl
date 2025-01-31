function pwaFunction(experi_name, X_variables, X_init_condits, model_params, error_type, noise_params, 
    ode_system, times_all, times_specific, data_all, data_specific, npts, seed_num)
    # Implementing an input measurement error models with an input mechanistic model of ODEs or PDEs
    # in a likelihood-based framework for estimation, identifiability analysis, 
    # and prediction. Note that the code is ideal for visualising a model with up to 6 state variables.
    # If more state variables are desired, add more colours to the colours array. The code will still
    # work without this update, however.

    # Dowling (2024)
    # Adapted from Murphy et al. (2023) 
    
    #######################################################################################
    ### Initialisation including plot options, and filepaths to save outputs
    fnt = Plots.font("sans-serif", 12)    
    mrkrsize = 2 
    mrkroutlinewidth = 0.1                                                                        # plot options
    global cur_colors = palette(:default)                                                       # plot options
    isdir(pwd() * "/ProfileWiseAnalysisOutput/"* experi_name *"/") || mkdir(pwd() * "/ProfileWiseAnalysisOutput/"* experi_name) # make folder to save figures if doesnt already exist
    filepath_save = [pwd() * "/ProfileWiseAnalysisOutput/"* experi_name *"/"]                   # location to save figures "//ProfileWiseAnalysisOutput//experi_name//"]
    mm_to_pts_scaling = 283.46/72;
    fig_2_across_size=(75.0*mm_to_pts_scaling,50.0*mm_to_pts_scaling);
    # fig_3_across_size=(50.0*mm_to_pts_scaling,33.3*mm_to_pts_scaling);
    # fig_4_across_size=(37.5*mm_to_pts_scaling,25.0*mm_to_pts_scaling);

    #######################################################################################
    ### Set random seed
    Random.seed!(seed_num)

    #######################################################################################
    ### Define values, names, initial conditions, and guesses
    Xi_string = join(X_variables.name, ", ")    # String of state variable(s)
    num_X = length(X_variables.name);           # Number of state variable(s)

    guess_r = model_params.guess                # Initial guess(es) for model parameter(s)
    guess_er = noise_params.guess               # Initial guess(es) for noise parameter(s)
    
    lb_r = model_params.lower_bound             # Lower bound(s) for model parameter(s)
    lb_er = noise_params.lower_bound            # Lower bound(s) for noise parameter(s)
         
    ub_r = model_params.upper_bound             # Upper bound(s) for model parameter(s)
    ub_er = noise_params.upper_bound            # Upper bound(s) for noise parameter(s)
    
    ICs = X_init_condits.true_value             # Initial condition(s) if known
    if isnothing(X_init_condits.guess)
        θG = [guess_r; guess_er]                # Combined initial guesses
        lb = [lb_r; lb_er];                     # Combined lower bounds
        ub = [ub_r; ub_er];                     # Combined upper bounds 
        lb = float.(lb); ub = float.(ub);       # Enforce bounds to be type float
        all_param_names = [                     # All parameter names
            model_params.name; noise_params.name]
        num_to_est =                            # Number of model and noise parameters to estimate
            length(guess_r) + length(guess_er)  
    else 
        guess_IC = X_init_condits.guess         # Initial guess(es) for initial condition(s)
        lb_IC = X_init_condits.lower_bound      # Lower bound(s) for initial condition(s)
        ub_IC = X_init_condits.upper_bound      # Upper bound(s) for initial condition(s)
        
        θG = [guess_IC; guess_r; guess_er]      # Combined initial guesses
        lb = [lb_IC; lb_r; lb_er];              # Combined lower bounds
        ub = [ub_IC; ub_r; ub_er];              # Combined upper bounds 
        lb = float.(lb); ub = float.(ub);       # Enforce bounds to be type float
        all_param_names = [                     # All parameter names
            X_init_condits.name; model_params.name; noise_params.name]
        num_to_est =                            # Number of model and noise parameters and ICs to estimate
            length(guess_r) + length(guess_er) + length(guess_IC)
    end
    
    num_known_parameters = 0;                   # Track the number of known parameters provided to the function
    
    if !isnothing(model_params.true_value)      # If the model parameters are known ...
        r_known = model_params.true_value       # ... count the number of model parameter values and add them to the total number of known paramters
        num_known_parameters = num_known_parameters + length(r_known)
    else                                        # Otherwise ...
        r_known = nothing                       # ... set the number of model parameter values to 'nothing' for future case-checking
    end

    if !isnothing(noise_params.true_value)      # If the error model parameters are known ...
        er_known = noise_params.true_value      # ... count the number of error model parameter values and add them to the total number of known parameters
        num_known_parameters = num_known_parameters +length(er_known)
    else                                        # Otherwise ...
        er_known = nothing                      # ... set the number of error model parameter values to 'nothing' for future case-checking
    end

    if !isnothing(X_init_condits.true_value)    # If the initial conditions are known ...
        ic_known = X_init_condits.true_value    # ... count the number of IC values and add them to the total number of known parameters
        num_known_parameters = num_known_parameters +length(ic_known)
    else                                        # Otherwise ...
        ic_known = nothing                      # ... set the number of IC values to 'nothing' for future case-checking
    end

    restricted_msrmnts_per_tp = true               # Check whether the measurements have been restricted for each timepoint: 1 = they have been restricted
    if typeof(data_specific) == Dict{Int64, Vector{Float64}}
        restricted_msrmnts_per_tp = false           # 0 = they have not been restricted, so times_specific and data_specific are dictionaries
    end

    #######################################################################################
    ### Generate time-points for plotting known values and best fits
    if typeof(times_all) == Dict{Int64, Vector{Float64}}    # If the complete data-set for tumour growth is a dictionary
        t_max = 0                                           # Initialise
        times_all_vec = []
        data_max = -Inf                                     
        for i in 1:num_X 
            times_all_i = times_all[i]                                 
            t_max = max(t_max, times_all_i[end])            # Find the largest time for each state variable times array
            times_all_vec = [times_all_vec;times_all_i[:]]  # Append all time-points to a vector
            data_max = max(data_max,maximum(data_all[i]))   # Find the largest data value for each state variable data array
        end
        unique_times_all = unique(times_all_vec)'           # Create a list of all of the unique time-points
    
    else                                                    # If the complete data-set for tumour growth is a regular matrix
        t_max = times_all[end]                              # Find the largest time
        unique_times_all = unique(times_all)'               # Create a list of all of the unique time-points
        data_max = maximum(data_all)
    end
    
    time_smooth =  LinRange(0.0,t_max*1.1, 201);

    #######################################################################################
    ### Confidence thresholds
    TH=-1.921; #95% confidence interval threshold (for confidence interval for model parameters, and confidence set for model solutions)
    TH_realisations = -2.51 # 97.5 confidence interval threshold (for model realisations)
    THalpha = 0.025; # 97.5% data realisations interval threshold

    #######################################################################################
    ### Plot presentation preferences
    # Plot limits - add a bit of leeway in the axes
    xmin =-0.01*t_max                        
    xmax = t_max + 0.1 * abs(t_max)
    ymin_dat = - 0.1 * abs(data_max)            # For data
    ymax_dat = data_max + 0.1 * abs(data_max)
    ymin_prof = -4;                             # For likelihood profiles
    ymax_prof = 0.1;
    yticks_prof = [-3,-2,-1,0]
    # Define colours - FEEL FREE TO ADD CUSTOM COLOURS HERE
    colour1 = :green2                           # Bright green
    colour2 = :magenta                          # Magenta
    colour3 = RGB(0.8, 0.7, 0.2);               # Mustard yellow
    colour4 = RGB(0.0, 1.0, 1.0);               # Cyan
    colour5 = RGB(0.90196, 0.90196, 0.98039);   # Lavender
    colour6 = RGB(1.0, 0.8, 0.6);               # Peach
    colour = [colour1, colour2, colour3, colour4, colour5, colour6]
    # Define more colours if necessary
    if num_X>length(colour)
        for i in 1:num_X-length(colour)
            new_colour = RGB(rand(), rand(), rand())   # Create a random new colour
            global colour = [colour, new_colour]             # Add it to the array of colours
        end
    end
    # Define slightly darker colours for the outlines of markers
    mrkroutlinecolour = [get_darker_colour(colour_i) for colour_i in colour]

    #######################################################################################
    ### Generate smooth data/true solution if true parameter parameter values are known
    if !isnothing(r_known) && !isnothing(ICs)
        data0_smooth = odesolver(time_smooth,r_known,ICs,ode_system,num_X,true)
        if size(data0_smooth)[1] > num_X
            data0_smooth = data0_smooth'    # Convert to series of row vectors if not already
        end
    end

    #######################################################################################
    ### Plot data - scatter plot

    # Plot data against true solution (if available)
    p0scatter = plot(                                                                       # Create figure for data (against true solution)
        xlab=L"t", ylab=L"%$Xi_string",
        legend=false,xlims=(xmin,xmax),ylims=(ymin_dat,ymax_dat),size=fig_2_across_size,
        titlefont=fnt,guidefont=fnt,tickfont=fnt,legendfont=fnt)   

    if !isnothing(r_known) && !isnothing(ICs)                                               # If the true parameter values are known
        for i in 1:num_X                                                                    # Plot true trajectories
            p0scatter=plot!(time_smooth,data0_smooth[i,:],lw=3,linecolor=colour[i])         # solid - Xi(t)
        end
    end

    for i in 1:num_X                                                                        # Plot data
        if restricted_msrmnts_per_tp    # If measurements have been restricted to the same number per time-point
            p0scatter=scatter!(times_specific[:],data_specific[i,:],markersize = mrkrsize,markercolor=colour[i], markerstrokewidth=mrkroutlinewidth, markerstrokecolor=mrkroutlinecolour[i]) # scatter - Xi(t)
        else                            # Otherwise if measurements have not been restricted to the same number per time-point
            times_specific_i = times_specific[i]
            data_specific_i = data_specific[i]
            #p0scatter=scatter!(times_specific_i[:],data_specific_i[:],markersize = mrkrsize,markercolor=colour[i], markerstrokewidth=mrkroutlinewidth, markerstrokecolor=mrkroutlinecolour[i]) # scatter - Xi(t)
            p0scatter=scatter!(times_specific[i],data_specific[i],markersize = mrkrsize,markercolor=colour[i], markerstrokewidth=mrkroutlinewidth, markerstrokecolor=mrkroutlinecolour[i]) # scatter - Xi(t)
        end
    end

    display(p0scatter)
    savefig(p0scatter,filepath_save[1] * "Fig0scatter" * ".pdf")

    #######################################################################################
    ### Parameter identifiability analysis - MLE, profiles
        
        # MLE

        # MLE optimisation 
        function funMLE(unknowns)
            # Define function whose inputs are the unknown parameters, and whose output is 
            # that of sum_loglikelihood
            if isnothing(ICs)                           # If the ICs are unknown ... 
                ode_ICs = unknowns[1:num_X]             # ... they will be the first values in the array
                unknown_pars = unknowns[num_X+1:end]    # Define the ODE parameters as the remaining values
            else
                ode_ICs = copy(ICs)
                unknown_pars = copy(unknowns)
            end
            return sum_loglikelihood(times_specific,unknown_pars,ode_ICs,ode_system,data_specific,error_type,num_X,restricted_msrmnts_per_tp)
        end
        (xopt,fopt)  = optimise(funMLE,θG,lb,ub) 

        # Storing MLE 
        rmle = Array{Float64}(undef,num_to_est);        # Initialise array for storing the running MLE estimates of parameters
        global fmle=fopt
        for i in 1:num_to_est
            global rmle[i] = xopt[i]
        end
        data0_smooth_MLE = odesolver(time_smooth,rmle,ICs,ode_system,num_X,true); # generate data using MLE parameter values for plotting

        # Plot model simulated at MLE and data
        p1 = plot(                                                                                  # Create figure for MLE best fit
            xlab=L"t", ylab=L"%$Xi_string", 
            legend=false,xlims=(xmin,xmax),ylims=(ymin_dat,ymax_dat),size=fig_2_across_size,
            titlefont=fnt,guidefont=fnt,tickfont=fnt,legendfont=fnt)

        for i in 1:num_X                                                                            # Plot true trajectories
            p1=plot!(time_smooth,data0_smooth_MLE[i,:],lw=3,linecolor=colour[i])                    # solid - Xi(t)
        end

        for i in 1:num_X                                                                            # Plot data
            if restricted_msrmnts_per_tp    # If measurements have been restricted to the same number per time-point
                p1=scatter!(times_specific[:],data_specific[i,:],markersize = mrkrsize,markercolor=colour[i], markerstrokewidth=mrkroutlinewidth, markerstrokecolor=mrkroutlinecolour[i]) # scatter - Xi(t)
            else                            # Otherwise if measurements have not been restricted to the same number per time-point
                times_specific_i = times_specific[i]
                data_specific_i = data_specific[i]
                #p1=scatter!(times_specific_i[:],data_specific_i[:],markersize = mrkrsize,markercolor=colour[i], markerstrokewidth=mrkroutlinewidth, markerstrokecolor=mrkroutlinecolour[i]) # scatter - Xi(t)
                p1=scatter!(times_specific[i],data_specific[i],markersize = mrkrsize,markercolor=colour[i], markerstrokewidth=mrkroutlinewidth, markerstrokecolor=mrkroutlinecolour[i]) # scatter - Xi(t)
            end
        end

        display(p1)
        savefig(p1,filepath_save[1] * "Fig1" * ".pdf")

    #######################################################################################
        ### Profiling (using the MLE as the first guess at each point)

        # Preallocate space for storing arrays with dimensions:
        #   num_to_est - the total number of model and error model parameters
        #   num_X - the total number of state variables
        #   length(time_smooth) - the number of timesteps at which the solution has been solved for
        max_ri_store = zeros(num_to_est,num_X,length(time_smooth))
        min_ri_store = zeros(num_to_est,num_X,length(time_smooth))
        max_realisations_ri_store = zeros(num_to_est,num_X,length(time_smooth))
        min_realisations_ri_store = zeros(num_to_est,num_X,length(time_smooth))
        ri_range_store = zeros(num_to_est,npts*2)
        nllri_store = zeros(num_to_est,npts*2)

        for i in 1:num_to_est
            # ri profile
            ri_min = lb[i]
            ri_max = ub[i]
            ri_range_lower = reverse(LinRange(ri_min,rmle[i],npts)) # OR ri_range_lower = ri_min:(rmle[i]-ri_min)/(npts-1):rmle[i]
            ri_range_upper = LinRange(rmle[i] + (ri_max-rmle[i])/npts,ri_max,npts) 

            non_ri_range_lower = zeros(num_to_est-1,npts)
            llri_lower = zeros(npts)                                # For storing the log-likelihood value associated with that value of ri (lower half of ri values)
            predict_ri_lower=zeros(num_X,length(time_smooth),npts)
            predict_ri_realisations_lower_lq=zeros(num_X,length(time_smooth),npts)
            predict_ri_realisations_lower_uq=zeros(num_X,length(time_smooth),npts)

            non_ri_range_upper = zeros(num_to_est-1,npts)
            llri_upper = zeros(npts)                                # For storing the log-likelihood value associated with that value of ri (upper half of ri values)
            predict_ri_upper=zeros(num_X,length(time_smooth),npts)
            predict_ri_realisations_upper_lq=zeros(num_X,length(time_smooth),npts)
            predict_ri_realisations_upper_uq=zeros(num_X,length(time_smooth),npts)

            # Lower and upper bounds for parameters that are not ri
            lb_i = copy(lb); deleteat!(lb_i,i);
            ub_i = copy(ub); deleteat!(ub_i,i);

            # Start at MLE and increase parameter (upper)
            for j in 1:npts
                function fun_upper(aa)
                    # Create an array of the unknown parameters
                    aa_copy = copy(aa)
                    # Add the single known parameter in the array in the appropriate order
                    insert!(aa_copy, i, ri_range_upper[j])
                    # Define the ICs separately to the other parameters
                    if isnothing(ICs)                           # If the ICs are unknown ... 
                        ode_ICs = aa_copy[1:num_X]             # ... they will be the first values in the array
                        unknown_pars = aa_copy[num_X+1:end]    # Define the ODE parameters as the remaining values
                    else
                        ode_ICs = copy(ICs)
                        unknown_pars = copy(aa_copy)
                    end
                    # Evaluate the sum of loglikelihoods for this combination of parameters
                    return sum_loglikelihood(times_specific,unknown_pars,ode_ICs,ode_system,data_specific,error_type,num_X,restricted_msrmnts_per_tp)
                end

                # Maximise the loglikelihood for each value of ri in the upper range to find the corresponding values of the non-ri parameters
                local θG_i=copy(rmle); deleteat!(θG_i,i) # Use MLE values as the guess for non-ri parameters
                
                local (xo,fo)=optimise(fun_upper,θG_i,lb_i,ub_i)
                non_ri_range_upper[:,j]=xo[:]
                llri_upper[j]=fo[1]
                
                # Find model solutions given the set value of ri, and the optimised non-ri values
                if error_type == "Normal" || error_type == "Lognormal"
                    if i != num_to_est
                        model_params_only = copy(non_ri_range_upper[1:end-1,j])
                        model_params_only = insert!(model_params_only,i,ri_range_upper[j])
                        Sd_j = non_ri_range_upper[end,j]
                    elseif i == num_to_est
                        model_params_only = non_ri_range_upper[:,j]
                        Sd_j = ri_range_upper[j]
                    end
                end
                modelmean = odesolver(time_smooth,model_params_only,ICs,ode_system,num_X,true)
                predict_ri_upper[:,:,j]=modelmean; 

                # Find α/2 and 1-α/2 quantiles of model realisations given the set value of ri, and the optimised non-ri values
                if error_type == "Normal" 
                    loop_data_dists=[Normal(mi,Sd_j) for mi in modelmean]; # Normal distribution about mean
                    predict_ri_realisations_upper_lq[:,:,j]= [quantile(data_dist, THalpha/2) for data_dist in loop_data_dists];
                    predict_ri_realisations_upper_uq[:,:,j]= [quantile(data_dist, 1-THalpha/2) for data_dist in loop_data_dists];
                elseif error_type == "Lognormal" 
                    loop_data_dists=[LogNormal(0,Sd_j) for mi in modelmean]; # LogNormal distribution about mean
                    predict_ri_realisations_upper_lq[:,:,j]= [quantile(data_dist, THalpha/2) for data_dist in loop_data_dists].*modelmean;
                    predict_ri_realisations_upper_uq[:,:,j]= [quantile(data_dist, 1-THalpha/2) for data_dist in loop_data_dists].*modelmean;
                end

                # If a new maximum sum of loglikelihoods has been discovered, update the MLE estimates
                if fo > fmle
                    println("new MLE r$i upper")
                    global fmle=fo
                    global rmle[i] = ri_range_upper[j]
                    indices = [1:i-1; i+1:num_to_est]
                    for l1 in indices
                        if l1>i
                            l2=l1-1   # Reduce index if parameter comes after the current profile's parameter
                        else
                            l2=l1
                        end
                        global rmle[l1] = non_ri_range_upper[l2,j]
                    end
                end
            end

            # Start at MLE and decrease parameter (lower)
            for j in 1:npts
                function fun_lower(aa)
                    # Create an array of the unknown parameters
                    aa_copy = copy(aa)
                    # Add the single known parameter in the array in the appropriate order
                    insert!(aa_copy, i, ri_range_lower[j])
                    # Define the ICs separately to the other parameters
                    if isnothing(ICs)                           # If the ICs are unknown ... 
                        ode_ICs = aa_copy[1:num_X]             # ... they will be the first values in the array
                        unknown_pars = aa_copy[num_X+1:end]    # Define the ODE parameters as the remaining values
                    else
                        ode_ICs = copy(ICs)
                        unknown_pars = copy(aa_copy)
                    end
                    # Evaluate the sum of loglikelihoods for this combination of parameters
                    return sum_loglikelihood(times_specific,unknown_pars,ode_ICs,ode_system,data_specific,error_type,num_X,restricted_msrmnts_per_tp)
                end
                
                # Find the non-ri parameter values that maximise the loglikelihood for each value of ri in the upper range
                local θG_i=copy(rmle); deleteat!(θG_i,i)
                local (xo,fo)=optimise(fun_lower,θG_i,lb_i,ub_i)
                non_ri_range_lower[:,j]=xo[:]
                llri_lower[j]=fo[1]
                
                # Find model solutions given the set value of ri, and the optimised non-ri values
                if error_type == "Normal" || error_type == "Lognormal"
                    if i != num_to_est
                        model_params_only = copy(non_ri_range_lower[1:end-1,j])
                        model_params_only = insert!(model_params_only,i,ri_range_lower[j])
                        Sd_j = non_ri_range_lower[end,j]
                    elseif i == num_to_est
                        model_params_only = non_ri_range_lower[:,j]
                        Sd_j = ri_range_lower[j]
                    end
                end
                modelmean = odesolver(time_smooth,model_params_only,ICs,ode_system,num_X,true)
                predict_ri_lower[:,:,j]=modelmean; 
    
                # Find α/2 and 1-α/2 quantiles of model realisations given the set value of ri, and the optimised non-ri values
                if error_type == "Normal" 
                    loop_data_dists=[Normal(mi,Sd_j) for mi in modelmean]; # Normal distribution about mean
                    predict_ri_realisations_lower_lq[:,:,j]= [quantile(data_dist, THalpha/2) for data_dist in loop_data_dists];
                    predict_ri_realisations_lower_uq[:,:,j]= [quantile(data_dist, 1-THalpha/2) for data_dist in loop_data_dists];
                elseif error_type == "Lognormal" 
                    loop_data_dists=[LogNormal(0,Sd_j) for mi in modelmean]; # LogNormal distribution about mean
                    predict_ri_realisations_lower_lq[:,:,j]= [quantile(data_dist, THalpha/2) for data_dist in loop_data_dists].*modelmean;
                    predict_ri_realisations_lower_uq[:,:,j]= [quantile(data_dist, 1-THalpha/2) for data_dist in loop_data_dists].*modelmean;
                end

                # If a new maximum sum of loglikelihoods has been discovered, update the MLE estimates
                if fo > fmle
                    println("new MLE r$i lower")
                    global fmle=fo
                    global rmle[i] = ri_range_lower[j]
                    indices = [1:i-1; i+1:num_to_est]
                    for l1 in indices
                        if l1>i
                            l2=l1-1   # Reduce index if parameter comes after the current profile's parameter
                        else
                            l2=l1
                        end
                        global rmle[l1] = non_ri_range_lower[l2,j]
                    end
                end
            end

            # combine the lower and upper
            ri_range_store[i,:] = [reverse(ri_range_lower); ri_range_upper]
            llri = [reverse(llri_lower); llri_upper] 
            nllri_store[i,:] = llri.-maximum(llri);

            max_ri = -Inf * ones(num_X,length(time_smooth))
            min_ri = Inf * ones(num_X,length(time_smooth))
            for j in 1:(npts)
                if (llri_lower[j].-maximum(llri)) >= TH
                    for k in eachindex(time_smooth)
                        for l in 1:num_X
                            max_ri[l,k]=max(predict_ri_lower[l,k,j],max_ri[l,k])
                            min_ri[l,k]=min(predict_ri_lower[l,k,j],min_ri[l,k])
                        end
                    end
                end
            end
            for j in 1:(npts)
                if (llri_upper[j].-maximum(llri)) >= TH
                    for k in eachindex(time_smooth)
                        for l in 1:num_X
                            max_ri[l,k]=max(predict_ri_upper[l,k,j],max_ri[l,k])
                            min_ri[l,k]=min(predict_ri_upper[l,k,j],min_ri[l,k])
                        end
                    end
                end
            end

            max_ri_store[i,:,:] = max_ri
            min_ri_store[i,:,:] = min_ri

            # combine the confidence sets for realisations
            max_realisations_ri = -Inf * ones(num_X,length(time_smooth))
            min_realisations_ri = Inf * ones(num_X,length(time_smooth))
            for j in 1:(npts)
                if (llri_lower[j].-maximum(llri)) >= TH_realisations
                    for k in eachindex(time_smooth)
                        for l in 1:num_X
                            max_realisations_ri[l,k]=max(predict_ri_realisations_lower_uq[l,k,j],max_realisations_ri[l,k])
                            min_realisations_ri[l,k]=min(predict_ri_realisations_lower_lq[l,k,j],min_realisations_ri[l,k])
                        end
                    end
                end
            end
            for j in 1:(npts)
                if (llri_upper[j].-maximum(llri)) >= TH_realisations
                    for k in eachindex(time_smooth)
                        for l in 1:num_X
                            max_realisations_ri[l,k]=max(predict_ri_realisations_upper_uq[l,k,j],max_realisations_ri[l,k])
                            min_realisations_ri[l,k]=min(predict_ri_realisations_upper_lq[l,k,j],min_realisations_ri[l,k]) 
                        end
                    end
                end
            end
            max_realisations_ri_store[i,:,:] = max_realisations_ri
            min_realisations_ri_store[i,:,:] = min_realisations_ri
        end
        
    #######################################################################################
        ### interpolate for smoother profile likelihoods
        interp_npts= 1001;

        # Recompute best fit 
        model_params_only = rmle[eachindex([guess_IC; guess_r])]                                        # Use MLE values of model parameters for best fit
        data0_smooth_MLE_recomputed = odesolver(time_smooth,model_params_only,ICs,ode_system,num_X,true);    # generate best fit
        if size(data0_smooth_MLE_recomputed)[1] > num_X
            data0_smooth_MLE_recomputed = data0_smooth_MLE_recomputed'                                  # Convert to series of row vectors if not already
        end

        # Plot model simulated at MLE and data
        p1_updated = plot(                                                                                  # Create figure for MLE best fit
            xlab=L"t", ylab=L"%$Xi_string", 
            legend=false,xlims=(xmin,xmax),ylims=(ymin_dat,ymax_dat),size=fig_2_across_size,
            titlefont=fnt,guidefont=fnt,tickfont=fnt,legendfont=fnt) 

        for i in 1:num_X                                                                                    # Plot best fit trajectories
            p1_updated=plot!(time_smooth,data0_smooth_MLE_recomputed[i,:],lw=3,linecolor=colour[i])         # solid - Xi(t)
        end

        for i in 1:num_X                                                                                    # Plot data
            if restricted_msrmnts_per_tp    # If measurements have been restricted to the same number per time-point
                p1_updated=scatter!(times_specific[:],data_specific[i,:],markersize = mrkrsize,markercolor=colour[i], markerstrokewidth=mrkroutlinewidth, markerstrokecolor=mrkroutlinecolour[i]) # scatter - Xi(t)
            else                            # Otherwise if measurements have not been restricted to the same number per time-point
                times_specific_i = times_specific[i]
                data_specific_i = data_specific[i]
                #p1_updated=scatter!(times_specific_i[:],data_specific_i[:],markersize = mrkrsize,markercolor=colour[i], markerstrokewidth=mrkroutlinewidth, markerstrokecolor=mrkroutlinecolour[i]) # scatter - Xi(t)
                p1_updated=scatter!(times_specific[i],data_specific[i],markersize = mrkrsize,markercolor=colour[i], markerstrokewidth=mrkroutlinewidth, markerstrokecolor=mrkroutlinecolour[i]) # scatter - Xi(t)
            end
        end

        display(p1_updated)
        savefig(p1_updated,filepath_save[1] * "Figp1_updated" * ".pdf")

        # Preallocate space for storing arrays
        interp_points_ri_range_store = zeros(num_to_est,interp_npts)
        interp_nllri_store = zeros(num_to_est,interp_npts)
        threshold_equality = 10^-13 # Threshold difference between an MLE and a bound within which equality is assumed
        for i in 1:num_to_est
            # Interpolate values of ri within the range
            ri_min = lb[i]
            ri_max = ub[i]
            indices_of_MLE = findall(isequal(maximum(nllri_store[i,:])),nllri_store[i,:])
            rmle_equals_rimin = ri_min > rmle[i] - threshold_equality && ri_min < rmle[i] + threshold_equality # Check if ri_min = rmle[i], which may have updated since ri_range_store[i,:] was created
            rmle_equals_rimax = ri_max > rmle[i] - threshold_equality && ri_max < rmle[i] + threshold_equality # Check if ri_max = rmle[i], which may have updated since ri_range_store[i,:] was created
            if rmle_equals_rimin && !rmle_equals_rimax
                index_of_MLE = maximum(indices_of_MLE)
                println("opt1")
                display(ri_range_store[i,:])
                display(nllri_store[i,:])
                display([ri_min,rmle[i],ri_max, index_of_MLE])
                interp_ri = LinearInterpolation(ri_range_store[i,index_of_MLE:end], nllri_store[i,index_of_MLE:end])
            elseif rmle_equals_rimax && !rmle_equals_rimin
                println("opt2")
                index_of_MLE = minimum(indices_of_MLE)
                interp_ri = LinearInterpolation(ri_range_store[i,1:index_of_MLE], nllri_store[i,1:index_of_MLE])
            elseif rmle_equals_rimin && rmle_equals_rimax
                println("opt3")
                index1_of_MLE = minimum(indices_of_MLE)
                index2_of_MLE = maximum(indices_of_MLE)
                interp_ri = LinearInterpolation(ri_range_store[i,[index1_of_MLE,index2_of_MLE]], nllri_store[i,[index1_of_MLE,index2_of_MLE]])
            else
                println("opt4")
                interp_ri = LinearInterpolation(ri_range_store[i,:], nllri_store[i,:])
            end
            interp_points_ri_range = LinRange(ri_min,ri_max,interp_npts)
            interp_nllri = interp_ri(interp_points_ri_range)

            interp_points_ri_range_store[i,:] = interp_points_ri_range
            interp_nllri_store[i,:] = interp_nllri

            param_name = all_param_names[i]

            # Plot the loglikelihood profile of parameter ri
            profile_i=plot(                                                                                
                interp_points_ri_range,interp_nllri,xlim=(ri_min,ri_max),ylim=(ymin_prof,ymax_prof),yticks=yticks_prof,
                xlab=L"%$param_name",ylab=L"\hat{\ell}_{p}",legend=false,lw=5,titlefont=fnt, guidefont=fnt, tickfont=fnt,size=fig_2_across_size,
                legendfont=fnt,linecolor=:deepskyblue3)
            profile_i=hline!([-1.92],lw=2,linecolor=:black,linestyle=:dot)
            profile_i=vline!([rmle[i]],lw=3,linecolor=:red)
            if num_known_parameters>0
                if !isnothing(ICs) && !isnothing(guess_IC) && i <= num_X # If the ICs are being estimated and are known
                    if i <= num_X
                        if !isnothing(r_known)    
                            profile_i=vline!([ic_known[i]],lw=3,linecolor=:rosybrown,linestyle=:dash)
                        end
                    elseif i>num_X && i<=length(r_known) 
                        if !isnothing(r_known)    
                            profile_i=vline!([r_known[i]],lw=3,linecolor=:rosybrown,linestyle=:dash)
                        end
                    else
                        if !isnothing(er_known)
                            m = i-length(r_known)
                            profile_i=vline!([er_known[m]],lw=3,linecolor=:rosybrown,linestyle=:dash)
                        end
                    end
                else # If the ICs are not being estimated
                    if i<=length(r_known) 
                        if !isnothing(r_known)    
                            profile_i=vline!([r_known[i]],lw=3,linecolor=:rosybrown,linestyle=:dash)
                        end
                    else
                        if !isnothing(er_known)
                            m = i-length(r_known)
                            profile_i=vline!([er_known[m]],lw=3,linecolor=:rosybrown,linestyle=:dash)
                        end
                    end
                end
            end

            display(profile_i)
            savefig(profile_i,filepath_save[1] * "Fig_profile" * param_name  * ".pdf")
        end
        
        # Plot residuals at MLE
        data0_MLE_recomputed_for_residuals = odesolver(unique_times_all,model_params_only,ICs,ode_system,num_X,true);           # generate data using known parameter values for plotting
        if size(data0_MLE_recomputed_for_residuals)[1] > num_X
            data0_MLE_recomputed_for_residuals = data0_MLE_recomputed_for_residuals'
        end
        # Initialize the main dictionary
        data_residuals = Dict{Int, Vector{Float64}}() # A dictionary where each key is the index of the state variable and the values are themselves dictionaries
        extrme_per_X = zeros(1,num_X)
        # Loop through state-variables and time points
        for j in 1:num_X
            data_residuals[j] = []                                  # Initialize empty vector for each dataset where each key is the index
            if typeof(times_all) == Dict{Int64, Vector{Float64}}    # If measurements have not been restricted to the same number per time-point   
                times_all_j = times_all[j]
                data_all_j = data_all[j]
                
            else                                                    # Otherwise measurements have been restricted to the same number per time-point
                times_all_j = copy(times_all)
                data_all_j = data_all[j,:]
            end
            for i in eachindex(unique_times_all)
                # Find the measurements at this unique time-point
                unique_tp = unique_times_all[i]
                unique_solution_at_tp = data0_MLE_recomputed_for_residuals[j,i]
                inds = findall(isequal(unique_tp), times_all_j)
                li = LinearIndices(times_all_j)
                cols_at_tp = li[inds]
                all_measurements_at_tp = data_all_j[cols_at_tp]
                if error_type == "Normal"
                    # Calculate additive residuals
                    data_residuals_at_tp = [measurement_at_tp .- unique_solution_at_tp for measurement_at_tp in all_measurements_at_tp]
                elseif error_type == "Lognormal"
                    # Calculate multiplicative residuals
                    data_residuals_at_tp =  [measurement_at_tp ./ unique_solution_at_tp for measurement_at_tp in all_measurements_at_tp]
                end
                # Store the residuals in the nested dictionary
                data_residuals[j] = [data_residuals[j];data_residuals_at_tp]
                # Store the maximum value
                extrme_per_X[1,j] = maximum(abs,data_residuals_at_tp)
            end
        end
        extrme = maximum(abs,extrme_per_X)  # Maximum value of the residuals for plotting
        
        if error_type == "Normal"           # Minimum value for plotting resdisuals
            res_ymin = -extrme-0.3*extrme   
        elseif error_type == "Lognormal"
            res_ymin = 0
        end
        p_residuals = plot(                 # Plot residuals                                                        # Create figure for residuals
            xlab=L"t", ylab=L"\hat{e}_{i}", size=fig_2_across_size,
            legend=false,xlims=(xmin,xmax),ylims=(res_ymin,extrme+0.3*extrme),
            titlefont=fnt,guidefont=fnt,tickfont=fnt,legendfont=fnt)   
        residuals_for_QQplot = []           # Initialise storage of residuals in appropriate format for QQ-plot
        for i in 1:num_X                                                                                            # Plot data
            residuals_for_QQplot = [residuals_for_QQplot;data_residuals[i]]                                         # Add to the residuals array
            if !(typeof(times_all) == Dict{Int64, Vector{Float64}}) # If measurements have been restricted to the same number per time-point
                p_residuals=scatter!(times_all[:],data_residuals[i],markersize =mrkrsize,markercolor=colour[i],markerstrokewidth=mrkroutlinewidth, markerstrokecolor=mrkroutlinecolour[i]) # scatter - Xi(t)    
            else                                                    # Otherwise if measurements have not been restricted to the same number per time-point
                times_all_i = times_all[i]
                #p_residuals=scatter!(times_all_i[:],data_residuals[i],markersize =mrkrsize,markercolor=colour[i],markerstrokewidth=mrkroutlinewidth, markerstrokecolor=mrkroutlinecolour[i]) # scatter - Xi(t)    
                p_residuals=scatter!(times_all[i],data_residuals[i],markersize =mrkrsize,markercolor=colour[i],markerstrokewidth=mrkroutlinewidth, markerstrokecolor=mrkroutlinecolour[i]) # scatter - Xi(t)    
            end
        end   
        residuals_for_QQplot = Float64.(residuals_for_QQplot)   # Convert values to floats for QQ-plot 

        if error_type == "Normal"           # Add a horizontal line to indicate where the residuals should be scattered around
            p_residuals=hline!([0],lw=2,linecolor=:black,linestyle=:dot)      
        elseif error_type == "Lognormal"
            p_residuals=hline!([1],lw=2,linecolor=:black,linestyle=:dot)      
        end   
        savefig(p_residuals,filepath_save[1] * "Figp_residuals"   * ".pdf")
        
        # Plot qq-plot
        if error_type == "Normal"
            p_residuals_qqplot = plot(                                                                              # QQ-plot with Normal percentiles
                qqplot(Normal,residuals_for_QQplot,lw=2,linecolor=:black,linestyle=:dot,markersize = mrkrsize,markercolor=:black),
                legend=false,xlab=L"\mathrm{Normal}",ylab=L"\mathrm{Residuals}",lw=3,titlefont=fnt, guidefont=fnt,size=fig_2_across_size,
                tickfont=fnt,ylims=(-15,15),yticks=[-10,0,10])
            display(p_residuals_qqplot)
            savefig(p_residuals_qqplot,filepath_save[1] * "Figp_residuals_qqplot"   * ".pdf")
        elseif error_type == "Lognormal"
            p_residuals_qqplot_lognormal = plot(                                                                    # QQ-plot with LogNormal percentiles
                qqplot(LogNormal,residuals_for_QQplot,lw=2,linecolor=:black,linestyle=:dot,markersize = mrkrsize,markercolor=:black),size=fig_2_across_size,
                legend=false,xlab=L"\mathrm{LogNormal}",ylab=L"\hat{e}_{i}",lw=3,titlefont=fnt, guidefont=fnt, tickfont=fnt)
            display(p_residuals_qqplot_lognormal)
            savefig(p_residuals_qqplot_lognormal,filepath_save[1] * "Figp_residuals_qqplot_lognormal"   * ".pdf")
        end

    #######################################################################################
        ### Compute the bounds of confidence interval
        
        # Preallocate space for lower and upper bounds of confidence intervals
        lb_CI = zeros(1,num_to_est)
        ub_CI = zeros(1,num_to_est)

        # Create column names for DataFrame
        column_names = [Symbol("$ii MLE") for ii in all_param_names]
        column_names = vcat(column_names, [Symbol("lb CI $ii") for ii in all_param_names])
        column_names = vcat(column_names, [Symbol("ub CI $ii") for ii in all_param_names])

        # Initialise DataFrame with one row and num_to_est columns (all filled with 'missing' initially)
        is_defined = 0
        if @isdefined(df_MLEBoundsAll) == 0 || isnothing(df_MLEBoundsAll)
            println("MLE TABLE NOT DEFINE")
            global df_MLEBoundsAll = DataFrame(column_names .=> fill(missing, length(column_names)))
        else
            println("MLE TABLE IS DEFINED")
            is_defined = 1
            global df_MLEBoundsAll_thisrow = DataFrame(column_names .=> fill(missing, length(column_names)))
        end

        # Parameter by parameter, fill up row of statistics
        for i in 1:num_to_est
            param_name = all_param_names[i]

            (lb_CI_ri,ub_CI_ri) = funcInterpCI(rmle[i],interp_points_ri_range_store[i,:],interp_nllri_store[i,:],TH)
            lb_CI[i]=lb_CI_ri
            ub_CI[i]=ub_CI_ri
            println("CI for parameter " * param_name * ":")
            println(round(lb_CI_ri; digits = 4))
            println(round(ub_CI_ri; digits = 4))

            if is_defined == 0
                df_MLEBoundsAll[!, Symbol(param_name * " MLE")] .= rmle[i]
                df_MLEBoundsAll[!, Symbol("lb CI " * param_name)] .= lb_CI[i]
                df_MLEBoundsAll[!, Symbol("ub CI " * param_name)] .= ub_CI[i]
            elseif is_defined == 1
                df_MLEBoundsAll_thisrow[!, Symbol(param_name * " MLE")] .= rmle[i]
                df_MLEBoundsAll_thisrow[!, Symbol("lb CI " * param_name)] .= lb_CI[i]
                df_MLEBoundsAll_thisrow[!, Symbol("ub CI " * param_name)] .= ub_CI[i]
                append!(df_MLEBoundsAll,df_MLEBoundsAll_thisrow)
            end
        end
        
        # Export MLE and bounds to csv (one file for all data) --- MLE ONLY 
        CSV.write(filepath_save[1] * "MLEBoundsALL.csv", df_MLEBoundsAll)
        
    #######################################################################################
        ### Confidence sets for model solutions

        # Define min and max y-values given the confidence intervals around the data
        ymin_datCI = copy(ymin_dat)                                                                     # For data with the CI
        ymax_datCI = copy(ymax_dat)
        for i in 1:num_to_est
            for j in 1:num_X
                yij_min = minimum(min_ri_store[i,j,:])      # Min value of best fit + confidence ribbon
                yij_max = maximum(max_ri_store[i,j,:])      # Max value of best fit + confidence ribbon
                ymin_datCI = min(yij_min, ymin_datCI)       # Compare against previous min
                ymax_datCI = max(yij_max, ymax_datCI)       # Compare against previous max
            end
        end
        # Add a bit of leeway in the axes
        ymin_datCI = ymin_datCI - 0.1 * abs(ymin_datCI)           
        ymax_datCI = ymax_datCI + 0.1 * abs(ymax_datCI)          

        for i in 1:num_to_est  
            param_name = all_param_names[i]
            
            # ri - plot model simulated at MLE, scatter data, and prediction interval.
            confmodel_i = plot(                                                                         # Create figure for parameter ri
                xlab=L"t", ylab=L"%$Xi_string",                                                                               
                legend=false,xlims=(xmin,xmax),ylims=(ymin_datCI,ymax_datCI),size=fig_2_across_size,
                titlefont=fnt,guidefont=fnt,tickfont=fnt,legendfont=fnt)   

            for j in 1:num_X                                                                            # Plot best fit trajectories
                confmodel_i=plot!(time_smooth,data0_smooth_MLE_recomputed[j,:],lw=3,linecolor=colour[j])# solid - Xi(t)
            end

            for j in 1:num_X                                                                                    # Plot data
                if restricted_msrmnts_per_tp    # If measurements have been restricted to the same number per time-point
                    confmodel_i=scatter!(times_specific[:],data_specific[j,:],markersize = mrkrsize,markercolor=colour[j], markerstrokewidth=mrkroutlinewidth, markerstrokecolor=mrkroutlinecolour[j]) # scatter Xi(t)
                else                            # Otherwise if measurements have not been restricted to the same number per time-point
                    times_specific_j = times_specific[j]
                    data_specific_j = data_specific[j]
                    #confmodel_i=scatter!(times_specific_j[:],data_specific_j[:],markersize = mrkrsize,markercolor=colour[j], markerstrokewidth=mrkroutlinewidth, markerstrokecolor=mrkroutlinecolour[j]) # scatter Xi(t)
                    confmodel_i=scatter!(times_specific[j],data_specific[j],markersize = mrkrsize,markercolor=colour[j], markerstrokewidth=mrkroutlinewidth, markerstrokecolor=mrkroutlinecolour[j]) # scatter Xi(t)
                end
            end

            for j in 1:num_X                                                                            # Plot the confidence ribbons
                confmodel_i=plot!(
                    time_smooth,data0_smooth_MLE_recomputed[j,:],w=0,c=colour[j],
                    ribbon=(
                        data0_smooth_MLE_recomputed[j,:].-min_ri_store[i,j,:],                          # Shaded region beneath line of best fit
                        max_ri_store[i,j,:].-data0_smooth_MLE_recomputed[j,:]                           # Shaded region above line of best fit
                    ),fillalpha=.2,lw=0,linecolor=colour[j])
            end

            display(confmodel_i)
            savefig(confmodel_i,filepath_save[1] * "Fig_confmodel" * param_name  * ".pdf")
        end
        
        # union
        max_overall = -Inf * ones(num_X,length(time_smooth))
        min_overall = Inf * ones(num_X,length(time_smooth))
        for i in 1:num_X
            for j in eachindex(time_smooth)
                max_overall[i,j]=maximum(max_ri_store[:,i,j])
                min_overall[i,j]=minimum(min_ri_store[:,i,j])
            end
        end

        # Plot model simulated at MLE, scatter data, and prediction interval for union of parameters.
        confmodel_u = plot(                                                                         # Create figure for union of parameters
            xlab=L"t", ylab=L"%$Xi_string",                                                                                  
            legend=false,xlims=(xmin,xmax),ylims=(ymin_datCI,ymax_datCI),size=fig_2_across_size,
            titlefont=fnt,guidefont=fnt,tickfont=fnt,legendfont=fnt)   

        for i in 1:num_X                                                                            # Plot best fit trajectories
            confmodel_u=plot!(time_smooth,data0_smooth_MLE_recomputed[i,:],lw=3,linecolor=colour[i])# solid - Xi(t)
        end

        for i in 1:num_X                                                                                    # Plot data
            if restricted_msrmnts_per_tp    # If measurements have been restricted to the same number per time-point
                confmodel_u=scatter!(times_specific[:],data_specific[i,:],markersize = mrkrsize,markercolor=colour[i], markerstrokewidth=mrkroutlinewidth, markerstrokecolor=mrkroutlinecolour[i]) # scatter Xi(t)
            else                            # Otherwise if measurements have not been restricted to the same number per time-point
                times_specific_i = times_specific[i]
                data_specific_i = data_specific[i]
                #confmodel_u=scatter!(times_specific_i[:],data_specific_i[:],markersize = mrkrsize,markercolor=colour[i], markerstrokewidth=mrkroutlinewidth, markerstrokecolor=mrkroutlinecolour[i]) # scatter Xi(t)
                confmodel_u=scatter!(times_specific[i],data_specific[i],markersize = mrkrsize,markercolor=colour[i], markerstrokewidth=mrkroutlinewidth, markerstrokecolor=mrkroutlinecolour[i]) # scatter Xi(t)
            end
        end

        for i in 1:num_X                                                                            # Plot the confidence ribbons
            confmodel_u=plot!(
                time_smooth,data0_smooth_MLE_recomputed[i,:],w=0,c=colour[i],
                ribbon=(
                    data0_smooth_MLE_recomputed[i,:].-min_overall[i,:],                             # Shaded region beneath line of best fit
                    max_overall[i,:].-data0_smooth_MLE_recomputed[i,:]                              # Shaded region above line of best fit
                ),fillalpha=.2,lw=0,linecolor=colour[i])
        end

        display(confmodel_u)
        savefig(confmodel_u,filepath_save[1] * "Fig_confmodelu"   * ".pdf")
    
    #######################################################################################
        ### Plot difference of confidence set for model solutions and the model solution at MLE 
        
        # Define min and max y-values given the confidence intervals
        ymin_CI1 = 100                                      # For the CI alone
        ymax_CI1 = 0
        for i in 1:num_to_est
            for j in 1:num_X
                CIij_min = minimum(min_ri_store[i,j,:].-data0_smooth_MLE_recomputed[j,:])     # Min value of confidence ribbon
                CIij_max = maximum(max_ri_store[i,j,:].-data0_smooth_MLE_recomputed[j,:])     # Max value of confidence ribbon
                ymin_CI1 = min(ymin_CI1, CIij_min)          # Compare against previous min
                ymax_CI1 = max(ymax_CI1, CIij_max)          # Comapre against previous max
            end
        end
        # Add a bit of leeway in the axes         
        ymin_CI1 = ymin_CI1 - 0.1 * abs(ymin_CI1)
        ymax_CI1 = ymax_CI1 + 0.1 * abs(ymax_CI1)

        # ri
        for i = 1:num_to_est
            param_name = all_param_names[i]
            confmodeldiff_i = plot(                                                                     # Create figure for parameter ri
                xlab=L"t", ylab=L"C_{y,0.95}^{%$param_name} - y(\hat{\theta})",                                                                                   
                legend=false,xlims=(xmin,xmax),ylims=(ymin_CI1,ymax_CI1),size=fig_2_across_size,
                titlefont=fnt,guidefont=fnt,tickfont=fnt,legendfont=fnt)   
            for j = 1:num_X
                confmodeldiff_i = plot!(                                                                # Plot the confidence ribbons
                    time_smooth,0 .* data0_smooth_MLE_recomputed[j,:],w=0,c=colour[j],
                    ribbon=(
                        data0_smooth_MLE_recomputed[j,:].-min_ri_store[i,j,:],                          # Shaded region beneath zero
                        max_ri_store[i,j,:].-data0_smooth_MLE_recomputed[j,:]                           # Shaded region above zero
                    ),fillalpha=.2,lw=0,linecolor=colour[j])
                if j==num_X                                                                             # Plot horizontal line at y=0
                    confmodeldiff_i = plot!(time_smooth, 0 .* data0_smooth_MLE_recomputed[j,:], lw=1, c=:black, linestyle=:dash)
                end
            end
            
            display(confmodeldiff_i)
            savefig(confmodeldiff_i,filepath_save[1] * "Fig_confmodeldiff" * param_name * ".pdf")
        end
        
        # union
        confmodeldiff_u = plot(                                                                         # Create figure for union of parameters
            xlab=L"t", ylab=L"C_{y,0.95} - y(\hat{\theta})",                                                                                 
            legend=false,xlims=(xmin,xmax),ylims=(ymin_CI1,ymax_CI1),size=fig_2_across_size,
            titlefont=fnt,guidefont=fnt,tickfont=fnt,legendfont=fnt) 
        for  i in 1:num_X
            confmodeldiff_u=plot!(                                                                      # Plot the confidence ribbons
                time_smooth, 0 .* data0_smooth_MLE_recomputed[i,:],w=0,c=colour[i],
                ribbon=(
                    data0_smooth_MLE_recomputed[i,:].-min_overall[i,:],                                 # Shaded region beneath line of best fit
                    max_overall[i,:].-data0_smooth_MLE_recomputed[i,:]                                  # Shaded region above line of best fit
                ),fillalpha=.2,lw=0,linecolor=colour[i])
            if i == num_X                                                                               # Plot horizontal line at y=0
                confmodeldiff_u=plot!(time_smooth,0 .* data0_smooth_MLE_recomputed[i,:],lw=1,c=:black,linestyle=:dash)
            end
        end
            
        display(confmodeldiff_u)
        savefig(confmodeldiff_u,filepath_save[1] * "Fig_confmodeldiffu"   * ".pdf")

    #######################################################################################
        ### Confidence sets for model realisations
            
        # Define min and max y-values given the confidence intervals around the realisations
        ymin_realCI = copy(ymin_dat)                                                                    # For data with the CI
        ymax_realCI = copy(ymax_dat)
        for i in 1:num_to_est
            for j in 1:num_X
                yij_min = minimum(min_realisations_ri_store[i,j,:])   # Min value of best fit + confidence ribbon
                yij_max = maximum(max_realisations_ri_store[i,j,:])   # Max value of best fit + confidence ribbon       
                ymin_realCI = min(yij_min, ymin_realCI)                                                 # Compare against previous min
                ymax_realCI = max(yij_max, ymax_realCI)                                                 # Compare against previous max
            end
        end
        # Add a bit of leeway in the axes
        ymin_realCI = ymin_realCI - 0.1 * abs(ymin_realCI)           
        ymax_realCI = ymax_realCI + 0.1 * abs(ymax_realCI)          

        for i in 1:num_to_est  
            param_name = all_param_names[i]
            
            # ri - plot model simulated at MLE, scatter data, and confidence set for model realisations.
            confreal_i = plot(                                                                          # Create figure for parameter ri
                xlab=L"t", ylab=L"%$Xi_string",                                                                        
                legend=false,xlims=(xmin,xmax),ylims=(ymin_realCI,ymax_realCI),size=fig_2_across_size,
                titlefont=fnt,guidefont=fnt,tickfont=fnt,legendfont=fnt)   

            for j in 1:num_X                                                                            # Plot best fit trajectories
                confreal_i=plot!(time_smooth,data0_smooth_MLE_recomputed[j,:],lw=3,linecolor=colour[j]) # solid - Xi(t)
            end

            for j in 1:num_X                                                                                    # Plot data
                if restricted_msrmnts_per_tp    # If measurements have been restricted to the same number per time-point
                    confreal_i=scatter!(times_specific[:],data_specific[j,:],markersize = mrkrsize,markercolor=colour[j], markerstrokewidth=mrkroutlinewidth, markerstrokecolor=mrkroutlinecolour[j]) # scatter Xi(t)
                else                            # Otherwise if measurements have not been restricted to the same number per time-point
                    times_specific_j = times_specific[j]
                    data_specific_j = data_specific[j]
                    #confreal_i=scatter!(times_specific_j[:],data_specific_j[:],markersize = mrkrsize,markercolor=colour[j], markerstrokewidth=mrkroutlinewidth, markerstrokecolor=mrkroutlinecolour[j]) # scatter Xi(t)
                    confreal_i=scatter!(times_specific[j],data_specific[j],markersize = mrkrsize,markercolor=colour[j], markerstrokewidth=mrkroutlinewidth, markerstrokecolor=mrkroutlinecolour[j]) # scatter Xi(t)
                end
            end

            for j in 1:num_X                                                                            # Plot the confidence ribbons
                confreal_i=plot!(
                    time_smooth,data0_smooth_MLE_recomputed[j,:],w=0,c=colour[j],
                    ribbon=(
                        data0_smooth_MLE_recomputed[j,:].-min_realisations_ri_store[i,j,:],             # Shaded region beneath line of best fit
                        max_realisations_ri_store[i,j,:].-data0_smooth_MLE_recomputed[j,:]              # Shaded region above line of best fit
                    ),fillalpha=.2)
            end

            display(confreal_i)
            savefig(confreal_i,filepath_save[1] * "Fig_confreal" * param_name  * ".pdf")
        end
        
        # union
        max_realisations_overall = -Inf * ones(num_X,length(time_smooth))
        min_realisations_overall = Inf * ones(num_X,length(time_smooth))
        for i in 1:num_X
            for j in eachindex(time_smooth)
                max_realisations_overall[i,j]=maximum(max_realisations_ri_store[:,i,j])
                min_realisations_overall[i,j]=minimum(min_realisations_ri_store[:,i,j])
            end
        end

        # Plot model simulated at MLE, scatter data, and confidence set for model realisations for union of parameters.
        confreal_u = plot(                                                                          # Create figure for union of parameters
            xlab=L"t", ylab=L"%$Xi_string",                                                                               
            legend=false,xlims=(xmin,xmax),ylims=(ymin_realCI,ymax_realCI),size=fig_2_across_size,
            titlefont=fnt,guidefont=fnt,tickfont=fnt,legendfont=fnt)   

        for i in 1:num_X                                                                            # Plot best fit trajectories
            confreal_u=plot!(time_smooth,data0_smooth_MLE_recomputed[i,:],lw=3,linecolor=colour[i]) # solid - Xi(t)
        end

        for i in 1:num_X                                                                                    # Plot data
            if restricted_msrmnts_per_tp    # If measurements have been restricted to the same number per time-point
                confreal_u=scatter!(times_specific[:],data_specific[i,:],markersize = mrkrsize,markercolor=colour[i], markerstrokewidth=mrkroutlinewidth, markerstrokecolor=mrkroutlinecolour[i]) # scatter Xi(t)
            else                            # Otherwise if measurements have not been restricted to the same number per time-point
                times_specific_i = times_specific[i]
                data_specific_i = data_specific[i]
                #confreal_u=scatter!(times_specific_i[:],data_specific_i[:],markersize = mrkrsize,markercolor=colour[i], markerstrokewidth=mrkroutlinewidth, markerstrokecolor=mrkroutlinecolour[i]) # scatter Xi(t)
                confreal_u=scatter!(times_specific[i],data_specific[i],markersize = mrkrsize,markercolor=colour[i], markerstrokewidth=mrkroutlinewidth, markerstrokecolor=mrkroutlinecolour[i]) # scatter Xi(t)
            end
        end

        for i in 1:num_X                                                                            # Plot the confidence ribbons
            confreal_u=plot!(
                time_smooth,data0_smooth_MLE_recomputed[i,:],w=0,c=colour[i],
                ribbon=(
                    data0_smooth_MLE_recomputed[i,:].-min_realisations_overall[i,:],                # Shaded region beneath line of best fit
                    max_realisations_overall[i,:].-data0_smooth_MLE_recomputed[i,:]                 # Shaded region above line of best fit
                ),fillalpha=.2)
        end

        display(confreal_u)
        savefig(confreal_u,filepath_save[1] * "Fig_confrealu"   * ".pdf")

    #######################################################################################
        # Plot difference of confidence sets for data realisations and model solution at MLE

        # Define min and max y-values given the confidence intervals
        ymin_CI2 = 100                                      # For the CI alone
        ymax_CI2 = 0
        for i in 1:num_to_est
            for j in 1:num_X
                CIij_min = minimum(min_realisations_ri_store[i,j,:].-data0_smooth_MLE_recomputed[j,:])  # Min value of confidence ribbon
                CIij_max = maximum(max_realisations_ri_store[i,j,:].-data0_smooth_MLE_recomputed[j,:])  # Max value of confidence ribbon
                ymin_CI2 = min(ymin_CI2, CIij_min)                                                      # Compare against previous min
                ymax_CI2 = max(ymax_CI2, CIij_max)                                                      # Comapre against previous max
            end
        end
        # Add a bit of leeway in the axes         
        ymin_CI2 = ymin_CI2 - 0.1 * abs(ymin_CI2)
        ymax_CI2 = ymax_CI2 + 0.1 * abs(ymax_CI2)

        # ri
        for i = 1:num_to_est
            param_name = all_param_names[i]
            confrealdiff_i = plot(                                                                     # Create figure for parameter ri
                xlab=L"t", ylab=L"C_{z_{i},0.95}^{%$param_name} - y(\hat{\theta})",                                                                                    
                legend=false,xlims=(xmin,xmax),ylims=(ymin_CI2,ymax_CI2),size=fig_2_across_size,
                titlefont=fnt,guidefont=fnt,tickfont=fnt,legendfont=fnt)   
            for j = 1:num_X
                confrealdiff_i = plot!(                                                                # Plot the confidence ribbons
                    time_smooth,0 .* data0_smooth_MLE_recomputed[j,:],w=0,c=colour[j],
                    ribbon=(
                        data0_smooth_MLE_recomputed[j,:].-min_realisations_ri_store[i,j,:],                          # Shaded region beneath zero
                        max_realisations_ri_store[i,j,:].-data0_smooth_MLE_recomputed[j,:]                           # Shaded region above zero
                    ),fillalpha=.2,lw=0,linecolor=colour[j])
                if j==num_X                                                                             # Plot horizontal line at y=0
                    confrealdiff_i = plot!(time_smooth, 0 .* data0_smooth_MLE_recomputed[j,:], lw=1, c=:black, linestyle=:dash)
                end
            end
            
            display(confrealdiff_i)
            savefig(confrealdiff_i,filepath_save[1] * "Fig_confrealdiff" * param_name * ".pdf")
        end
        
        # union
        confrealdiff_u = plot(                                                                         # Create figure for union of parameters
            xlab=L"t", ylab=L"C_{z_{i},0.95} - y(\hat{\theta})",                                                                                
            legend=false,xlims=(xmin,xmax),ylims=(ymin_CI2,ymax_CI2),size=fig_2_across_size,
            titlefont=fnt,guidefont=fnt,tickfont=fnt,legendfont=fnt) 
        for  i in 1:num_X
            confrealdiff_u=plot!(                                                                      # Plot the confidence ribbons
                time_smooth, 0 .* data0_smooth_MLE_recomputed[i,:],w=0,c=colour[i],
                ribbon=(
                    data0_smooth_MLE_recomputed[i,:].-min_realisations_overall[i,:],                                 # Shaded region beneath line of best fit
                    max_realisations_overall[i,:].-data0_smooth_MLE_recomputed[i,:]                                  # Shaded region above line of best fit
                ),fillalpha=.2,lw=0,linecolor=colour[i])
            if i == num_X                                                                               # Plot horizontal line at y=0
                confrealdiff_u=plot!(time_smooth,0 .* data0_smooth_MLE_recomputed[i,:],lw=1,c=:black,linestyle=:dash)
            end
        end
            
        display(confrealdiff_u)
        savefig(confrealdiff_u,filepath_save[1] * "Fig_confrealdiffu"   * ".pdf")

    #######################################################################################
        # If the following are defined, delete their definitions
        if @isdefined(df_MLEBoundsAll)
            df_MLEBoundsAll = nothing
        elseif @isdefined(df_MLEBoundsAll_thisrow)
            df_MLEBoundsAll_thisrow = nothing
        end
end
