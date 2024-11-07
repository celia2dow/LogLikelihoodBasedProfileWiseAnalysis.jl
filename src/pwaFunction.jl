function pwaFunction(time, X_array, ode_system, model_params, error_type, noise_params, seed_num, data, npts)
    # Implementing an input measurement error models with an input mechanistic model of ODEs or PDEs
    # in a likelihood-based framework for estimation, identifiability analysis, 
    # and prediction. Note that the code is ideal for visualising a model with up to 6 state variables.
    # If more state variables are desired, add more colours to the colours array. The code will still
    # work without this update, however.

    # Dowling (2024)
    # Adapted from Murphy et al. (2023) 
    
    #######################################################################################
    ### Initialisation including plot options, and filepaths to save outputs
    #pyplot()                                                                                    # plot options
    fnt = Plots.font("sans-serif", 28)                                                          # plot options
    global cur_colors = palette(:default)                                                       # plot options
    isdir(pwd() * "/ProfileWiseAnalysisOutput/") || mkdir(pwd() * "/ProfileWiseAnalysisOutput") # make folder to save figures if doesnt already exist
    filepath_save = [pwd() * "/ProfileWiseAnalysisOutput/"]                                     # location to save figures "//ProfileWiseAnalysisOutput//"]

    #######################################################################################
    ### Generate time-points for plotting known values and best-fits
    t_max = time[end]
    time_smooth =  LinRange(0.0,t_max*1.1, 201);

    #######################################################################################
    ### Set random seed
    Random.seed!(seed_num)

    #######################################################################################
    ### Define values, names, initial conditions, and guesses
    ICs = X_array.initial_values                # Initial conditions
    num_x = length(ICs);                        # Number of state variables
    Xi_string = join(X_array.names)             # String of state variables
    guess_r = model_params.initial_values       # Initial guesses for model parameters
    guess_er = noise_params.initial_values      # Initial guesses for noise parameters
    θG = [guess_r guess_er]                     # Combined initial guesses
    lb_r = model_params.lower_bounds            # Lower bounds
    lb_er = noise_params.lower_bounds
    lb = [lb_r lb_er];        
    ub_r = model_params.upper_bounds            # Upper bounds
    ub_er = noise_params.upper_bounds
    ub = [ub_r ub_er];        
    num_r = length(guess_r) + length(guess_er)  # Number of model and noise parameters combined
    if !isnothing(model_params.true_values)
        r = model_params.true_values            # Model parameter values
    else
        r = nothing
    end
    if !isnothing(noise_params.true_values)
        er = noise_params.true_values           # Error model parameter values
    else
        er = nothing
    end

    #######################################################################################
    ### Confidence thresholds
    TH=-1.921; #95% confidence interval threshold (for confidence interval for model parameters, and confidence set for model solutions)
    TH_realisations = -2.51 # 97.5 confidence interval threshold (for model realisations)
    THalpha = 0.025; # 97.5% data realisations interval threshold

    #######################################################################################
    ### Plot presentation preferences
    # Plot limits
    xmin=-0.1*t_max                        
    xmax= t_max+0.1*t_max
    ymin_dat= -0.1*maximum(data)
    ymax_dat= maximum(data) + 0.1*maximum(data)
    # Define colours - FEEL FREE TO ADD CUSTOM COLOURS HERE
    colour1 = :green2                           # Bright green
    colour2 = :magenta                          # Magenta
    colour3 = RGB(0.8, 0.7, 0.2);               # Mustard yellow
    colour4 = RGB(0.0, 1.0, 1.0);               # Cyan
    colour5 = RGB(0.90196, 0.90196, 0.98039);   # Lavender
    colour6 = RGB(1.0, 0.8, 0.6);               # Peach
    colour = [colour1, colour2, colour3, colour4, colour5, colour6]
    # Define more colours if necessary
    if num_x>length(colour)
        for i in 1:num_x-length(colour)
            new_colour = RGB(rand(), rand(), rand())   # Create a random new colour
            colour = [colour, new_colour]             # Add it to the array of colours
        end
    end

    #######################################################################################
    ### Generate smooth data/true solution if true parameter parameter values are known
    if !isnothing(r)
        data0_smooth = odesolver(time_smooth,r,ICs,ode_system)
    end

    #######################################################################################
    ### Plot data - scatter plot

    # Plot data against true solution (if available)
    p0scatter = plot(                                                                                   # Create figure
        xlab=L"t", ylab=L"%$Xi_string",
        legend=false,xlims=(xmin,xmax),ylims=(ymin_dat,ymax_dat),
        guidefont=fnt,tickfont=fnt,legendfont=fnt)   
    if !isnothing(r)                                                                                    # If the true parameter values are known
        for i in 1:num_x                                                                                # Plot true trajectories
            p0scatter=plot!(time_smooth,data0_smooth[i,:],lw=3,linecolor=colour[i])                     # solid - Xi(t)
        end
    end
    for i in 1:num_x                                                                                    # Plot data
        p0scatter=scatter!(time,data[i,:],markersize = 8,markercolor=colour[i],markerstrokewidth=0)     # scatter - Xi(t)
    end
    display(p0scatter)
    savefig(p0scatter,filepath_save[1] * "Fig0scatter" * ".pdf")
    savefig(p0scatter,filepath_save[1] * "Fig0scatter" * ".png")

    #######################################################################################
    ### Parameter identifiability analysis - MLE, profiles
        
        # MLE

        # MLE optimisation 
        function funMLE(unknowns)
            # Define function whose inputs are the unknown parameters, and whose output is 
            # that of sum_loglikelihood
            return sum_loglikelihood(time,unknowns,ICs,ode_system,data,error_type)
        end
        (xopt,fopt)  = optimise(funMLE,θG,lb,ub) 

        # Storing MLE 
        rmle = [0 for i in 1:num_r];    # Initialise array for storing the running MLE estimates of parameters
        global fmle=fopt
        for i in 1:num_r
            global rmle[i] = xopt[i]
        end

        data0_smooth_MLE = model(time_smooth,rmle,ICs,ode_system); # generate data using MLE parameter values for plotting

        # Plot model simulated at MLE and data
        p1 = plot(                                                                                   # Create figure
            xlab=L"t", ylab=L"%$Xi_string",
            legend=false,xlims=(xmin,xmax),ylims=(ymin_dat,ymax_dat),
            guidefont=fnt,tickfont=fnt,legendfont=fnt)   
        for i in 1:num_x                                                                             # Plot true trajectories
            p1=plot!(time_smooth,data0_smooth_MLE[i,:],lw=3,linecolor=colour[i])                     # solid - Xi(t)
        end
        for i in 1:num_x                                                                             # Plot data
            p1=scatter!(time,data[i,:],markersize = 8,markercolor=colour[i],markerstrokewidth=0)     # scatter - Xi(t)
        end
        display(p1)
        savefig(p1,filepath_save[1] * "Fig1" * ".pdf")
        savefig(p1,filepath_save[1] * "Fig1" * ".png")

    #######################################################################################
        # Profiling (using the MLE as the first guess at each point)
        for i in 1:num_r
            ri_min = lb[i]
            ri_max = ub[i]
            ri_range_lower = reverse(LinRange(ri_min,rmle[i],npts))
            ri_range_upper = reverse(LinRange(rmle[i] + (ri_max-rmle[i])/npts,ri_max,npts))

            not_ri_range_lower = zeros(num_r-1,npts)
            llri_lower = zeros(npts)
            predict_ri_lower=zeros(num_r-1,length(time_smooth),npts)
            predict_ri_realisations_lower_lq=zeros(num_r-1,length(time_smooth),npts)
            predict_ri_realisations_lower_uq=zeros(num_r-1,length(time_smooth),npts)

            not_ri_range_upper = zeros(num_r-1,npts)
            llri_upper = zeros(npts)
            predict_ri_upper=zeros(num_r-1,length(time_smooth),npts)
            predict_ri_realisations_upper_lq=zeros(num_r-1,length(time_smooth),npts)
            predict_ri_realisations_upper_uq=zeros(num_r-1,length(time_smooth),npts)

            # Lower and upper bounds for parameters that are not ri
            lb_i = lb; lb_i = deleteat!(lb,i);
            ub_i = ub; ub_i = deleteat!(ub,i);

            # Start at mle and increase parameter (upper)
            for j in 1:npts
                function fun(aa)
                    # Create an array of the unknown parameters
                    params_combined = aa[1]
                    for k = 2:num_r-1
                        params_combined = [params_combined, aa[k]]
                    end
                    # Add the single known parameter in the array in the appropriate order
                    params_combined = insert!(params_combined, i, ri_range_upper[j])
                    return sum_loglikelihood(time,params_combined,ICs,ode_system,data,error_type)
                end
                local θG_i=rmle; θG_i = deleteat!(θG_i,i)
                local (xo,fo)=optimise(fun,θG_i,lb_i,ub_i)
                not_ri_range_upper[:,j]=xo[:]
                llri_upper[j]=fo[1]
                model_params_only = not_ri_range_upper[:,j]
                if error_type == "Normal" || error_type == "Lognormal" && i != num_r
                    model_params_only = model_params_only[1:end-1]
                    model_params_only = insert!(model_params_only,i,ri_range_upper[j])
                end
                modelmean = odesolver(time_smooth,model_params_only,ICs,ode_system)
                predict_r1_upper[:,:,j]=modelmean; 

                loop_data_dists=[LogNormal(0,nr1range_upper[2,j]) for mi in modelmean]; 
                predict_r1_realisations_upper_lq[:,:,j]= [quantile(data_dist, THalpha/2) for data_dist in loop_data_dists].*modelmean;
                predict_r1_realisations_upper_uq[:,:,j]= [quantile(data_dist, 1-THalpha/2) for data_dist in loop_data_dists].*modelmean;

                if fo > fmle
                    println("new MLE r1 upper")
                    global fmle=fo
                    global r1mle=r1range_upper[j]
                    global r2mle=nr1range_upper[1,j]
                    global Sdmle=nr1range_upper[2,j]
                end
            end
        end
        
        
        # start at mle and increase parameter (upper)
        for j in 1:nptss
            function fun(aa)
                return error(data,[r1range_upper[j],aa[1],aa[2]])
            end
            lb1=[lb[2],lb[3]];
            ub1=[ub[2],ub[3]];
            local θG1=[r2mle,Sdmle]
            local (xo,fo)=optimise(fun,θG1,lb1,ub1)
            nr1range_upper[:,j]=xo[:]
            llr1_upper[j]=fo[1]

            modelmean = model(t_smooth,[r1range_upper[j],nr1range_upper[1,j]]);
            predict_r1_upper[:,:,j]=modelmean;

            loop_data_dists=[LogNormal(0,nr1range_upper[2,j]) for mi in modelmean]; 
            predict_r1_realisations_upper_lq[:,:,j]= [quantile(data_dist, THalpha/2) for data_dist in loop_data_dists].*modelmean;
            predict_r1_realisations_upper_uq[:,:,j]= [quantile(data_dist, 1-THalpha/2) for data_dist in loop_data_dists].*modelmean;

            if fo > fmle
                println("new MLE r1 upper")
                global fmle=fo
                global r1mle=r1range_upper[j]
                global r2mle=nr1range_upper[1,j]
                global Sdmle=nr1range_upper[2,j]
            end
        end
        
        # start at mle and decrease parameter (lower)
        for j in 1:nptss
            function fun1a(aa)
                return error(data,[r1range_lower[j],aa[1],aa[2]])
            end
            lb1=[lb[2],lb[3]];
            ub1=[ub[2],ub[3]];
            local θG1=[r2mle, Sdmle]
            local (xo,fo)=optimise(fun1a,θG1,lb1,ub1)
            nr1range_lower[:,j]=xo[:]
            llr1_lower[j]=fo[1]
            
            modelmean = model(t_smooth,[r1range_lower[j],nr1range_lower[1,j]]);
            predict_r1_lower[:,:,j]=modelmean;

            loop_data_dists=[LogNormal(0,nr1range_lower[2,j]) for mi in modelmean]; 
            predict_r1_realisations_lower_lq[:,:,j]= [quantile(data_dist, THalpha/2) for data_dist in loop_data_dists].*modelmean;
            predict_r1_realisations_lower_uq[:,:,j]= [quantile(data_dist, 1-THalpha/2) for data_dist in loop_data_dists].*modelmean;

            if fo > fmle
                println("new MLE r2 lower")
                global fmle = fo
                global r1mle=r1range_lower[j]
                global r2mle=nr1range_lower[1,j]
                global Sdmle=nr1range_lower[2,j]
            end
        end
        
        # combine the lower and upper
        r1range = [reverse(r1range_lower); r1range_upper]
        nr1range = [reverse(nr1range_lower); nr1range_upper ]
        llr1 = [reverse(llr1_lower); llr1_upper] 
        nllr1=llr1.-maximum(llr1);

        max_r1=zeros(2,length(t_smooth))
        min_r1=1000*ones(2,length(t_smooth))
        for j in 1:(nptss)
            if (llr1_lower[j].-maximum(llr1)) >= TH
                for j in 1:length(t_smooth)
                    max_r1[1,j]=max(predict_r1_lower[1,j,j],max_r1[1,j])
                    min_r1[1,j]=min(predict_r1_lower[1,j,j],min_r1[1,j])
                    max_r1[2,j]=max(predict_r1_lower[2,j,j],max_r1[2,j])
                    min_r1[2,j]=min(predict_r1_lower[2,j,j],min_r1[2,j])
                end
            end
        end
        for j in 1:(nptss)
            if (llr1_upper[j].-maximum(llr1)) >= TH
                for j in 1:length(t_smooth)
                    max_r1[1,j]=max(predict_r1_upper[1,j,j],max_r1[1,j])
                    min_r1[1,j]=min(predict_r1_upper[1,j,j],min_r1[1,j])
                    max_r1[2,j]=max(predict_r1_upper[2,j,j],max_r1[2,j])
                    min_r1[2,j]=min(predict_r1_upper[2,j,j],min_r1[2,j])
                end
            end
        end
        # combine the confidence sets for realisations
        max_realisations_r1=zeros(2,length(t_smooth))
        min_realisations_r1=1000*ones(2,length(t_smooth))
        for j in 1:(nptss)
            if (llr1_lower[j].-maximum(llr1)) >= TH_realisations
                for j in 1:length(t_smooth)
                    max_realisations_r1[1,j]=max(predict_r1_realisations_lower_uq[1,j,j],max_realisations_r1[1,j])
                    min_realisations_r1[1,j]=min(predict_r1_realisations_lower_lq[1,j,j],min_realisations_r1[1,j])
                    max_realisations_r1[2,j]=max(predict_r1_realisations_lower_uq[2,j,j],max_realisations_r1[2,j])
                    min_realisations_r1[2,j]=min(predict_r1_realisations_lower_lq[2,j,j],min_realisations_r1[2,j])
                end
            end
        end
        for j in 1:(nptss)
            if (llr1_upper[j].-maximum(llr1)) >= TH_realisations
                for j in 1:length(t_smooth)
                    max_realisations_r1[1,j]=max(predict_r1_realisations_upper_uq[1,j,j],max_realisations_r1[1,j])
                    min_realisations_r1[1,j]=min(predict_r1_realisations_upper_lq[1,j,j],min_realisations_r1[1,j]) 
                    max_realisations_r1[2,j]=max(predict_r1_realisations_upper_uq[2,j,j],max_realisations_r1[2,j])
                    min_realisations_r1[2,j]=min(predict_r1_realisations_upper_lq[2,j,j],min_realisations_r1[2,j])
                end
            end
        end
end
