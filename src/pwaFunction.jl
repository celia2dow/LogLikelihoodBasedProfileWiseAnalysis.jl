function pwaFunction(model_params, noise_params, ICs, times, seed_num)
    # Implementing an input measurement error models with an input mechanistic model of ODEs or PDEs
    # in a likelihood-based framework for estimation, identifiability analysis, 
    # and prediction

    # Dowling (2024)
    # Adapted from Murphy et al. (2023) 
    
    #######################################################################################
    ## Initialisation including plot options, and filepaths to save outputs
    pyplot()                                                                                    # plot options
    fnt = Plots.font("sans-serif", 28)                                                          # plot options
    global cur_colors = palette(:default)                                                       # plot options
    isdir(pwd() * "/ProfileWiseAnalysisOutput/") || mkdir(pwd() * "/ProfileWiseAnalysisOutput") # make folder to save figures if doesnt already exist
    filepath_save = [pwd() * "/ProfileWiseAnalysisOutput/"]                                     # location to save figures "//ProfileWiseAnalysisOutput//"]

    #######################################################################################
    ## Generate time-points for plotting known values and best-fits
    t_max = times[end]
    t_smooth =  LinRange(0.0,t_max*1.1, 201);

    #######################################################################################
    ## Set random seed
    Random.seed(seed_num)

    #######################################################################################
    ### Plot data - scatter plot

    xmin=-0.5;
    xmax= t_max+0.5;
    ymin= 0;
    ymax= 330;

    colour1=:green2
    colour2=:magenta

    # plot data
    p0scatter=scatter(t,data[1,:],xlims=(xmin,xmax),ylims=(ymin,ymax),legend=false,markersize = 12,markercolor=colour2,xlab=L"t",ylab=L"R(t)",lw=3,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,markerstrokewidth=0) # scatter R(t)
    display(p0scatter)
    savefig(p0scatter,filepath_save[1] * "Fig0scatter" * ".pdf")
    savefig(p0scatter,filepath_save[1] * "Fig0scatter" * ".png")

    #######################################################################################
    ### Parameter identifiability analysis - MLE, profiles

        # initial guesses for MLE search
        lambda_g=1.05; rr1=lambda_g; # convert to pre-existing variable names
        zeta_g=1.1; rr2=zeta_g; # convert to pre-existing variable names
        R_d_g=70;  rr3=R_d_g
        Sdinit=8;
        
        
        TH=-1.921; #95% confidence interval threshold (for confidence interval for model parameters, and confidence set for model solutions)
        TH_realisations = -2.51 # 97.5 confidence interval threshold (for model realisations)
        THalpha = 0.025; # 97.5% data realisations interval threshold

        function error(data,a)
            y=zeros(length(t))
            y=model(t,a);
            e=0;
            data_dists=[Normal(mi,a[4]) for mi in model(t,a[1:3])]; # Normal noise model
            e+=sum([loglikelihood(data_dists[i],data[i]) for i in 1:length(data_dists)])
            return e
        end
        
        function fun(a)
            return error(data,a)
        end

        function optimise(fun,θ₀,lb,ub) 
            tomax = (θ,∂θ) -> fun(θ)
            opt = Opt(:LN_NELDERMEAD,length(θ₀))
            opt.max_objective = tomax
            opt.lower_bounds = lb       # Lower bound
            opt.upper_bounds = ub       # Upper bound
            opt.maxtime = 30.0; # maximum time in seconds
            res = optimize(opt,θ₀)
            return res[[2,1]]
        end
    
        #######################################################################################
        # MLE

        
        θG = [rr1,rr2,rr3,Sdinit] # first guess

        # lower and upper bounds for parameter estimation
        lambda_lb = 0.9; r1_lb = lambda_lb; # convert to pre-existing variable names
        lambda_ub = 1.1; r1_ub = lambda_ub; # convert to pre-existing variable names
        zeta_lb = 1; r2_lb = zeta_lb; # convert to pre-existing variable names
        zeta_ub = 2; r2_ub = zeta_ub; # convert to pre-existing variable names
        R_d_lb = 60; r3_lb = R_d_lb; # convert to pre-existing variable names
        R_d_ub = 90; r3_ub = R_d_ub; # convert to pre-existing variable names
        Sd_lb = 1.0
        Sd_ub = 10.0

        lb=[r1_lb,r2_lb,r3_lb,Sd_lb];
        ub=[r1_ub,r2_ub,r3_ub,Sd_ub];

        # MLE optimisation
        (xopt,fopt)  = optimise(fun,θG,lb,ub)

        # storing MLE 
        global fmle=fopt
        global r1mle=xopt[1]
        global r2mle=xopt[2]
        global r3mle=xopt[3]
        global Sdmle=xopt[4]

        data0_smooth_MLE = model(t_smooth,[r1mle,r2mle,r3mle,C0]); # model solution at MLE for plotting
        if size(data0_smooth_MLE)[1] > num_x
            data0_smooth_MLE = data0_smooth_MLE'
        end
        # plot model simulated at MLE and data
        p1=plot(t_smooth,data0_smooth_MLE[1,:],xlab=L"t",ylab=L"R(t)",lw=3,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=colour2) # solid - R(t)
        p1=scatter!(t,data[1,:],xlims=(xmin,xmax),ylims=(ymin,ymax),legend=false,markersize = 8,markercolor=colour2,markerstrokewidth=0) # scatter c2
        
        display(p1)
        savefig(p1,filepath_save[1] * "Fig1" * ".pdf")
        savefig(p1,filepath_save[1] * "Fig1" * ".png")


        #######################################################################################
            # Profiling (using the MLE as the first guess at each point)
        
            nptss=20;

            #Profile r1

            r1min=lb[1]
            r1max=ub[1]
            r1range_lower=reverse(LinRange(r1min,r1mle,nptss))
            r1range_upper=LinRange(r1mle + (r1max-r1mle)/nptss,r1max,nptss)
            
            nr1range_lower=zeros(3,nptss) # NOT r1, so r2, r3 and Sd
            llr1_lower=zeros(nptss)
            nllr1_lower=zeros(nptss) # Not necessary?
            predict_r1_lower=zeros(1,length(t_smooth),nptss)
            predict_r1_realisations_lower_lq=zeros(1,length(t_smooth),nptss)
            predict_r1_realisations_lower_uq=zeros(1,length(t_smooth),nptss)
            
            nr1range_upper=zeros(3,nptss) # NOT r1, so r2, r3 and Sd 
            llr1_upper=zeros(nptss)
            nllr1_upper=zeros(nptss) # Not necessary?
            predict_r1_upper=zeros(1,length(t_smooth),nptss)
            predict_r1_realisations_upper_lq=zeros(1,length(t_smooth),nptss)
            predict_r1_realisations_upper_uq=zeros(1,length(t_smooth),nptss)

            
            # start at mle and increase parameter (upper)
            for i in 1:nptss
                function fun1(aa)
                    return error(data,[r1range_upper[i],aa[1],aa[2],aa[3]]) 
                end
                lb1=[lb[2],lb[3],lb[4]]; 
                ub1=[ub[2],ub[3],ub[4]]; 
                local θG1=[r2mle,r3mle,Sdmle] 
                local (xo,fo)=optimise(fun1,θG1,lb1,ub1)
                nr1range_upper[:,i]=xo[:]
                llr1_upper[i]=fo[1]
                predict_r1_upper[:,:,i]=model(t_smooth,[r1range_upper[i],nr1range_upper[1,i],nr1range_upper[2,i]]) 

                loop_data_dists=[Normal(mi,nr1range_upper[3,i]) for mi in model(t_smooth,[r1range_upper[i],nr1range_upper[1,i],nr1range_upper[2,i]])]; # normal distribution about mean
                predict_r1_realisations_upper_lq[:,:,i]= [quantile(data_dist, THalpha/2) for data_dist in loop_data_dists];
                predict_r1_realisations_upper_uq[:,:,i]= [quantile(data_dist, 1-THalpha/2) for data_dist in loop_data_dists];
                if fo > fmle
                    println("new MLE r1 upper")
                    global fmle=fo
                    global r1mle=r1range_upper[i]
                    global r2mle=nr1range_upper[1,i] 
                    global r3mle=nr1range_upper[2,i] 
                    global Sdmle=nr1range_upper[3,i] 
                end
            end
            
            # start at mle and decrease parameter (lower)
            for i in 1:nptss
                function fun1a(aa)
                    return error(data,[r1range_lower[i],aa[1],aa[2],aa[3]]) 
                end
                lb1=[lb[2],lb[3],lb[4]]; 
                ub1=[ub[2],ub[3],ub[4]]; 
                local θG1=[r2mle, r3mle, Sdmle] 
                local (xo,fo)=optimise(fun1a,θG1,lb1,ub1)
                nr1range_lower[:,i]=xo[:]
                llr1_lower[i]=fo[1]
                predict_r1_lower[:,:,i]=model(t_smooth,[r1range_lower[i],nr1range_lower[1,i],nr1range_lower[2,i]]); 

                loop_data_dists=[Normal(mi,nr1range_lower[3,i]) for mi in model(t_smooth,[r1range_lower[i],nr1range_lower[1,i],nr1range_lower[2,i]])]; # normal distribution about mean
                predict_r1_realisations_lower_lq[:,:,i]= [quantile(data_dist, THalpha/2) for data_dist in loop_data_dists];
                predict_r1_realisations_lower_uq[:,:,i]= [quantile(data_dist, 1-THalpha/2) for data_dist in loop_data_dists];

                if fo > fmle
                    println("new MLE r2 lower")
                    global fmle = fo
                    global r1mle=r1range_lower[i]
                    global r2mle=nr1range_lower[1,i]
                    global r3mle=nr1range_lower[2,i] 
                    global Sdmle=nr1range_lower[3,i] 
                end
            end
            
            # combine the lower and upper
            r1range = [reverse(r1range_lower); r1range_upper]
            nr1range = [reverse(nr1range_lower); nr1range_upper ]
            llr1 = [reverse(llr1_lower); llr1_upper] 
            nllr1=llr1.-maximum(llr1);
            
            max_r1=zeros(1,length(t_smooth))
            min_r1=1000*ones(1,length(t_smooth))
            for i in 1:(nptss)
                if (llr1_lower[i].-maximum(llr1)) >= TH
                    for j in 1:length(t_smooth)
                        max_r1[1,j]=max(predict_r1_lower[1,j,i],max_r1[1,j])
                        min_r1[1,j]=min(predict_r1_lower[1,j,i],min_r1[1,j])
                    end
                end
            end
            for i in 1:(nptss)
                if (llr1_upper[i].-maximum(llr1)) >= TH
                    for j in 1:length(t_smooth)
                        max_r1[1,j]=max(predict_r1_upper[1,j,i],max_r1[1,j])
                        min_r1[1,j]=min(predict_r1_upper[1,j,i],min_r1[1,j]) 
                    end
                end
            end
            # combine the confidence sets for realisations
            max_realisations_r1=zeros(1,length(t_smooth))
            min_realisations_r1=1000*ones(1,length(t_smooth))
            for i in 1:(nptss)
                if (llr1_lower[i].-maximum(llr1)) >= TH_realisations
                    for j in 1:length(t_smooth)
                        max_realisations_r1[1,j]=max(predict_r1_realisations_lower_uq[1,j,i],max_realisations_r1[1,j])
                        min_realisations_r1[1,j]=min(predict_r1_realisations_lower_lq[1,j,i],min_realisations_r1[1,j])
                    end
                end
            end
            for i in 1:(nptss)
                if (llr1_upper[i].-maximum(llr1)) >= TH_realisations
                    for j in 1:length(t_smooth)
                        max_realisations_r1[1,j]=max(predict_r1_realisations_upper_uq[1,j,i],max_realisations_r1[1,j])
                        min_realisations_r1[1,j]=min(predict_r1_realisations_upper_lq[1,j,i],min_realisations_r1[1,j]) 
                    end
                end
            end


            #Profile r2
            r2min=lb[2]
            r2max=ub[2]
            r2range_lower=reverse(LinRange(r2min,r2mle,nptss))
            r2range_upper=LinRange(r2mle + (r2max-r2mle)/nptss,r2max,nptss)
            
            nr2range_lower=zeros(3,nptss) # NOT r2, so r1, r3 and Sd
            llr2_lower=zeros(nptss)
            nllr2_lower=zeros(nptss) # Not necessary?
            predict_r2_lower=zeros(1,length(t_smooth),nptss)
            predict_r2_realisations_lower_lq=zeros(1,length(t_smooth),nptss)
            predict_r2_realisations_lower_uq=zeros(1,length(t_smooth),nptss)
            
            nr2range_upper=zeros(3,nptss) # NOT r2, so r1, r3 and Sd
            llr2_upper=zeros(nptss)
            nllr2_upper=zeros(nptss) # Not necessary?
            predict_r2_upper=zeros(1,length(t_smooth),nptss)
            predict_r2_realisations_upper_lq=zeros(1,length(t_smooth),nptss)
            predict_r2_realisations_upper_uq=zeros(1,length(t_smooth),nptss)
            
            # start at mle and increase parameter (upper)
            for i in 1:nptss
                function fun2(aa)
                    return error(data,[aa[1],r2range_upper[i],aa[2],aa[3]])
                end
                lb1=[lb[1],lb[3],lb[4]];
                ub1=[ub[1],ub[3],ub[4]];
                local θG1=[r1mle, r3mle, Sdmle]
                local (xo,fo)=optimise(fun2,θG1,lb1,ub1)
                nr2range_upper[:,i]=xo[:]
                llr2_upper[i]=fo[1]
                predict_r2_upper[:,:,i]=model(t_smooth,[nr2range_upper[1,i],r2range_upper[i],nr2range_upper[2,i]])

                loop_data_dists=[Normal(mi,nr2range_upper[3,i]) for mi in model(t_smooth,[nr2range_upper[1,i],r2range_upper[i],nr2range_upper[2,i]])]; # normal distribution about mean
                predict_r2_realisations_upper_lq[:,:,i]= [quantile(data_dist, THalpha/2) for data_dist in loop_data_dists];
                predict_r2_realisations_upper_uq[:,:,i]= [quantile(data_dist, 1-THalpha/2) for data_dist in loop_data_dists];

                if fo > fmle
                    println("new MLE r2 upper")
                    global fmle = fo
                    global r1mle=nr2range_upper[1,i]
                    global r2mle=r2range_upper[i]
                    global r3mle=nr2range_upper[2,i]
                    global Sdmle=nr2range_upper[3,i]
                end
            end
            
            # start at mle and decrease parameter (lower)
            for i in 1:nptss
                function fun2a(aa)
                    return error(data,[aa[1],r2range_lower[i],aa[2],aa[3]])
                end
                lb1=[lb[1],lb[3],lb[4]];
                ub1=[ub[1],ub[3],ub[4]];
                local θG1=[r1mle, r3mle, Sdmle]
                local (xo,fo)=optimise(fun2a,θG1,lb1,ub1)
                nr2range_lower[:,i]=xo[:]
                llr2_lower[i]=fo[1]
                predict_r2_lower[:,:,i]=model(t_smooth,[nr2range_lower[1,i],r2range_lower[i],nr2range_lower[2,i]])

                loop_data_dists=[Normal(mi,nr2range_lower[3,i]) for mi in model(t_smooth,[nr2range_lower[1,i],r2range_lower[i],nr2range_lower[2,i]])]; # normal distribution about mean
                predict_r2_realisations_lower_lq[:,:,i]= [quantile(data_dist, THalpha/2) for data_dist in loop_data_dists];
                predict_r2_realisations_lower_uq[:,:,i]= [quantile(data_dist, 1-THalpha/2) for data_dist in loop_data_dists];
                

                if fo > fmle
                    println("new MLE r2 lower")
                    global fmle = fo
                    global r1mle=nr2range_lower[1,i]
                    global r2mle=r2range_lower[i]
                    global r3mle=nr2range_lower[2,i]
                    global Sdmle=nr2range_lower[3,i]
                end
            end
            
            # combine the lower and upper
            r2range = [reverse(r2range_lower);r2range_upper]
            nr2range = [reverse(nr2range_lower); nr2range_upper ]
            llr2 = [reverse(llr2_lower); llr2_upper]     
            nllr2=llr2.-maximum(llr2);

            max_r2=zeros(1,length(t_smooth))
            min_r2=1000*ones(1,length(t_smooth))
            for i in 1:(nptss)
                if (llr2_lower[i].-maximum(llr2)) >= TH
                    for j in 1:length(t_smooth)
                        max_r2[1,j]=max(predict_r2_lower[1,j,i],max_r2[1,j])
                        min_r2[1,j]=min(predict_r2_lower[1,j,i],min_r2[1,j])
                    end
                end
            end
            for i in 1:(nptss)
                if (llr2_upper[i].-maximum(llr2)) >= TH
                    for j in 1:length(t_smooth)
                        max_r2[1,j]=max(predict_r2_upper[1,j,i],max_r2[1,j])
                        min_r2[1,j]=min(predict_r2_upper[1,j,i],min_r2[1,j])
                    end
                end
            end
            # combine the confidence sets for realisations
            max_realisations_r2=zeros(1,length(t_smooth))
            min_realisations_r2=1000*ones(1,length(t_smooth))
            for i in 1:(nptss)
                if (llr2_lower[i].-maximum(llr2)) >= TH_realisations
                    for j in 1:length(t_smooth)
                        max_realisations_r2[1,j]=max(predict_r2_realisations_lower_uq[1,j,i],max_realisations_r2[1,j])
                        min_realisations_r2[1,j]=min(predict_r2_realisations_lower_lq[1,j,i],min_realisations_r2[1,j])
                    end
                end
            end
            for i in 1:(nptss)
                if (llr2_upper[i].-maximum(llr2)) >= TH_realisations
                    for j in 1:length(t_smooth)
                        max_realisations_r2[1,j]=max(predict_r2_realisations_upper_uq[1,j,i],max_realisations_r2[1,j])
                        min_realisations_r2[1,j]=min(predict_r2_realisations_upper_lq[1,j,i],min_realisations_r2[1,j]) 
                    end
                end
            end
            
            #Profile r3
            r3min=lb[3]
            r3max=ub[3]
            r3range_lower=reverse(LinRange(r3min,r3mle,nptss))
            r3range_upper=LinRange(r3mle + (r3max-r3mle)/nptss,r3max,nptss)
            
            nr3range_lower=zeros(3,nptss) # NOT r3, so r1, r2 and Sd
            llr3_lower=zeros(nptss)
            nllr3_lower=zeros(nptss) # Not necessary?
            predict_r3_lower=zeros(1,length(t_smooth),nptss)
            predict_r3_realisations_lower_lq=zeros(1,length(t_smooth),nptss)
            predict_r3_realisations_lower_uq=zeros(1,length(t_smooth),nptss)
            
            nr3range_upper=zeros(3,nptss) # NOT r3, so r1, r2 and Sd
            llr3_upper=zeros(nptss)
            nllr3_upper=zeros(nptss) # Not necessary?
            predict_r3_upper=zeros(1,length(t_smooth),nptss)
            predict_r3_realisations_upper_lq=zeros(1,length(t_smooth),nptss)
            predict_r3_realisations_upper_uq=zeros(1,length(t_smooth),nptss)
            
            # start at mle and increase parameter (upper)
            for i in 1:nptss
                function fun3(aa)
                    return error(data,[aa[1],aa[2],r3range_upper[i],aa[3]])
                end
                lb1=[lb[1],lb[2],lb[4]];
                ub1=[ub[1],ub[2],ub[4]];
                local θG1=[r1mle, r2mle, Sdmle]
                local (xo,fo)=optimise(fun3,θG1,lb1,ub1)
                nr3range_upper[:,i]=xo[:]
                llr3_upper[i]=fo[1]
                predict_r3_upper[:,:,i]=model(t_smooth,[nr3range_upper[1,i],nr3range_upper[2,i],r3range_upper[i]])

                loop_data_dists=[Normal(mi,nr3range_upper[3,i]) for mi in model(t_smooth,[nr3range_upper[1,i],nr3range_upper[2,i],r3range_upper[i]])]; # normal distribution about mean
                predict_r3_realisations_upper_lq[:,:,i]= [quantile(data_dist, THalpha/2) for data_dist in loop_data_dists];
                predict_r3_realisations_upper_uq[:,:,i]= [quantile(data_dist, 1-THalpha/2) for data_dist in loop_data_dists];

                if fo > fmle
                    println("new MLE r3 upper")
                    global fmle = fo
                    global r1mle=nr3range_upper[1,i]
                    global r2mle=nr3range_upper[2,i]
                    global r3mle=r3range_upper[i]
                    global Sdmle=nr3range_upper[3,i]
                end
            end
            
            # start at mle and decrease parameter (lower)
            for i in 1:nptss
                function fun3a(aa)
                    return error(data,[aa[1],aa[2],r3range_lower[i],aa[3]])
                end
                lb1=[lb[1],lb[2],lb[4]];
                ub1=[ub[1],ub[2],ub[4]];
                local θG1=[r1mle, r2mle, Sdmle]
                local (xo,fo)=optimise(fun3a,θG1,lb1,ub1)
                nr3range_lower[:,i]=xo[:]
                llr3_lower[i]=fo[1]
                predict_r3_lower[:,:,i]=model(t_smooth,[nr3range_lower[1,i],nr3range_lower[2,i],r3range_lower[i]])

                loop_data_dists=[Normal(mi,nr3range_lower[3,i]) for mi in model(t_smooth,[nr3range_lower[1,i],nr3range_lower[2,i],r3range_lower[i]])]; # normal distribution about mean
                predict_r3_realisations_lower_lq[:,:,i]= [quantile(data_dist, THalpha/2) for data_dist in loop_data_dists];
                predict_r3_realisations_lower_uq[:,:,i]= [quantile(data_dist, 1-THalpha/2) for data_dist in loop_data_dists];
                

                if fo > fmle
                    println("new MLE r3 lower")
                    global fmle = fo
                    global r1mle=nr3range_lower[1,i]
                    global r2mle=nr3range_lower[2,i]
                    global r3mle=r3range_lower[i]
                    global Sdmle=nr3range_lower[3,i]
                end
            end
            
            # combine the lower and upper
            r3range = [reverse(r3range_lower);r3range_upper]
            nr3range = [reverse(nr3range_lower); nr3range_upper ]
            llr3 = [reverse(llr3_lower); llr3_upper]     
            nllr3=llr3.-maximum(llr3);

            max_r3=zeros(1,length(t_smooth))
            min_r3=1000*ones(1,length(t_smooth))
            for i in 1:(nptss)
                if (llr3_lower[i].-maximum(llr3)) >= TH
                    for j in 1:length(t_smooth)
                        max_r3[1,j]=max(predict_r3_lower[1,j,i],max_r3[1,j])
                        min_r3[1,j]=min(predict_r3_lower[1,j,i],min_r3[1,j])
                    end
                end
            end
            for i in 1:(nptss)
                if (llr3_upper[i].-maximum(llr3)) >= TH
                    for j in 1:length(t_smooth)
                        max_r3[1,j]=max(predict_r3_upper[1,j,i],max_r3[1,j])
                        min_r3[1,j]=min(predict_r3_upper[1,j,i],min_r3[1,j])
                    end
                end
            end
            # combine the confidence sets for realisations
            max_realisations_r3=zeros(1,length(t_smooth))
            min_realisations_r3=1000*ones(1,length(t_smooth))
            for i in 1:(nptss)
                if (llr3_lower[i].-maximum(llr3)) >= TH_realisations
                    for j in 1:length(t_smooth)
                        max_realisations_r3[1,j]=max(predict_r3_realisations_lower_uq[1,j,i],max_realisations_r3[1,j])
                        min_realisations_r3[1,j]=min(predict_r3_realisations_lower_lq[1,j,i],min_realisations_r3[1,j])
                    end
                end
            end
            for i in 1:(nptss)
                if (llr3_upper[i].-maximum(llr3)) >= TH_realisations
                    for j in 1:length(t_smooth)
                        max_realisations_r3[1,j]=max(predict_r3_realisations_upper_uq[1,j,i],max_realisations_r3[1,j])
                        min_realisations_r3[1,j]=min(predict_r3_realisations_upper_lq[1,j,i],min_realisations_r3[1,j]) 
                    end
                end
            end
            

            #Profile Sd
            Sdmin=lb[4]
            Sdmax=ub[4]
            Sdrange_lower=reverse(LinRange(Sdmin,Sdmle,nptss))
            Sdrange_upper=LinRange(Sdmle + (Sdmax-Sdmle)/nptss,Sdmax,nptss)
            
            nSdrange_lower=zeros(3,nptss) # NOT Sd, so r1, r2, and r3
            llSd_lower=zeros(nptss)
            nllSd_lower=zeros(nptss) # Not necessary?
            predict_Sd_lower=zeros(1,length(t_smooth),nptss)
            predict_Sd_realisations_lower_lq=zeros(1,length(t_smooth),nptss)
            predict_Sd_realisations_lower_uq=zeros(1,length(t_smooth),nptss)

            nSdrange_upper=zeros(3,nptss) # NOT Sd, so r1, r2, and r3
            llSd_upper=zeros(nptss)
            nllSd_upper=zeros(nptss) # Not necessary?
            predict_Sd_upper=zeros(1,length(t_smooth),nptss)
            predict_Sd_realisations_upper_lq=zeros(1,length(t_smooth),nptss)
            predict_Sd_realisations_upper_uq=zeros(1,length(t_smooth),nptss)
            
            # start at mle and increase parameter (upper)
            for i in 1:nptss
                function fun4(aa)
                    return error(data,[aa[1],aa[2],aa[3],Sdrange_upper[i]])
                end
                lb1=[lb[1],lb[2],lb[3]];
                ub1=[ub[1],ub[2],ub[3]];
                local θG1=[r1mle,r2mle,r3mle]    
                local (xo,fo)=optimise(fun4,θG1,lb1,ub1)
                nSdrange_upper[:,i]=xo[:]
                llSd_upper[i]=fo[1]
                predict_Sd_upper[:,:,i]=model(t_smooth,[nSdrange_upper[1,i],nSdrange_upper[2,i],nSdrange_upper[3,i]])

                loop_data_dists=[Normal(mi,Sdrange_upper[i]) for mi in model(t_smooth,[nSdrange_upper[1,i],nSdrange_upper[2,i],nSdrange_upper[3,i]])]; # normal distribution about mean
                predict_Sd_realisations_upper_lq[:,:,i]= [quantile(data_dist, THalpha/2) for data_dist in loop_data_dists];
                predict_Sd_realisations_upper_uq[:,:,i]= [quantile(data_dist, 1-THalpha/2) for data_dist in loop_data_dists];
                
                if fo > fmle
                    println("new MLE Sd upper")
                    global fmle = fo
                    global r1mle=nSdrange_upper[1,i]
                    global r2mle=nSdrange_upper[2,i]
                    global r3mle=nSdrange_upper[3,i]
                    global Sdmle=Sdrange_upper[i]
                end
            end
            
            # start at mle and decrease parameter (lower)
            for i in 1:nptss
                function fun4a(aa)
                    return error(data,[aa[1],aa[2],aa[3],Sdrange_lower[i]])
                end
                lb1=[lb[1],lb[2],lb[3]];
                ub1=[ub[1],ub[2],ub[3]];
                local θG1=[r1mle,r2mle,r3mle]      
                local (xo,fo)=optimise(fun4a,θG1,lb1,ub1)
                nSdrange_lower[:,i]=xo[:]
                llSd_lower[i]=fo[1]
                predict_Sd_lower[:,:,i]=model(t_smooth,[nSdrange_upper[1,i],nSdrange_upper[2,i],nSdrange_upper[3,i]])

                loop_data_dists=[Normal(mi,Sdrange_lower[i]) for mi in model(t_smooth,[nSdrange_lower[1,i],nSdrange_lower[2,i],nSdrange_upper[3,i]])]; # normal distribution about mean
                predict_Sd_realisations_lower_lq[:,:,i]= [quantile(data_dist, THalpha/2) for data_dist in loop_data_dists];
                predict_Sd_realisations_lower_uq[:,:,i]= [quantile(data_dist, 1-THalpha/2) for data_dist in loop_data_dists];

                if fo > fmle
                    println("new MLE Sd lower")
                    global fmle = fo
                    global r1mle=nSdrange_lower[1,i]
                    global r2mle=nSdrange_lower[2,i]
                    global r3mle=nSdrange_lower[3,1]
                    global Sdmle=Sdrange_lower[i]
                end
            end
            
            # combine the lower and upper
            Sdrange = [reverse(Sdrange_lower);Sdrange_upper]
            nSdrange = [reverse(nSdrange_lower); nSdrange_upper ]
            llSd = [reverse(llSd_lower); llSd_upper] 
            nllSd=llSd.-maximum(llSd)

            max_Sd=zeros(1,length(t_smooth))
            min_Sd=1000*ones(1,length(t_smooth))
            for i in 1:(nptss)
                if (llSd_lower[i].-maximum(llSd)) >= TH
                    for j in 1:length(t_smooth)
                        max_Sd[1,j]=max(predict_Sd_lower[1,j,i],max_Sd[1,j])
                        min_Sd[1,j]=min(predict_Sd_lower[1,j,i],min_Sd[1,j])
                    end
                end
            end
            for i in 1:(nptss)
                if (llSd_upper[i].-maximum(llSd)) >= TH
                    for j in 1:length(t_smooth)
                        max_Sd[1,j]=max(predict_Sd_upper[1,j,i],max_Sd[1,j])
                        min_Sd[1,j]=min(predict_Sd_upper[1,j,i],min_Sd[1,j])
                    end
                end
            end

            # combine the confidence sets for realisations
            max_realisations_Sd=zeros(1,length(t_smooth))
            min_realisations_Sd=1000*ones(1,length(t_smooth))
            for i in 1:(nptss)
                if (llSd_lower[i].-maximum(llSd)) >= TH_realisations
                    for j in 1:length(t_smooth)
                        max_realisations_Sd[1,j]=max(predict_Sd_realisations_lower_uq[1,j,i],max_realisations_Sd[1,j])
                        min_realisations_Sd[1,j]=min(predict_Sd_realisations_lower_lq[1,j,i],min_realisations_Sd[1,j])
                    end
                end
            end
            for i in 1:(nptss)
                if (llSd_upper[i].-maximum(llSd)) >= TH_realisations
                    for j in 1:length(t_smooth)
                        max_realisations_Sd[1,j]=max(predict_Sd_realisations_upper_uq[1,j,i],max_realisations_Sd[1,j])
                        min_realisations_Sd[1,j]=min(predict_Sd_realisations_upper_lq[1,j,i],min_realisations_Sd[1,j]) 
                    end
                end
            end


            ##############################################################################################################
            # Plot Profile Likelihoods
        
            combined_plot_linewidth = 2;

            # interpolate for smoother profile likelihoods
            interp_nptss= 1001;
            
            # r1
            interp_points_r1range =  LinRange(r1min,r1max,interp_nptss)
            interp_r1 = LinearInterpolation(r1range,nllr1)
            interp_nllr1 = interp_r1(interp_points_r1range)
            
            profile1=plot(interp_points_r1range,interp_nllr1,xlim=(r1min,r1max),ylim=(-4,0.1),yticks=[-3,-2,-1,0],xlab=L"λ",ylab=L"\hat{\ell}_{p}",legend=false,lw=5,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=:deepskyblue3)
            profile1=hline!([-1.92],lw=2,linecolor=:black,linestyle=:dot)
            profile1=vline!([r1mle],lw=3,linecolor=:red)
            profile1=vline!([r1],lw=3,linecolor=:rosybrown,linestyle=:dash)

            # r2
            interp_points_r2range =  LinRange(r2min,r2max,interp_nptss)
            interp_r2 = LinearInterpolation(r2range,nllr2)
            interp_nllr2 = interp_r2(interp_points_r2range)
            
            profile2=plot(interp_points_r2range,interp_nllr2,xlim=(r2min,r2max),ylim=(-4,0.1),yticks=[-3,-2,-1,0],xlab=L"ζ",ylab=L"\hat{\ell}_{p}",legend=false,lw=5,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=:deepskyblue3)
            profile2=hline!([-1.92],lw=2,linecolor=:black,linestyle=:dot)
            profile2=vline!([r2mle],lw=3,linecolor=:red)
            profile2=vline!([r2],lw=3,linecolor=:rosybrown,linestyle=:dash)

            # r3
            interp_points_r3range =  LinRange(r3min,r3max,interp_nptss)
            interp_r3 = LinearInterpolation(r3range,nllr3)
            interp_nllr3 = interp_r3(interp_points_r3range)
            
            profile3=plot(interp_points_r3range,interp_nllr3,xlim=(r3min,r3max),ylim=(-4,0.1),yticks=[-3,-2,-1,0],xlab=L"R_{d}",ylab=L"\hat{\ell}_{p}",legend=false,lw=5,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=:deepskyblue3)
            profile3=hline!([-1.92],lw=2,linecolor=:black,linestyle=:dot)
            profile3=vline!([r3mle],lw=3,linecolor=:red)
            profile3=vline!([r3],lw=3,linecolor=:rosybrown,linestyle=:dash)
            
            # Sd
            interp_points_Sdrange =  LinRange(Sdmin,Sdmax,interp_nptss)
            interp_Sd = LinearInterpolation(Sdrange,nllSd)
            interp_nllSd = interp_Sd(interp_points_Sdrange)
            
            profile4=plot(interp_points_Sdrange,interp_nllSd,xlim=(Sdmin,Sdmax),ylim=(-4,0.1),yticks=[-3,-2,-1,0],xlab=L"\sigma_{N}",ylab=L"\hat{\ell}_{p}",legend=false,lw=5,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=:deepskyblue3)
            profile4=hline!([-1.92],lw=2,linecolor=:black,linestyle=:dot)
            profile4=vline!([Sdmle],lw=3,linecolor=:red)
            profile4=vline!([Sd],lw=3,linecolor=:rosybrown,linestyle=:dash)

            # Recompute best fit
            data0_smooth_MLE_recomputed = model(t_smooth,[r1mle,r2mle,r3mle,C0]); # generate data using known parameter values for plotting
            if size(data0_smooth_MLE_recomputed)[1] > num_x
                data0_smooth_MLE_recomputed = data0_smooth_MLE_recomputed'
            end

            # plot model simulated at MLE and data
            p1_updated=plot(t_smooth,data0_smooth_MLE_recomputed[1,:],xlab=L"t",ylab=L"R(t)",lw=3,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=colour2) # solid - R(t)
            p1_updated=scatter!(t,data[1,:],xlims=(xmin,xmax),ylims=(ymin,ymax),legend=false,markersize = 8,markercolor=colour2,markerstrokewidth=0) # scatter R(t)

        
            # save figures        
            display(p1_updated)
            savefig(p1_updated,filepath_save[1] * "Figp1_updated" * ".pdf")
            savefig(p1_updated,filepath_save[1] * "Figp1_updated" * ".png")

            display(profile1)
            savefig(profile1,filepath_save[1] * "Fig_profile1"   * ".pdf")
            savefig(profile1,filepath_save[1] * "Fig_profile1"   * ".png")
            
            display(profile2)
            savefig(profile2,filepath_save[1] * "Fig_profile2"   * ".pdf")
            savefig(profile2,filepath_save[1] * "Fig_profile2"   * ".png")
            
            display(profile3)
            savefig(profile3,filepath_save[1] * "Fig_profile3"   * ".pdf")
            savefig(profile3,filepath_save[1] * "Fig_profile3"   * ".png")

            display(profile4)
            savefig(profile4,filepath_save[1] * "Fig_profile4"   * ".pdf")
            savefig(profile4,filepath_save[1] * "Fig_profile4"   * ".png")


            # plot residuals at MLE
            t_residuals = t;
            data0_MLE_recomputed_for_residuals = model(t,[r1mle,r2mle,r3mle,C0]); # generate data using known parameter values for plotting
            if size(data0_MLE_recomputed_for_residuals)[1] > num_x
                data0_MLE_recomputed_for_residuals = data0_MLE_recomputed_for_residuals'
            end
            data_residuals = data-data0_MLE_recomputed_for_residuals
            extrme = maximum(abs,data_residuals)
            p_residuals=scatter(t,data_residuals[:],xlims=(xmin,xmax),ylims=(-extrme,extrme),legend=false,markersize =15,markercolor=colour2,xlab=L"t",ylab=L"\hat{e}_{i}",lw=3,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,markerstrokewidth=0) # scatter R(t)      
            p_residuals=hline!([0],lw=2,linecolor=:black,linestyle=:dot)

            display(p_residuals)
            savefig(p_residuals,filepath_save[1] * "Figp_residuals"   * ".pdf")
            savefig(p_residuals,filepath_save[1] * "Figp_residuals"   * ".png")

            # plot qq-plot
            p_residuals_qqplot = plot(qqplot(Normal,data_residuals[:],lw=2,linecolor=:black,linestyle=:dot,markersize = 8,markercolor=:black),legend=false,xlab=L"\mathrm{Normal}",ylab=L"\mathrm{Residuals}",lw=3,titlefont=fnt, guidefont=fnt, tickfont=fnt,ylims=(-15,15),yticks=[-10,0,10])
            display(p_residuals_qqplot)
            savefig(p_residuals_qqplot,filepath_save[1] * "Figp_residuals_qqplot"   * ".pdf")
            savefig(p_residuals_qqplot,filepath_save[1] * "Figp_residuals_qqplot"   * ".png")

            
            #######################################################################################
            # compute the bounds of confidence interval
            
            function fun_interpCI(mle,interp_points_range,interp_nll,TH)
                # find bounds of CI
                range_minus_mle = interp_points_range - mle*ones(length(interp_points_range),1)
                abs_range_minus_mle = broadcast(abs, range_minus_mle)
                findmin_mle = findmin(abs_range_minus_mle)
            
                # find closest value to CI threshold intercept
                value_minus_threshold = interp_nll - TH*ones(length(interp_nll),1)
                abs_value_minus_threshold = broadcast(abs, value_minus_threshold)
                lb_CI_tmp = findmin(abs_value_minus_threshold[1:findmin_mle[2][1]])
                ub_CI_tmp = findmin(abs_value_minus_threshold[findmin_mle[2][1]:length(abs_value_minus_threshold)])
                lb_CI = interp_points_range[lb_CI_tmp[2][1]]
                ub_CI = interp_points_range[findmin_mle[2][1]-1 + ub_CI_tmp[2][1]]
            
                return lb_CI,ub_CI
            end
            
            # r1
            (lb_CI_r1,ub_CI_r1) = fun_interpCI(r1mle,interp_points_r1range,interp_nllr1,TH)
            println(round(lb_CI_r1; digits = 4))
            println(round(ub_CI_r1; digits = 4))
            
            # r2
            (lb_CI_r2,ub_CI_r2) = fun_interpCI(r2mle,interp_points_r2range,interp_nllr2,TH)
            println(round(lb_CI_r2; digits = 4))
            println(round(ub_CI_r2; digits = 4))

            # r3
            (lb_CI_r3,ub_CI_r3) = fun_interpCI(r3mle,interp_points_r3range,interp_nllr3,TH)
            println(round(lb_CI_r3; digits = 4))
            println(round(ub_CI_r3; digits = 4))
            
            # Sd
            (lb_CI_Sd,ub_CI_Sd) = fun_interpCI(Sdmle,interp_points_Sdrange,interp_nllSd,TH)
            println(round(lb_CI_Sd; digits = 3))
            println(round(ub_CI_Sd; digits = 3))
            
            # Export MLE and bounds to csv (one file for all data) -- -MLE ONLY 
            if @isdefined(df_MLEBoundsAll) == 0
                println("not defined")
                global df_MLEBoundsAll = DataFrame( r1mle=r1mle, lb_CI_r1=lb_CI_r1, ub_CI_r1=ub_CI_r1, r2mle=r2mle, lb_CI_r2=lb_CI_r2, ub_CI_r2=ub_CI_r2, r3mle=r3mle, lb_CI_r3=lb_CI_r3, ub_CI_r3=ub_CI_r3, Sdmle=Sdmle, lb_CI_Sd=lb_CI_Sd, ub_CI_Sd=ub_CI_Sd)
            else 
                println("DEFINED")
                global df_MLEBoundsAll_thisrow = DataFrame( r1mle=r1mle, lb_CI_r1=lb_CI_r1, ub_CI_r1=ub_CI_r1, r2mle=r2mle, lb_CI_r2=lb_CI_r2, ub_CI_r2=ub_CI_r2, r3mle=r3mle, lb_CI_r3=lb_CI_r3, ub_CI_r3=ub_CI_r3, Sdmle=Sdmle, lb_CI_Sd=lb_CI_Sd, ub_CI_Sd=ub_CI_Sd)
                append!(df_MLEBoundsAll,df_MLEBoundsAll_thisrow)
            end
            
            CSV.write(filepath_save[1] * "MLEBoundsALL.csv", df_MLEBoundsAll)

    #########################################################################################
            # Confidence sets for model solutions
        
            ymle_smooth =  model(t_smooth,[r1mle,r2mle,r3mle,C0]); # MLE
            if size(ymle_smooth)[1] > num_x
                ymle_smooth = ymle_smooth'
            end

            # r1 - plot model simulated at MLE, scatter data, and prediction interval.
            confmodel1=plot(t_smooth,ymle_smooth[1,:],xlab=L"t",ylab=L"R(t)",lw=3,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=colour2) # solid - R(t)
            confmodel1=scatter!(t,data[1,:],xlims=(xmin,xmax),ylims=(ymin,ymax),legend=false,markersize = 8,markercolor=colour2,msw=0) # scatter R(t)
            confmodel1=plot!(t_smooth,ymle_smooth[1,:],w=0,c=colour2,ribbon=(ymle_smooth[1,:].-min_r1[1,:],max_r1[1,:].-ymle_smooth[1,:]),fillalpha=.2)
            
            display(confmodel1)
            savefig(confmodel1,filepath_save[1] * "Fig_confmodel1_lambda"   * ".pdf")
            savefig(confmodel1,filepath_save[1] * "Fig_confmodel1_lambda"   * ".png")

            # r2 - plot model simulated at MLE, scatter data, and prediction interval.
            confmodel2=plot(t_smooth,ymle_smooth[1,:],xlab=L"t",ylab=L"R(t)",lw=3,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=colour2) # solid - R(t)
            confmodel2=scatter!(t,data[1,:],xlims=(xmin,xmax),ylims=(ymin,ymax),legend=false,markersize = 8,markercolor=colour2,msw=0) # scatter R(t)
            confmodel2=plot!(t_smooth,ymle_smooth[1,:],w=0,c=colour2,ribbon=(ymle_smooth[1,:].-min_r2[1,:],max_r2[1,:].-ymle_smooth[1,:]),fillalpha=.2)
            
            display(confmodel2)
            savefig(confmodel2,filepath_save[1] * "Fig_confmodel2_zeta"   * ".pdf")
            savefig(confmodel2,filepath_save[1] * "Fig_confmodel2_zeta"   * ".png")
            
            # r3 - plot model simulated at MLE, scatter data, and prediction interval.
            confmodel3=plot(t_smooth,ymle_smooth[1,:],xlab=L"t",ylab=L"R(t)",lw=3,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=colour2) # solid - R(t)
            confmodel3=scatter!(t,data[1,:],xlims=(xmin,xmax),ylims=(ymin,ymax),legend=false,markersize = 8,markercolor=colour2,msw=0) # scatter R(t)
            confmodel3=plot!(t_smooth,ymle_smooth[1,:],w=0,c=colour2,ribbon=(ymle_smooth[1,:].-min_r3[1,:],max_r3[1,:].-ymle_smooth[1,:]),fillalpha=.2)
            
            display(confmodel3)
            savefig(confmodel3,filepath_save[1] * "Fig_confmodel3_Rd"   * ".pdf")
            savefig(confmodel3,filepath_save[1] * "Fig_confmodel3_Rd"   * ".png")

            # sd

            confmodel4=plot(t_smooth,ymle_smooth[1,:],xlab=L"t",ylab=L"R(t)",lw=3,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=colour2) # solid - R(t)
            confmodel4=scatter!(t,data[1,:],xlims=(xmin,xmax),ylims=(ymin,ymax),legend=false,markersize = 8,markercolor=colour2,msw=0) # scatter R(t)
            confmodel4=plot!(t_smooth,ymle_smooth[1,:],w=0,c=colour2,ribbon=(ymle_smooth[1,:].-min_Sd[1,:],max_Sd[1,:].-ymle_smooth[1,:]),fillalpha=.2)

            display(confmodel4)
            savefig(confmodel4,filepath_save[1] * "Fig_confmodel4_Sd"   * ".pdf")
            savefig(confmodel4,filepath_save[1] * "Fig_confmodel4_Sd"   * ".png")

            # union
            max_overall=zeros(1,length(t_smooth))
            min_overall=1000*ones(1,length(t_smooth))
            for j in 1:length(t_smooth)
                max_overall[1,j]=max(max_r1[1,j],max_r2[1,j],max_r3[1,j])
                min_overall[1,j]=min(min_r1[1,j],min_r2[1,j],min_r3[1,j])
            end
            confmodelu=plot(t_smooth,ymle_smooth[1,:],xlab=L"t",ylab=L"R(t)",lw=3,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=colour2) # solid - R(t)
            confmodelu=scatter!(t,data[1,:],xlims=(xmin,xmax),ylims=(ymin,ymax),legend=false,markersize = 8,markercolor=colour2,msw=0) # scatter R(t)
            confmodelu=plot!(t_smooth,ymle_smooth[1,:],w=0,c=colour2,ribbon=(ymle_smooth[1,:].-min_overall[1,:],max_overall[1,:].-ymle_smooth[1,:]),fillalpha=.2)

            display(confmodelu)
            savefig(confmodelu,filepath_save[1] * "Fig_confmodelu"   * ".pdf")
            savefig(confmodelu,filepath_save[1] * "Fig_confmodelu"   * ".png")
        



    #########################################################################################
            # Plot difference of confidence set for model solutions and the model solution at MLE 

            ydiffmin=-4.5
            ydiffmax= 4.5

            # r1
            confmodeldiff1 =plot(t_smooth,0 .*ymle_smooth[1,:],w=0,c=colour2,ribbon=(ymle_smooth[1,:].-min_r1[1,:],max_r1[1,:].-ymle_smooth[1,:]),ylims=(ydiffmin,ydiffmax),legend=false,fillalpha=.2,xlab=L"t",ylab=L"C_{y,0.95}^{λ} - y(\hat{\theta})",lw=0,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=colour2)
            confmodeldiff1 =plot!(t_smooth,0 .*ymle_smooth[1,:],lw=1,c=:black,linestyle=:dash)
            
            display(confmodeldiff1)
            savefig(confmodeldiff1,filepath_save[1] * "Fig_confmodeldiff1_lambda"   * ".pdf")
            savefig(confmodeldiff1,filepath_save[1] * "Fig_confmodeldiff1_lambda"   * ".png")

            # r2
            confmodeldiff2 =plot(t_smooth,0 .*ymle_smooth[1,:],w=0,c=colour2,ribbon=(ymle_smooth[1,:].-min_r2[1,:],max_r2[1,:].-ymle_smooth[1,:]),ylims=(ydiffmin,ydiffmax),legend=false,fillalpha=.2,xlab=L"t",ylab=L"C_{y,0.95}^{ζ} - y(\hat{\theta})",lw=0,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=colour2)
            confmodeldiff2 =plot!(t_smooth,0 .*ymle_smooth[1,:],lw=1,c=:black,linestyle=:dash)
            
            display(confmodeldiff2)
            savefig(confmodeldiff2,filepath_save[1] * "Fig_confmodeldiff2_zeta"   * ".pdf")
            savefig(confmodeldiff2,filepath_save[1] * "Fig_confmodeldiff2_zeta"   * ".png")

            # r3
            confmodeldiff3 =plot(t_smooth,0 .*ymle_smooth[1,:],w=0,c=colour2,ribbon=(ymle_smooth[1,:].-min_r3[1,:],max_r3[1,:].-ymle_smooth[1,:]),ylims=(ydiffmin,ydiffmax),legend=false,fillalpha=.2,xlab=L"t",ylab=L"C_{y,0.95}^{R_{d}} - y(\hat{\theta})",lw=0,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=colour2)
            confmodeldiff3 =plot!(t_smooth,0 .*ymle_smooth[1,:],lw=1,c=:black,linestyle=:dash)
            
            display(confmodeldiff3)
            savefig(confmodeldiff3,filepath_save[1] * "Fig_confmodeldiff3_Rd"   * ".pdf")
            savefig(confmodeldiff3,filepath_save[1] * "Fig_confmodeldiff3_Rd"   * ".png")
            
            # sd
            confmodeldiff4=plot(t_smooth,0 .*ymle_smooth[1,:],w=0,c=colour2,ribbon=(ymle_smooth[1,:].-min_Sd[1,:],max_Sd[1,:].-ymle_smooth[1,:]),ylims=(ydiffmin,ydiffmax),legend=false,fillalpha=.2,xlab=L"t",ylab=L"C_{y,0.95}^{\sigma_{\mathrm{N}}} - y(\hat{\theta})",lw=0,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=colour2)
            confmodeldiff4=plot!(t_smooth,0 .*ymle_smooth[1,:],lw=1,c=:black,linestyle=:dash)
    
            display(confmodeldiff4)
            savefig(confmodeldiff4,filepath_save[1] * "Fig_confmodeldiff4_Sd"   * ".pdf")
            savefig(confmodeldiff4,filepath_save[1] * "Fig_confmodeldiff4_Sd"   * ".png")

            # union
            confmodeldiffu=plot(t_smooth,0 .*ymle_smooth[1,:],w=0,c=colour2,ribbon=(ymle_smooth[1,:].-min_overall[1,:],max_overall[1,:].-ymle_smooth[1,:]),ylims=(ydiffmin,ydiffmax),legend=false,fillalpha=.2,xlab=L"t",ylab=L"C_{y,0.95} - y(\hat{\theta})",lw=0,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=colour2)
            confmodeldiffu=plot!(t_smooth,0 .*ymle_smooth[1,:],lw=1,c=:black,linestyle=:dash)

            display(confmodeldiffu)
            savefig(confmodeldiffu,filepath_save[1] * "Fig_confmodeldiffu"   * ".pdf")
            savefig(confmodeldiffu,filepath_save[1] * "Fig_confmodeldiffu"   * ".png")

        


        
        
    #########################################################################################
            # Confidence sets for model realisations
                
            # r1
            confreal1=plot(t_smooth,ymle_smooth[1,:],xlab=L"t",ylab=L"R(t)",lw=3,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=colour2) # solid - R(t)
            confreal1=scatter!(t,data[1,:],xlims=(xmin,xmax),ylims=(ymin,ymax),legend=false,markersize = 8,markercolor=colour2,msw=0) # scatter R(t)
            confreal1=plot!(t_smooth,ymle_smooth[1,:],w=0,c=colour2,ribbon=(ymle_smooth[1,:].-min_realisations_r1[1,:],max_realisations_r1[1,:].-ymle_smooth[1,:]),fillalpha=.2)
        
            display(confreal1)
            savefig(confreal1,filepath_save[1] * "Fig_confreal1"   * ".pdf")
            savefig(confreal1,filepath_save[1] * "Fig_confreal1"   * ".png")

            # r2
            confreal2=plot(t_smooth,ymle_smooth[1,:],xlab=L"t",ylab=L"R(t)",lw=3,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=colour2) # solid - R(t)
            confreal2=scatter!(t,data[1,:],xlims=(xmin,xmax),ylims=(ymin,ymax),legend=false,markersize = 8,markercolor=colour2,msw=0) # scatter R(t)
            confreal2=plot!(t_smooth,ymle_smooth[1,:],w=0,c=colour2,ribbon=(ymle_smooth[1,:].-min_realisations_r2[1,:],max_realisations_r2[1,:].-ymle_smooth[1,:]),fillalpha=.2)
            
            display(confreal2)
            savefig(confreal2,filepath_save[1] * "Fig_confreal2"   * ".pdf")
            savefig(confreal2,filepath_save[1] * "Fig_confreal2"   * ".png")

            # r3
            confreal3=plot(t_smooth,ymle_smooth[1,:],xlab=L"t",ylab=L"R(t)",lw=3,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=colour2) # solid - R(t)
            confreal3=scatter!(t,data[1,:],xlims=(xmin,xmax),ylims=(ymin,ymax),legend=false,markersize = 8,markercolor=colour2,msw=0) # scatter R(t)
            confreal3=plot!(t_smooth,ymle_smooth[1,:],w=0,c=colour2,ribbon=(ymle_smooth[1,:].-min_realisations_r3[1,:],max_realisations_r3[1,:].-ymle_smooth[1,:]),fillalpha=.2)
        
            display(confreal3)
            savefig(confreal3,filepath_save[1] * "Fig_confreal3"   * ".pdf")
            savefig(confreal3,filepath_save[1] * "Fig_confreal3"   * ".png")

            #Sd
            confreal4=plot(t_smooth,ymle_smooth[1,:],xlab=L"t",ylab=L"R(t)",lw=3,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=colour2) # solid - R(t)
            confreal4=scatter!(t,data[1,:],xlims=(xmin,xmax),ylims=(ymin,ymax),legend=false,markersize = 8,markercolor=colour2,msw=0) # scatter R(t)
            confreal4=plot!(t_smooth,ymle_smooth[1,:],w=0,c=colour2,ribbon=(ymle_smooth[1,:].-min_realisations_Sd[1,:],max_realisations_Sd[1,:].-ymle_smooth[1,:]),fillalpha=.2)
            
            display(confreal4)
            savefig(confreal4,filepath_save[1] * "Fig_confreal4"   * ".pdf")
            savefig(confreal4,filepath_save[1] * "Fig_confreal4"   * ".png")


            # overall - plot model simulated at MLE, scatter data, and prediction interval (union of confidence sets for realisations)
            max_realisations_overall=zeros(1,length(t_smooth))
            min_realisations_overall=1000*ones(1,length(t_smooth))
            for j in 1:length(t_smooth)
                max_realisations_overall[1,j]=max(max_realisations_r1[1,j],max_realisations_r2[1,j],max_realisations_r3[1,j],max_realisations_Sd[1,j])
                min_realisations_overall[1,j]=min(min_realisations_r1[1,j],min_realisations_r2[1,j],min_realisations_r3[1,j],min_realisations_Sd[1,j])
            end

            confrealu=plot(t_smooth,ymle_smooth[1,:],xlab=L"t",ylab=L"R(t)",lw=3,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=colour2) # solid - R(t)
            confrealu=scatter!(t,data[1,:],xlims=(xmin,xmax),ylims=(ymin,ymax),legend=false,markersize = 8,markercolor=colour2,msw=0) # scatter R(t)
            confrealu=plot!(t_smooth,ymle_smooth[1,:],w=0,c=colour2,ribbon=(ymle_smooth[1,:].-min_realisations_overall[1,:],max_realisations_overall[1,:].-ymle_smooth[1,:]),fillalpha=.2)
            
            display(confrealu)
            savefig(confrealu,filepath_save[1] * "Fig_confrealu"   * ".pdf")
            savefig(confrealu,filepath_save[1] * "Fig_confrealu"   * ".png")

        ##############################
            # PLot difference of confidence sets for data realisations and model solution at MLE

            ynoisediffmin=-15
            ynoisediffmax=15

            # r1
            confrealdiff1=plot(t_smooth,0 .* ymle_smooth[1,:],w=0,c=colour2,ribbon=(ymle_smooth[1,:].-min_realisations_r1[1,:],max_realisations_r1[1,:].-ymle_smooth[1,:]),ylims=(ynoisediffmin,ynoisediffmax),legend=false,fillalpha=.2,xlab=L"t",ylab=L"C_{z_{i},0.95}^{λ} - y(\hat{\theta})",lw=0,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=colour2)
            confrealdiff1=plot!(t_smooth,0 .*ymle_smooth[1,:],lw=1,c=:black,linestyle=:dash)

            display(confrealdiff1)
            savefig(confrealdiff1,filepath_save[1] * "Fig_confrealdiff1_lambda"   * ".pdf")
            savefig(confrealdiff1,filepath_save[1] * "Fig_confrealdiff1_lambda"   * ".png")

            # r2
            confrealdiff2=plot(t_smooth,0 .* ymle_smooth[1,:],w=0,c=colour2,ribbon=(ymle_smooth[1,:].-min_realisations_r2[1,:],max_realisations_r2[1,:].-ymle_smooth[1,:]),ylims=(ynoisediffmin,ynoisediffmax),legend=false,fillalpha=.2,xlab=L"t",ylab=L"C_{z_{i},0.95}^{ζ} - y(\hat{\theta})",lw=0,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=colour2)
            confrealdiff2=plot!(t_smooth,0 .*ymle_smooth[1,:],lw=1,c=:black,linestyle=:dash)

            display(confrealdiff2)
            savefig(confrealdiff2,filepath_save[1] * "Fig_confrealdiff2_zeta"   * ".pdf")
            savefig(confrealdiff2,filepath_save[1] * "Fig_confrealdiff2_zeta"   * ".png")

            # r3
            confrealdiff3=plot(t_smooth,0 .* ymle_smooth[1,:],w=0,c=colour2,ribbon=(ymle_smooth[1,:].-min_realisations_r3[1,:],max_realisations_r3[1,:].-ymle_smooth[1,:]),ylims=(ynoisediffmin,ynoisediffmax),legend=false,fillalpha=.2,xlab=L"t",ylab=L"C_{z_{i},0.95}^{R_{d}} - y(\hat{\theta})",lw=0,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=colour2)
            confrealdiff3=plot!(t_smooth,0 .*ymle_smooth[1,:],lw=1,c=:black,linestyle=:dash)

            display(confrealdiff3)
            savefig(confrealdiff3,filepath_save[1] * "Fig_confrealdiff3_Rd"   * ".pdf")
            savefig(confrealdiff3,filepath_save[1] * "Fig_confrealdiff3_Rd"   * ".png")

            #Sd
            confrealdiff4=plot(t_smooth,0 .* ymle_smooth[1,:],w=0,c=colour2,ribbon=(ymle_smooth[1,:].-min_realisations_Sd[1,:],max_realisations_Sd[1,:].-ymle_smooth[1,:]),ylims=(ynoisediffmin,ynoisediffmax),legend=false,fillalpha=.2,xlab=L"t",ylab=L"C_{z_{i},0.95}^{\sigma_{\mathrm{N}}} - y(\hat{\theta})",lw=0,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=colour2)
            confrealdiff4=plot!(t_smooth,0 .*ymle_smooth[1,:],lw=1,c=:black,linestyle=:dash)

            display(confrealdiff4)
            savefig(confrealdiff4,filepath_save[1] * "Fig_confrealdiff4_Sd"   * ".pdf")
            savefig(confrealdiff4,filepath_save[1] * "Fig_confrealdiff4_Sd"   * ".png")


            # overall - plot model simulated at MLE, scatter data, and prediction interval (union of confidence sets for realisations)

            confrealdiffu=plot(t_smooth,0 .* ymle_smooth[1,:],w=0,c=colour2,ribbon=(ymle_smooth[1,:].-min_realisations_overall[1,:],max_realisations_overall[1,:].-ymle_smooth[1,:]),ylims=(ynoisediffmin,ynoisediffmax),legend=false,fillalpha=.2,xlab=L"t",ylab=L"C_{z_{i},0.95} - y(\hat{\theta})",lw=0,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=colour2)
            confrealdiffu=plot!(t_smooth,0 .*ymle_smooth[1,:],lw=1,c=:black,linestyle=:dash)

            display(confrealdiffu)
            display(confrealdiffu) # Not necessary?
            savefig(confrealdiffu,filepath_save[1] * "Fig_confrealdiffu"   * ".pdf")
            savefig(confrealdiffu,filepath_save[1] * "Fig_confrealdiffu"   * ".png")

    print([r1mle, r2mle, r3mle, Sdmle]) 
