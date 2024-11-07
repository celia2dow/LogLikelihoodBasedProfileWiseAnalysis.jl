function funcInterpCI(mle,interp_points_range,interp_nll,TH)
    # Function for interpolating a confidence interval

    # Dowling (2024)
    # Adapted from Murphy et al. (2023) 
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