#######################################################################################
## Functions involved in extracting and storing the appropriate data for the analysis
#######################################################################################

function storeMultiVars(args ...)
    # Stores the given arrays in the one object: a matrix if time-point restrictions have been implemented,
    # or a dictionary if time-point restrictions have not been implemented
    
    # Dowling (2024)
    num_args = length(args) # Check the number of arguments
    if num_args == 2     
        array_var1, array_var2 = args
        if array_var1 == array_var2
            # If the rows in the matrix are identical, only return one
            all_arrays = copy(array_var1)
        elseif length(array_var1) == length(array_var2)
            # If time-point restrictions have been implemented, store in matrix
            all_arrays = [array_var1; array_var2] 
        else
            # Otherwise, store in a dictionary
            all_arrays = Dict{Int, Vector{Float64}}()
        end
    elseif num_args == 3
        array_var1, array_var2, array_var3 = args
        if array_var1 == array_var2
            # If the rows in the matrix are identical, only return one
            all_arrays = copy(array_var1)
        elseif length(array_var1) == length(array_var2) && length(array_var3) == length(array_var2)
            # If time-point restrictions have been implemented, store in matrix
            all_arrays = [array_var1; array_var2; array_var3]
        else
            # Otherwise, store in a dictionary
            all_arrays = Dict{Int, Vector{Float64}}()
        end
    else
        # Issue warning
        println("storeMultiVars: should have 2-3 inputs corresponding to distinct state variables.\n")
    end

    # If the combination of arrays is to be stored in a dictionary, fill in the dictionary
    if typeof(all_arrays) == Dict{Int64, Vector{Float64}}
        for j in 1:num_args
            if j == 1
                all_arrays[j] = array_var1[:]
            elseif j == 2
                all_arrays[j] = array_var2[:]
            elseif j == 3
                all_arrays[j] = array_var3[:]
            end
        end
    end

    return all_arrays
end

function loadExtractData(args ...)
    # Loads and extracts data at the desired user-defined time-points.
    
    # Dowling (2024)
    num_args = length(args) # Check the number of arguments
    if num_args == 11      
        file_name_all,file_name_specific,restrict_measures_per_tp,initial_cell_number,day_growth_begins,desired_time_points,type_of_data,core_name,variables_included,seed_num,restricted_num = args
    elseif num_args == 10 && typeof(args[2]) == String   # Define the second argument as the second file name
        file_name_all,file_name_specific,initial_cell_number,day_growth_begins,desired_time_points,type_of_data,core_name,variables_included,seed_num,restricted_num = args
    elseif num_args == 10 && typeof(args[2]) != String   # Define the second argument as the restricted time-points check
        file_name_all,restrict_measures_per_tp,initial_cell_number,day_growth_begins,desired_time_points,type_of_data,core_name,variables_included,seed_num,restricted_num = args
    else
        # Issue warning
        println("loadExtractData: should have either of the following sets of 10-11 inputs:\n")
        println("\tfile_name_all,file_name_specific,restrict_measures_per_tp,initial_cell_number,day_growth_begins,desired_time_points,type_of_data\n")
        println("\tfile_name_all,file_name_specific,initial_cell_number,day_growth_begins,desired_time_points,type_of_data\n")
        println("\tfile_nam_all,restrict_measures_per_tp,initial_cell_number,day_growth_begins,desired_time_points,type_of_data\n")
    end

    # Load all* of the data (*depending on the number of arguments, this may be the raw data, or tidy data without time-points being specifically extracted)
    @suppress begin                                                             # Suppress warnings about empty columns
        dataFile = DataFrame(skipmissing(CSV.File(file_name_all)))                  # Load data
        @warn("Warnings arising from loading the CSV file are ignored.")
    end    
    dataFile = dataFile[!, map(x->!all(ismissing, x), eachcol(dataFile))]       # Remove all columns with entirely missing entries
    
    # Depending on the data-type loaded, proceed to filter and tidy the measurements
    if type_of_data == "raw"
        # Extract data from the beginning of tumour growth for the desired cell seeding number
        dataInitSeedSpecified = sort(                                           # Extract all rows with data from experiments seeded 
            filter(row -> row.InitialCondition ∈ [initial_cell_number]                  # with 5000 cells from day 4, and sort data by the day 
            && row.Day >= day_growth_begins,                                    # of the measurement, then the size of the measurement
            dataFile),[:Day,:Radius]) 
        data_all = dataInitSeedSpecified.Radius'                                # Extract all of the radius data (one row per state variable)
        times_all = dataInitSeedSpecified.Day'                                  # Extract all of the times at which radius data is measured (one row)
        times_all = times_all .- day_growth_begins                              # Subtract 4 from all times to remove days of tumour formation
    else
        # Assume data has already been extracted from the beginning of tumour growth for the desired cell seeding number
        data_all = dataFile.DATA_ALL'                                           # Ensure one row per state variable
        times_all = dataFile.TIMES_ALL'                                         # Ensure one row 
    end

    # Depending on the data-type loaded, proceed to load or extract the measurements at the desired time-points
    if type_of_data == "extracted_tps"
        # Load the second file containing the previously extracted measurements at the desired time-points
        dataFileSpecific = DataFrame(skipmissing(CSV.File(file_name_specific)))
        data_specific = dataFileSpecific.DATA_SPECIFIC'                         # Ensure one row per state variable
        times_specific = dataFileSpecific.TIMES_SPECIFIC'                       # Ensure one row 
    else
        # Extract data at the given specific time-points defined in desired_time_points
        data_specific = [data_all[i] for i in eachindex(times_all) if times_all[i] ∈ desired_time_points]'  # Ensure one row per state variable
        times_specific = [x for x in times_all if x ∈ desired_time_points]'                                 # Ensure one row 
    end

    # If specified, reduce the number of measurements at each time-point to the minimum number of measurements made at any time-point
    if @isdefined(restrict_measures_per_tp) # Check if the restriction time-points check is defined (or if the data is desired to be used as is)
        if restrict_measures_per_tp == 1    # If the check is defined and restriction of measurement numbers at each time-point is desired
            num_of_timepoints_per_time = [count(==(i),times_all) for i in desired_time_points]
            restricted_num_potential = minimum(num_of_timepoints_per_time)
            if restricted_num==0    # If a restricted number is not given, assign it as the smallest number
                restricted_num = restricted_num_potential
            else                    # If a restricted number is given, assign it as the smaller of the two
                restricted_num = min(restricted_num_potential,restricted_num)
            end
            data_specific_restricted = []
            for time_point in desired_time_points
                data_specific_restricted_at_tp = sample(data_all[findall(x->x==time_point,times_all)],restricted_num,replace=false)
                size(data_specific)
                data_specific_restricted = [data_specific_restricted; data_specific_restricted_at_tp]
            end
            data_specific = data_specific_restricted'                # Ensure one row per state variable
            times_specific_restricted = [x .* ones(1,restricted_num) for x in desired_time_points]
            times_specific = reduce(hcat,times_specific_restricted)  # Ensure one row 
        end
    end

    # Depending on the data-type loaded, save the resulting vectors in the current directory
    # Save the tidy and extracted data as CSVs - comment out if it is already saved
    if type_of_data == "raw"
        # Save the tidy data and times at all valid time-points for the desired cell seeding number
        df_data_all = DataFrame(DATA_ALL = vec(data_all), TIMES_ALL = vec(times_all))
        CSV.write(core_name * variables_included * "_tumourGrowthData.csv", df_data_all) 
    end

    # Save the specific data and times at the desired valid time-points for the desired cell seeding number
    df_data_specific = DataFrame(DATA_SPECIFIC = vec(data_specific), TIMES_SPECIFIC = vec(times_specific))
    if @isdefined(restrict_measures_per_tp) # Save the tidy data and times at specific time-points, restricted by the number at each time-point
        if restrict_measures_per_tp == 1
            CSV.write(core_name * variables_included * "_tumourGrowthData_specificTimePoints_restrictNumMsrmntsPerTp.csv", df_data_specific)
        else                                # Save the tidy data and times at specific time-points
            CSV.write(core_name * variables_included * "_tumourGrowthData_specificTimePoints_rawNumMsrmntsPerTp.csv", df_data_specific)
        end
    end

    return data_all, times_all, data_specific, times_specific, restricted_num
end