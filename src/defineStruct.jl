# Create a struct that can store information in the following form:
#       store_variables([name1, name 2, ...], [true_value1, true_value2, ...], [lower_bound1, lower_bound2, ...], [upper_bound1, upper_bound2,...], [initial_value1, initial_value2, ...])
# Or if true values are unknown:
#       store_variables([name1, name 2, ...], [lower_bound1, lower_bound2, ...], [upper_bound1, upper_bound2,...], [guess_value1, initial_value2, ...])
# Or if true values and bounds are unnecessary:
#       store_variables([name1, name 2, ...], [guess_value1, initial_value2, ...])

struct store_variables # DO NOT EDIT
    names::Union{Vector{String}, String}                    # String for names
    true_values::Union{Vector{Float64}, Float64, Int64, Nothing}   # True values if known
    lower_bounds::Union{Vector{Float64}, Float64, Int64, Nothing}  # Lower bounds
    upper_bounds::Union{Vector{Float64}, Float64, Int64, Nothing}  # Upper bounds
    initial_values::Union{Vector{Float64}, Float64, Int64}         # Initial guesses or initial conditions
    # Constructor function for all five fields
    function store_variables(names::Union{Vector{String}, String}, 
        true_values::Union{Vector{Float64}, Float64, Int64, Nothing},
        lower_bounds::Union{Vector{Float64}, Float64, Int64, Nothing}, 
        upper_bounds::Union{Vector{Float64}, Float64, Int64, Nothing},  
        initial_values::Union{Vector{Float64}, Float64, Int64})
        return new(names, true_values, lower_bounds, upper_bounds, initial_values)
    end
    # Constructor function for four fields
    function store_variables(names::Union{Vector{String}, String}, 
        lower_bounds::Union{Vector{Float64}, Float64, Int64, Nothing}, 
        upper_bounds::Union{Vector{Float64}, Float64, Int64, Nothing},  
        initial_values::Union{Vector{Float64}, Float64, Int64})
        return new(names, nothing, lower_bounds, upper_bounds, initial_values) # Assign nothing to true_values
    end
    # Constructor function for two fields (default true_values to nothing)
    function store_variables(names::Union{Vector{String}, String}, 
        initial_values::Union{Vector{Float64}, Float64, Int64})
        return new(names, nothing, nothing, nothing, initial_values)  # Assign nothing to true_values and bounds
    end
end

function defineStruct()
    # Define function to create the struct
    return store_variables
end