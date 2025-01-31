# Create a struct that can store information in the following form:
#       store_variables([name1, name 2, ...], [true_value1, true_value2, ...], [lower_bound1, lower_bound2, ...], [upper_bound1, upper_bound2,...], [guess1, guess2, ...])
# Or if true values are unknown:
#       store_variables([name1, name 2, ...], [lower_bound1, lower_bound2, ...], [upper_bound1, upper_bound2,...], [guess1, guess2, ...])
# Or if true values and bounds are unnecessary:
#       store_variables([name1, name 2, ...], [guess1, guess2, ...])
# Or if true values, bounds and guesses are unnecessary:
#       store_variables([name1, name 2, ...])

struct store_variables # DO NOT EDIT
    name::Union{Vector{String}, String}                     # String for names
    true_value::Union{Vector{Float64}, Vector{Int64}, Float64, Int64, Nothing}      # True values if known
    lower_bound::Union{Vector{Float64}, Vector{Int64}, Float64, Int64, Nothing}     # Lower bounds to guesses
    upper_bound::Union{Vector{Float64}, Vector{Int64}, Float64, Int64, Nothing}     # Upper bounds to guesses
    guess::Union{Vector{Float64}, Vector{Int64}, Float64, Int64, Nothing}                    # Initial guess
    # Constructor function for all five fields
    function store_variables(name::Union{Vector{String}, String}, 
        true_value::Union{Vector{Float64}, Vector{Int64}, Float64, Int64, Nothing},
        lower_bound::Union{Vector{Float64}, Vector{Int64}, Float64, Int64, Nothing}, 
        upper_bound::Union{Vector{Float64}, Vector{Int64},Float64, Int64, Nothing},  
        guess::Union{Vector{Float64}, Vector{Int64}, Float64, Int64, Nothing})
        return new(name, true_value, lower_bound, upper_bound, guess)
    end
    # Constructor function for four fields
    function store_variables(name::Union{Vector{String}, String}, 
        lower_bound::Union{Vector{Float64}, Vector{Int64}, Float64, Int64, Nothing}, 
        upper_bound::Union{Vector{Float64}, Vector{Int64}, Float64, Int64, Nothing},  
        guess::Union{Vector{Float64}, Vector{Int64}, Float64, Int64, Nothing})
        return new(name, nothing, lower_bound, upper_bound, guess) # Assign nothing to true_value
    end
    # Constructor function for two fields
    function store_variables(name::Union{Vector{String}, String}, 
        guess::Union{Vector{Float64}, Vector{Int64}, Float64, Int64, Nothing})
        return new(name, nothing, nothing, nothing, guess)  # Assign nothing to true_value and bounds
    end
    # Constructor function for one field 
    function store_variables(name::Union{Vector{String}, String})
        return new(name, nothing, nothing, nothing, nothing)  # Assign nothing to true_value, bounds and initiial_values
    end
end

function defineStruct()
    # Define function to create the struct
    return store_variables
end