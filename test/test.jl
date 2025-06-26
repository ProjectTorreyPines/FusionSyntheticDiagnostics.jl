
using Pkg
using Printf

pkg_name = "FusionSyntheticDiagnostics"

valid_arguments = [
    ("ifo", "Test interferometer"),
    ("langmuir", "Test langmuir probes"),
    ("magnetics", "Test magnetic diagnostics"),
    ("gas", "Test gas injection actuators"),
    ("h", "Show this help message and exit"),
    ("help", "Show this help message and exit"),
]

function get_options_string(
    valid_arguments::Vector{Tuple{String, String}};
    indent::Int=4,
)
    keys_arr = [k for (k, v) ∈ valid_arguments]
    values_arr = [v for (k, v) ∈ valid_arguments]

    # Determine the maximum width of each column (keys and values)
    max_width_keys = maximum(length.(string.(keys_arr)))
    max_width_values = maximum(length.(string.(values_arr)))

    # Create the format string with indentation
    format_str = repeat(" ", indent) * "%-$(max_width_keys)s %-$(max_width_values)s\n"

    # Create a Format object
    fmt = Printf.Format(format_str)

    # Build the combined string
    combined_string = ""
    for (key, value) ∈ valid_arguments
        combined_string *= Printf.format(fmt, string(key), string(value))
    end

    return combined_string
end

function check_args(
    args::Vector{String},
    valid_arguments::Vector{Tuple{String, String}},
)
    valid_keys = Set(t[1] for t ∈ valid_arguments)

    any_valid = false
    for arg ∈ args
        if arg in valid_keys
            any_valid = true
        else
            println("Error: Invalid argument found: \"", arg, "\"")
            println()
            return false
        end
    end

    if !any_valid && !isempty(args)
        println("Error: No valid arguments found.")
        return false
    end
    return true
end

usage_string = "Usage (from inside $(pkg_name).jl): \njulia --project test/test.jl "
usage_string *= "[" * join([k for (k, v) ∈ valid_arguments], "] [") * "]\n\n"
usage_string *= "Run tests. Default is all tests.\n\nOptional arguments:\n"
usage_string *= get_options_string(valid_arguments)

if "h" in ARGS || "help" in ARGS || !check_args(ARGS, valid_arguments)
    println(usage_string)
    exit()
end
try
    Pkg.test("$(pkg_name)"; test_args=ARGS)
catch e
    println("Error occured.")
    println(e)
    println(usage_string)
end
