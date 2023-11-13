using Printf

# List all the .jl files in the folder
function list_jl_files(path)
    return filter(f -> endswith(f, ".jl"), readdir(path))
end

# Time the execution of a .jl file
function time_execution(file_path)
    t = @elapsed include(file_path)
    return t
end

# Run and time all .jl files in a folder
function run_and_time_all(path)
    files = list_jl_files(path)
    timings = Dict{String,Float64}()

    for file in files
        timings[file] = time_execution(joinpath(path, file))
    end

    return timings
end

# Display the timings in a table
function display_timings(timings)
    println("+", "-"^50, "+", "-"^10, "+")
    @printf("| %-48s | %-8s |\n", "File Name", "Time (s)")
    println("+", "-"^50, "+", "-"^10, "+")

    for (file, time) in timings
        @printf("| %-48s | %8.4f |\n", file, time)
    end

    println("+", "-"^50, "+", "-"^10, "+")
end

# Example usage:
folder_path = joinpath(@__DIR__)  # Replace with your folder's path
timings = run_and_time_all(folder_path)
display_timings(timings)
