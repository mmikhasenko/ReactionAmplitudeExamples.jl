using Pkg

# Shared configuration that matches CI
const PLUTO_VERSION = "0.3.2-0.3"
const CACHE_DIR = "test_pluto_state_cache"


# from https://github.com/JuliaPluto/PlutoSliderServer.jl/issues/73
"""
Given a selection of notebooks, check individually for each notebook if there are 
any cells that errored.
If it is the case print out, for each notebook, the content of the failing cells
and the output messages.
"""
function check_for_failed_notebooks(result::NamedTuple)
    failed_notebooks = Dict{String, Vector}()
    # notebook session is a `NotebookSession` from PlutoSliderServer.jl.
    for notebook_session in result.notebook_sessions
        # Check for every notebook that no cell errored.
        # State is a large JSON style `Dict` containing all the informations about the ran notebook.
        # You can find the definition in Pluto.jl/src/webserver/Dynamic.jl/notebook_to_js.
        state = notebook_session.run.original_state
        errored_cells = findall(cell -> cell["errored"], state["cell_results"])
        isempty(errored_cells) && continue
        failed_notebooks[notebook_session.path] = [
            (input = state["cell_inputs"][id]["code"], output = state["cell_results"][id]["output"]["body"][:msg]) for
            id in sort(errored_cells; by = id -> findfirst(==(id), state["cell_order"]))
        ]
    end
    if !isempty(failed_notebooks)
        io = IOBuffer()
        for (key, cells) in pairs(failed_notebooks)
            printstyled(IOContext(io, :color => true), "$key:\n"; bold = true, color = :green, underline = true)
            for (input, output) in cells
                printstyled(IOContext(io, :color => true), "• $input"; color = :blue)
                print(io, " => ")
                printstyled(IOContext(io, :color => true), "$output\n"; color = :red)
            end
            println(io)
        end
        error_msgs = String(take!(io))
        error(
            "The following Pluto notebook",
            length(failed_notebooks) > 1 ? "s" : "",
            " failed to run successfully: $(keys(failed_notebooks))\n\n",
            error_msgs,
        )
    end
end



# Create a test environment like in CI
tmp_env = mktempdir()
Pkg.activate(tmp_env)

# Use exact same package specs as CI
Pkg.add([
    Pkg.PackageSpec(name = "PlutoSliderServer", version = PLUTO_VERSION),
])

using PlutoSliderServer: PlutoSliderServer

# Create a test cache directory
mkpath(CACHE_DIR)

println("Testing notebook exports...")
# Use same export configuration as CI
PlutoSliderServer.export_directory(
    dirname(@__DIR__);  # Go up one level to reach project root
    Export_cache_dir = CACHE_DIR,
    Export_baked_notebookfile = false,
    Export_baked_state = false,
    on_ready = check_for_failed_notebooks,
)

println("\nTest completed! If no errors occurred, notebooks should work in CI.")
println("You can find the test cache in: ", CACHE_DIR)

# Cleanup
rm(CACHE_DIR; recursive = true, force = true)
println("Test cache directory cleaned up.")
