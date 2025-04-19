using Pkg

# Shared configuration that matches CI
const PLUTO_VERSION = "0.3.2-0.3"
const CACHE_DIR = "test_pluto_state_cache"

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
)

println("\nTest completed! If no errors occurred, notebooks should work in CI.")
println("You can find the test cache in: ", CACHE_DIR)

# Cleanup
rm(CACHE_DIR; recursive = true, force = true)
println("Test cache directory cleaned up.")
