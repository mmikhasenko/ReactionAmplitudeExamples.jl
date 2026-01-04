#!/usr/bin/env julia
# build_and_serve.jl
#
# Combined script to build and serve the local index
# This generates index.html and README.html, then starts a local server

using Pkg
using Sockets
using HTTP

# Activate the project environment
Pkg.activate(".")

const PORT = 8000
const ROOT_DIR = dirname(dirname(abspath(@__FILE__)))

println("=== Build and Serve Script ===")
println("This script generates index.html and README.html, then starts a local server\n")

# Step 1: Generate index.html
println("Step 1: Generating index.html...")
try
    include(joinpath(@__DIR__, "..", "services", "list_notebooks.jl"))
    println("✓ index.html generated successfully\n")
catch e
    @error "Failed to generate index.html" exception = e
    exit(1)
end

# Step 2: Check if pandoc is available and generate README.html
println("Step 2: Generating README.html...")
readme_md = joinpath(ROOT_DIR, "README.md")
readme_html = joinpath(ROOT_DIR, "README.html")

if isfile(readme_md)
    # Try to use pandoc if available, otherwise create a simple HTML wrapper
    try
        run(`pandoc $readme_md -s -o $readme_html`)
        println("✓ README.html generated using pandoc\n")
    catch
        # Fallback: create a simple HTML wrapper
        println("⚠️  pandoc not found, creating simple HTML wrapper...")
        md_content = read(readme_md, String)
        # Simple markdown to HTML conversion (very basic)
        html_content = """
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>README - Reaction Amplitude Examples</title>
            <link rel="stylesheet" href="styles.css">
            <style>
                body {
                    max-width: 800px;
                    margin: 40px auto;
                    padding: 20px;
                    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
                    line-height: 1.6;
                }
                pre {
                    background: #f5f5f5;
                    padding: 15px;
                    border-radius: 4px;
                    overflow-x: auto;
                }
                code {
                    background: #f5f5f5;
                    padding: 2px 6px;
                    border-radius: 3px;
                }
                a {
                    color: #3182ce;
                }
            </style>
        </head>
        <body>
            <pre><code>$(escape_string(md_content))</code></pre>
            <p><em>Note: This is a basic preview. For full markdown rendering, install pandoc.</em></p>
        </body>
        </html>
        """
        write(readme_html, html_content)
        println("✓ README.html created (basic version)\n")
    end
else
    @warn "README.md not found, skipping README.html generation"
end

# Step 3: Start the server
println("=== Starting Local Server ===")
println("Serving from: $ROOT_DIR")
println("Open your browser to: http://localhost:$PORT/index.html")
println("\nPress Ctrl+C to stop the server\n")

function handler(req::HTTP.Request)
    local_path = req.target == "/" ? "/index.html" : req.target

    # Remove leading slash
    file_path = joinpath(ROOT_DIR, lstrip(local_path, '/'))

    # Security: ensure we're serving from the root directory
    file_path = abspath(file_path)
    if !startswith(file_path, ROOT_DIR)
        return HTTP.Response(403, "Forbidden")
    end

    if isfile(file_path)
        content = read(file_path)
        content_type = HTTP.sniff(content)
        return HTTP.Response(200, ["Content-Type" => content_type]; body=content)
    elseif isdir(file_path)
        # Try index.html in directory
        index_path = joinpath(file_path, "index.html")
        if isfile(index_path)
            content = read(index_path)
            return HTTP.Response(200, ["Content-Type" => "text/html"]; body=content)
        end
        return HTTP.Response(404, "Not Found")
    else
        return HTTP.Response(404, "Not Found")
    end
end

try
    HTTP.serve(handler, "127.0.0.1", PORT)
catch e
    if isa(e, InterruptException)
        println("\n\nServer stopped.")
    else
        rethrow(e)
    end
end

