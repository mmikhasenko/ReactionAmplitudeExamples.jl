module NotebookIndexer

export generate_index_html

"""
    is_pluto_notebook(path2file::String)

Check if a file is a Pluto notebook by examining its first line.
"""
function is_pluto_notebook(path2file)
    try
        l = readline(path2file)
        return l == "### A Pluto.jl notebook ###"
    catch e
        @warn "Could not read file $path2file: $e"
        return false
    end
end

"""
    extract_md_blocks(content::AbstractString)

Extract all markdown blocks from a Pluto notebook content.
Returns an array of strings containing the markdown content.
"""
function extract_md_blocks(content::AbstractString)
    r = r"md\"\"\"\r?\n(.*?)\"\"\""s
    return [m.captures[1] for m in eachmatch(r, content)]
end

"""
    parse_heading(heading::AbstractString)

Parse a markdown heading to extract its level and title.
Returns a tuple of (level::Int, title::String).
"""
function parse_heading(heading::AbstractString)
    firstheading = split(heading, '\n') |> first
    m = match(r"^(#+) (.*)\r?$", firstheading)
    if isnothing(m)
        return 1, "Untitled"  # Default values for invalid headings
    end
    level = length(m.captures[1])
    title = strip(m.captures[2])
    return level, title
end

"""
    build_index_structure(content::AbstractString)

Build a hierarchical structure of headings from notebook content.
Returns an array of Pairs representing the document structure.
"""
function build_index_structure(content::AbstractString)
    # Extract all markdown blocks
    blocks = extract_md_blocks(content)
    filter!(blocks) do block
        occursin(r"^(#+) \r?.*$"s, block)
    end
    # Extract headings from blocks and their levels
    headings = parse_heading.(blocks)

    # Initialize root of the structure
    root = Pair("root", [])
    current_struct = [root]

    for (level, title) in headings
        while length(current_struct) < level
            # Fill missing levels with "Untitled"
            new_untitled = Pair("Untitled", [])
            push!(last(current_struct).second, new_untitled)
            push!(current_struct, new_untitled)
        end
        while length(current_struct) > level
            pop!(current_struct)
        end

        # Add current title to the structure
        new_entry = Pair(title, [])
        push!(last(current_struct).second, new_entry)

        if level < length(current_struct)
            current_struct[level] = new_entry
        else
            push!(current_struct, new_entry)
        end
    end

    return root.second
end

"""
    extend_itemized_list(notebooks, html_template)

Generate HTML content for the notebook index by processing each notebook
and extracting its structure.
"""
function extend_itemized_list(notebooks, html_template)
    println("Processing $(length(notebooks)) notebooks...")
    index_end = first(findfirst("</ul>", html_template)) - 1
    list_items = ""
    untitled_count = 0

    for (i, notebook) in enumerate(notebooks)
        println("[$i/$(length(notebooks))] Processing $notebook...")
        notebook_html = replace(notebook, ".jl" => ".html")
        index = build_index_structure(read(notebook, String))

        # Extract titles and subtitles from the index
        if length(index) == 0
            index = [notebook => []]
        end
        length(index) != 1 && error("Fishy index, not `Title=>[Subs...] in $(notebook)`\n\nindex = $(index)\n\nlength(index) = $(length(index))")
        title = first(index[1])
        subtitles = last(index[1])

        # Check if notebook is untitled
        is_untitled = title == "Untitled" || title == notebook
        if is_untitled
            untitled_count += 1
            println("⚠️  Warning: $notebook has no title. Please add a markdown cell with a # Title and brief description.")
        end

        subtitle_list_items = ""
        for subtitle in subtitles
            subtitle_name = first(subtitle)
            subtitle_list_items *= """<li>$(subtitle_name)</li>"""
        end

        # Add warning style for untitled notebooks
        title_class = is_untitled ? "untitled-warning" : "title"
        warning_msg = is_untitled ? """
            <span class="warning-message">
                ⚠️ This notebook needs a title and description. 
                Please add a markdown cell with a # Title and brief description.
            </span>""" : ""

        # incorporate to html
        list_items *= """
        <li>
            <a href="notebooks/$(notebook_html)" target="notebook-frame">$(notebook)</a>
            <div class="notebook-details">
                <span class="$(title_class)">$(title)</span>
                $(warning_msg)
                <ul class="subtitles">$(subtitle_list_items)</ul>
            </div>
        </li>
        """
    end

    println("Generating HTML output...")
    if untitled_count > 0
        println("\n⚠️  Found $untitled_count untitled notebook(s). Please add titles and descriptions to improve navigation.")
    end

    extended_html_string = string(html_template[1:index_end-1], list_items, html_template[index_end:end])
    return extended_html_string
end

"""
    generate_index_html()

Main function to generate the index.html file for the notebook collection.
"""
function generate_index_html()
    println("\n=== Notebook Index Generator ===")

    # Verify we're in the notebooks directory
    notebooks_dir = joinpath(@__DIR__, "..", "notebooks")
    isdir(notebooks_dir) || error("Notebooks directory not found: $notebooks_dir")

    # Get all Pluto notebooks
    try
        println("\nScanning for Pluto notebooks...")
        files = readdir(notebooks_dir)
        notebooks = filter(f -> is_pluto_notebook(joinpath(notebooks_dir, f)), files)

        if isempty(notebooks)
            @warn "No Pluto notebooks found in directory"
            return
        end
        println("Found $(length(notebooks)) Pluto notebooks")

        # Generate the base HTML template
        println("\nGenerating HTML template...")
        html_template = """
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>$(notebooks_dir)</title>
            <link rel="stylesheet" href="styles.css">
            <style>
                body {
                    margin: 0;
                    padding: 0;
                    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Oxygen, Ubuntu, Cantarell, "Open Sans", "Helvetica Neue", sans-serif;
                }
                .container {
                    width: 100%;
                    height: 100vh;
                }
                .index-panel {
                    position: relative;
                    width: 300px;
                    min-width: 50px;
                    height: 100vh;
                    overflow-y: auto;
                    float: left;
                    padding: 20px;
                    box-sizing: border-box;
                    background-color: #f8f9fa;
                    border-right: 1px solid #dee2e6;
                    transition: width 0.3s;
                }
                .index-panel.collapsed {
                    width: 30px;
                    padding: 10px 5px;
                }
                .index-panel .content {
                    opacity: 1;
                    transition: opacity 0.3s;
                    width: 100%;
                }
                .index-panel.collapsed .content {
                    opacity: 0;
                }
                .resize-handle {
                    position: absolute;
                    right: 0;
                    top: 0;
                    width: 4px;
                    height: 100%;
                    background-color: rgba(0, 0, 0, 0.1);
                    cursor: col-resize;
                }
                .resize-handle:hover {
                    background-color: rgba(0, 0, 0, 0.2);
                }
                .toggle-button {
                    position: absolute;
                    right: -15px;
                    top: 50%;
                    transform: translateY(-50%);
                    width: 15px;
                    height: 30px;
                    background-color: #f8f9fa;
                    border: 1px solid #dee2e6;
                    border-left: none;
                    border-radius: 0 4px 4px 0;
                    cursor: pointer;
                    z-index: 1000;
                    display: flex;
                    align-items: center;
                    justify-content: center;
                }
                .toggle-button::after {
                    content: "◄";
                    font-size: 10px;
                }
                .index-panel.collapsed .toggle-button::after {
                    content: "►";
                }
                .iframe-container {
                    margin-left: 300px;
                    height: 100vh;
                    transition: margin-left 0.3s;
                }
                .index-panel.collapsed + .iframe-container {
                    margin-left: 30px;
                }
                .notebook-display {
                    width: 100%;
                    height: 100%;
                    border: none;
                }
                .untitled-warning {
                    color: #856404;
                    background-color: #fff3cd;
                    padding: 2px 8px;
                    border-radius: 4px;
                    margin-right: 8px;
                }
                .warning-message {
                    color: #856404;
                    font-size: 0.9em;
                    font-style: italic;
                    display: block;
                    margin-top: 4px;
                }
                .notebook-details {
                    margin-left: 20px;
                    margin-bottom: 10px;
                }
                .subtitles {
                    margin: 5px 0;
                    padding-left: 20px;
                    list-style-type: disc;
                }
                .title {
                    font-weight: bold;
                }
                h3 {
                    margin-top: 0;
                }
            </style>
        </head>
        <body>
            <div class="container">
                <div class="index-panel">
                    <div class="resize-handle"></div>
                    <button class="toggle-button"></button>
                    <div class="content">
                        <h3>Index</h3>
                        <ul>
                            <li><a href="README.html" target="notebook-frame">README</a></li>
                        </ul>
                    </div>
                </div>
                <div class="iframe-container">
                    <iframe name="notebook-frame" class="notebook-display" src="README.html"></iframe>
                </div>
            </div>
            <script>
                // Sidebar resize functionality
                const resizeHandle = document.querySelector('.resize-handle');
                const sidebar = document.querySelector('.index-panel');
                const container = document.querySelector('.container');
                const iframeContainer = document.querySelector('.iframe-container');
                let isResizing = false;

                function handleMouseMove(e) {
                    if (!isResizing) return;
                    e.preventDefault();
                    // Calculate new width based on mouse position
                    const newWidth = Math.max(50, Math.min(600, e.clientX));
                    sidebar.style.width = newWidth + 'px';
                    iframeContainer.style.marginLeft = newWidth + 'px';
                }

                function stopResize(e) {
                    if (!isResizing) return;
                    isResizing = false;
                    window.removeEventListener('mousemove', handleMouseMove);
                    window.removeEventListener('mouseup', stopResize);
                    window.removeEventListener('mouseleave', stopResize);
                }

                resizeHandle.addEventListener('mousedown', (e) => {
                    e.preventDefault();
                    isResizing = true;
                    window.addEventListener('mousemove', handleMouseMove);
                    window.addEventListener('mouseup', stopResize);
                    window.addEventListener('mouseleave', stopResize);
                });

                // Sidebar collapse functionality
                const toggleButton = document.querySelector('.toggle-button');
                
                toggleButton.addEventListener('click', () => {
                    sidebar.classList.toggle('collapsed');
                    const isCollapsed = sidebar.classList.contains('collapsed');
                    iframeContainer.style.marginLeft = isCollapsed ? '30px' : '300px';
                    
                    // Store the current state
                    localStorage.setItem('sidebarCollapsed', isCollapsed);
                });

                // Restore sidebar state on page load
                window.addEventListener('load', () => {
                    const isCollapsed = localStorage.getItem('sidebarCollapsed') === 'true';
                    if (isCollapsed) {
                        sidebar.classList.add('collapsed');
                        iframeContainer.style.marginLeft = '30px';
                    }
                });
            </script>
        </body>
        </html>
        """

        # Generate and write the index
        extended_html = extend_itemized_list(joinpath.(notebooks_dir, notebooks), html_template)
        output_path = joinpath(@__DIR__, "..", "index.html")
        write(output_path, extended_html)
        println("\n✓ Successfully generated index at $output_path")
        return output_path

    catch e
        @error "Failed to generate index" exception = e
        rethrow(e)
    end
end

end # module

# Run the index generation when script is executed directly
NotebookIndexer.generate_index_html()