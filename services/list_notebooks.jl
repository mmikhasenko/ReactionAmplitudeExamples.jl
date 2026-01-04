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
    html_escape(s::AbstractString)

Escape special HTML characters for use in HTML attributes.
"""
function html_escape(s::AbstractString)
    s = replace(s, "&" => "&amp;")
    s = replace(s, "\"" => "&quot;")
    s = replace(s, "'" => "&#39;")
    s = replace(s, "<" => "&lt;")
    s = replace(s, ">" => "&gt;")
    return s
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
        # Extract just the filename from the full path
        notebook_filename = basename(notebook)
        notebook_html = replace(notebook_filename, ".jl" => ".html")
        index = build_index_structure(read(notebook, String))

        # Extract titles and subtitles from the index
        if length(index) == 0
            index = [notebook_filename => []]
        end
        length(index) != 1 && error("Fishy index, not `Title=>[Subs...] in $(notebook)`\n\nindex = $(index)\n\nlength(index) = $(length(index))")
        title = first(index[1])
        subtitles = last(index[1])

        # Check if notebook is untitled
        is_untitled = title == "Untitled" || title == notebook_filename
        if is_untitled
            untitled_count += 1
            println("⚠️  Warning: $notebook_filename has no title. Please add a markdown cell with a # Title and brief description.")
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

        # Prepare searchable text for search functionality
        subtitle_text = join([s.first for s in subtitles], ", ")
        search_description = isempty(subtitle_text) ? title : subtitle_text

        # incorporate to html
        list_items *= """
        <li>
            <a href="notebooks/$(notebook_html)" target="notebook-frame" data-title="$(html_escape(title))" data-description="$(html_escape(search_description))">$(notebook_filename)</a>
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
        # MathJax configuration (using raw string to avoid $ escaping issues)
        mathjax_config = raw"""<script>
window.MathJax = {
    tex: {
        inlineMath: [['$', '$'], ['\\(', '\\)']],
        displayMath: [['$$', '$$'], ['\\[', '\\]']],
        processEscapes: true,
        processEnvironments: true
    },
    options: {
        skipHtmlTags: ['script', 'noscript', 'style', 'textarea', 'pre']
    }
};
</script>
<script src="https://polyfill.io/v3/polyfill.min.js?features=es6"></script>
<script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>"""
        html_template = """
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>Pluto Notebooks - Reaction Amplitude Examples</title>
            <link rel="stylesheet" href="styles.css">
            <!-- MathJax for LaTeX rendering -->
            """ * mathjax_config * """
        </head>
        <body>
            <div class="container">
                <div class="index-panel">
                    <div class="resize-handle"></div>
                    <button class="toggle-button" aria-label="Toggle sidebar"></button>
                    <div class="content">
                        <div class="sidebar-header">
                            <h3>Notebook Index</h3>
                            <p>Browse and search $(length(notebooks)) notebooks</p>
                        </div>
                        <div class="search-container">
                            <span class="search-icon">🔍</span>
                            <input type="text" class="search-box" id="searchBox" placeholder="Search notebooks..." autocomplete="off">
                            <button class="search-clear" id="searchClear" aria-label="Clear search">×</button>
                        </div>
                        <ul id="notebookList">
                            <li><a href="README.html" target="notebook-frame" data-title="README" data-description="Project documentation">README</a></li>
                        </ul>
                    </div>
                </div>
                <div class="iframe-container">
                    <iframe name="notebook-frame" class="notebook-display" src="README.html" title="Notebook viewer" onerror="handleIframeError(this)"></iframe>
                </div>
            </div>
            <script>
                (function() {
                // Sidebar resize functionality
                const resizeHandle = document.querySelector('.resize-handle');
                const sidebar = document.querySelector('.index-panel');
                const iframeContainer = document.querySelector('.iframe-container');
                let isResizing = false;

                function handleMouseMove(e) {
                    if (!isResizing) return;
                    e.preventDefault();
                    const newWidth = Math.max(50, Math.min(600, e.clientX));
                    sidebar.style.width = newWidth + 'px';
                }

                function stopResize(e) {
                    if (!isResizing) return;
                    isResizing = false;
                    window.removeEventListener('mousemove', handleMouseMove);
                    window.removeEventListener('mouseup', stopResize);
                    window.removeEventListener('mouseleave', stopResize);
                }

                    if (resizeHandle) {
                resizeHandle.addEventListener('mousedown', (e) => {
                    e.preventDefault();
                    isResizing = true;
                    window.addEventListener('mousemove', handleMouseMove);
                    window.addEventListener('mouseup', stopResize);
                    window.addEventListener('mouseleave', stopResize);
                });
                    }

                // Sidebar collapse functionality
                const toggleButton = document.querySelector('.toggle-button');
                
                    if (toggleButton) {
                toggleButton.addEventListener('click', () => {
                    sidebar.classList.toggle('collapsed');
                    const isCollapsed = sidebar.classList.contains('collapsed');
                    localStorage.setItem('sidebarCollapsed', isCollapsed);
                });
                    }

                // Restore sidebar state on page load
                window.addEventListener('load', () => {
                    const isCollapsed = localStorage.getItem('sidebarCollapsed') === 'true';
                        if (isCollapsed && sidebar) {
                        sidebar.classList.add('collapsed');
                        }
                    });

                    // Search functionality
                    const searchBox = document.getElementById('searchBox');
                    const searchClear = document.getElementById('searchClear');
                    const notebookList = document.getElementById('notebookList');
                    const allNotebookItems = [];

                    // Collect all notebook items (excluding README)
                    function collectNotebookItems() {
                        const items = notebookList.querySelectorAll('li:not(:first-child)');
                        items.forEach(item => {
                            const link = item.querySelector('a');
                            if (link) {
                                const title = link.getAttribute('data-title') || link.textContent;
                                const description = link.getAttribute('data-description') || '';
                                const filename = link.textContent;
                                const details = item.querySelector('.notebook-details');
                                const subtitles = details ? Array.from(details.querySelectorAll('.subtitles li')).map(li => li.textContent).join(' ') : '';
                                
                                allNotebookItems.push({
                                    element: item,
                                    searchText: (filename + ' ' + title + ' ' + description + ' ' + subtitles).toLowerCase()
                                });
                            }
                        });
                    }

                    function filterNotebooks(searchTerm) {
                        const term = searchTerm.toLowerCase().trim();
                        let visibleCount = 0;

                        allNotebookItems.forEach(({element, searchText}) => {
                            if (term === '' || searchText.includes(term)) {
                                element.classList.remove('no-match');
                                visibleCount++;
                            } else {
                                element.classList.add('no-match');
                            }
                        });

                        // Show/hide clear button
                        if (searchClear) {
                            if (term.length > 0) {
                                searchClear.classList.add('visible');
                            } else {
                                searchClear.classList.remove('visible');
                            }
                        }

                        // Show empty state if no results
                        let emptyState = notebookList.querySelector('.empty-state');
                        if (visibleCount === 0 && term.length > 0) {
                            if (!emptyState) {
                                emptyState = document.createElement('li');
                                emptyState.className = 'empty-state';
                                emptyState.innerHTML = '<div class="empty-state-icon">📚</div><div class="empty-state-text">No notebooks found matching "' + term + '"</div>';
                                notebookList.appendChild(emptyState);
                            }
                        } else if (emptyState) {
                            emptyState.remove();
                        }
                    }

                    if (searchBox) {
                        searchBox.addEventListener('input', (e) => {
                            filterNotebooks(e.target.value);
                        });

                        // Keyboard shortcuts
                        searchBox.addEventListener('keydown', (e) => {
                            if (e.key === 'Escape') {
                                searchBox.value = '';
                                filterNotebooks('');
                                searchBox.blur();
                            }
                        });
                    }

                    if (searchClear) {
                        searchClear.addEventListener('click', () => {
                            if (searchBox) {
                                searchBox.value = '';
                                filterNotebooks('');
                                searchBox.focus();
                            }
                        });
                    }

                    // Active link highlighting
                    function setActiveLink(href) {
                        const links = notebookList.querySelectorAll('a');
                        links.forEach(link => {
                            if (link.getAttribute('href') === href) {
                                link.classList.add('active');
                            } else {
                                link.classList.remove('active');
                            }
                        });
                    }

                    // Update active link when iframe loads
                    const iframe = document.querySelector('iframe[name="notebook-frame"]');
                    if (iframe) {
                        iframe.addEventListener('load', () => {
                            try {
                                const iframeSrc = iframe.contentWindow.location.pathname;
                                const relativePath = iframeSrc.split('/').pop() || 'README.html';
                                setActiveLink(relativePath);
                            } catch (e) {
                                // Cross-origin or other error, use src attribute
                                const src = iframe.getAttribute('src');
                                if (src) {
                                    setActiveLink(src);
                                }
                            }
                        });
                    }

                    // Handle clicks on notebook links
                    notebookList.addEventListener('click', (e) => {
                        const link = e.target.closest('a');
                        if (link) {
                            setActiveLink(link.getAttribute('href'));
                        }
                    });

                    // Initialize
                    collectNotebookItems();
                    
                    // Handle iframe load errors
                    function handleIframeError(iframe) {
                        // This will be called if iframe fails to load
                        console.warn('Iframe failed to load:', iframe.src);
                    }
                    
                    // Check if iframe loaded successfully
                    const iframe = document.querySelector('iframe[name="notebook-frame"]');
                    if (iframe) {
                        iframe.addEventListener('load', function() {
                            try {
                                // Try to access iframe content to check if it loaded
                                const iframeDoc = iframe.contentDocument || iframe.contentWindow.document;
                                // If we can access it, check for 404
                                if (iframeDoc.body && iframeDoc.body.textContent.includes('404')) {
                                    showNotFoundMessage(iframe.src);
                                }
                            } catch (e) {
                                // Cross-origin or other error - might be 404
                                // Check after a delay
                                setTimeout(() => {
                                    try {
                                        const iframeDoc = iframe.contentDocument || iframe.contentWindow.document;
                                        if (!iframeDoc || !iframeDoc.body) {
                                            showNotFoundMessage(iframe.src);
                                        }
                                    } catch (err) {
                                        // Can't access - might be cross-origin, which is fine
                                    }
                                }, 1000);
                            }
                        });
                        
                        iframe.addEventListener('error', function() {
                            showNotFoundMessage(iframe.src);
                        });
                    }
                    
                    function showNotFoundMessage(url) {
                        const iframe = document.querySelector('iframe[name="notebook-frame"]');
                        if (iframe) {
                            const filename = url.split('/').pop() || url;
                            const errorHtml = '<!DOCTYPE html><html><head><meta charset="UTF-8"><style>body{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif;display:flex;align-items:center;justify-content:center;height:100vh;margin:0;background:#f8f9fa;color:#4a5568}.error-container{text-align:center;padding:40px;max-width:600px}.error-icon{font-size:64px;margin-bottom:20px}h1{color:#1a202c;margin-bottom:16px}p{line-height:1.6;margin-bottom:12px}code{background:#e2e8f0;padding:2px 6px;border-radius:3px;font-family:monospace}</style></head><body><div class="error-container"><div class="error-icon">📄</div><h1>Notebook Not Found</h1><p>The notebook <code>' + filename + '</code> could not be found.</p><p><strong>Local Testing:</strong> Notebook HTML files are generated by PlutoSliderServer in the CI workflow. They won\'t exist locally unless you export them manually.</p><p><strong>On GitHub Pages:</strong> If you see this on the live site, the notebook may not have been exported successfully.</p></div></body></html>';
                            iframe.srcdoc = errorHtml;
                        }
                    }
                })();
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