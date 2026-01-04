using Graphs
using GraphPlot

struct DecayTopology{G}
    topology::G
    names::Dict{Any,Int}
end

node2string(n) = replace(string(n), "Any" => "")

ordered_node_names(dt::DecayTopology) =
    first.(sort(collect(dt.names), by=x -> x[2]))

# build graph from bracket representation,
# Ex: [[1,2],3]
# 
#       0
#       |
#   [[1,2],3]
#      /  \
#     /    3
#   [1,2]
#   /  \
#  /    2
# 1
function parse_topology!(g::DiGraph, node_map::Dict{Any,Int}, topology, parent::Int=0)
    # This function recursively processes the topology
    # If it's an array, it needs to create a new vertex for this node
    if topology isa Array
        # Create a new node for the current structure
        add_vertex!(g)
        current_node = nv(g)
        node_map[topology] = current_node

        # Link it to its parent if there is one
        if parent != 0
            add_edge!(g, parent, current_node)
        end

        # Process each child
        for item in topology
            parse_topology!(g, node_map, item, current_node)
        end
    else
        # It's a leaf node, check if it exists, if not, add it
        if !haskey(node_map, topology)
            add_vertex!(g)
            node_map[topology] = nv(g)
        end
        # Link leaf to the current parent node
        if parent != 0
            add_edge!(g, parent, node_map[topology])
        end
    end
end

function DecayTopology(topology::Array{Any,1})
    g = DiGraph()
    node_map = Dict{Any,Int}()

    parse_topology!(g, node_map, topology)
    DecayTopology(g, node_map)
end

using Test
# test
dt = DecayTopology([[1, 2], 3])
@test Set(node2string.(keys(dt.names))) == Set(["[1, 2]", "1", "2", "3", "[[1, 2], 3]"])
# compute adjacency matrix
@test adjacency_matrix(g) ==
      [0 1 0 0 1
    0 0 1 1 0
    0 0 0 0 0
    0 0 0 0 0
    0 0 0 0 0]

dt = DecayTopology([[2, 1], 3])
@test Set(node2string.(keys(dt.names))) == Set(["[2, 1]", "1", "2", "3", "[[2, 1], 3]"])
@test adjacency_matrix(g) ==
      [0 1 0 0 1
    0 0 1 1 0
    0 0 0 0 0
    0 0 0 0 0
    0 0 0 0 0]

let
    names = ordered_node_names(dt)
    nodelabel = node2string.(names)
    gplot(dt.topology; nodelabel)
end

dt = DecayTopology([[[1, 2], 3], 4])
gplot(g; nodelabel=labels(dt))

g = dt.topology
neighbors(g, 3)

# # get children of a node
# function decay_products(dt::DecayTopology, node_name)
#     index = dt.names[node_name]
#     nei = neighbors(g, index)
#     invnames = Dict(v => k for (k, v) in dt.names)
#     getindex.(Ref(invnames), collect(nei))
# end

# decay_products(dt::DecayTopology, [2, 1])
# decay_products(dt::DecayTopology, [[2, 1], 3]) == [2, 1], 3
