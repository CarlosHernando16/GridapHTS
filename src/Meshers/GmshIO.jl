# ──────────────────────────────────────────────
# Gmsh Mesh I/O Utilities
# ──────────────────────────────────────────────
#
# Helpers for loading and working with Gmsh meshes via GridapGmsh.
# ──────────────────────────────────────────────

"""
    load_gmsh_model(mesh_file::AbstractString) -> DiscreteModel

Load a Gmsh .msh file and return a Gridap `DiscreteModel`.

# Arguments
- `mesh_file::AbstractString`: Path to .msh file

# Returns
- `model::DiscreteModel`: Gridap discrete model

# Example
```julia
model = load_gmsh_model(datadir("meshes", "tape_2d", "rectangular.msh"))
```

# Notes
- Expects Gmsh format version 4.x
- Physical groups in the .msh file become Gridap boundary/domain tags
"""
function load_gmsh_model(mesh_file::AbstractString)
    isfile(mesh_file) || throw(ArgumentError(
        "Mesh file not found: $mesh_file"))

    return GmshDiscreteModel(mesh_file)
end

"""
    setup_cartesian_model(domain, partition) -> DiscreteModel

Create a simple Cartesian (structured) discrete model.

# Arguments
- `domain`: Tuple defining domain bounds, e.g., `(0, 1, 0, 1)` for 2D
- `partition`: Tuple defining number of cells, e.g., `(20, 20)`

# Returns
- `model::CartesianDiscreteModel`

# Example
```julia
model = setup_cartesian_model((0, 1, 0, 1), (50, 50))
```
"""
function setup_cartesian_model(domain, partition)
    return CartesianDiscreteModel(domain, partition)
end

"""
    mesh_info(model) -> Dict{Symbol,Any}

Extract basic information from a discrete model.

# Returns
Dictionary with:
- `:num_cells`: Number of cells
- `:num_vertices`: Number of vertices
- `:dimension`: Spatial dimension
- `:labels`: Available boundary/domain labels
"""
function mesh_info(model)
    topo = get_grid_topology(model)
    D = num_cell_dims(model)

    info = Dict{Symbol,Any}(
        :dimension    => D,
        :num_cells    => num_cells(model),
        :num_vertices => num_vertices(model),
    )

    # Extract label information
    labels = get_face_labeling(model)
    info[:labels] = labels

    return info
end
