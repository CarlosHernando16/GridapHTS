# ──────────────────────────────────────────────
# Utility Functions
# ──────────────────────────────────────────────

"""
    norm_safe(v::VectorValue; eps=EPS_REG)

Compute the Euclidean norm of a `VectorValue` with regularization
to avoid division-by-zero issues in power-law computations.

# Arguments
- `v::VectorValue`: Input vector
- `eps`: Regularization parameter (default: `EPS_REG`)

# Returns
- `Float64`: Regularized norm, `sqrt(v⋅v + eps)`
"""
function norm_safe(v::VectorValue; eps=EPS_REG)
    return sqrt(v ⋅ v + eps)
end

"""
    norm_safe(s::Real; eps=EPS_REG)

Scalar version: returns `sqrt(s^2 + eps)`.
"""
function norm_safe(s::Real; eps=EPS_REG)
    return sqrt(s^2 + eps)
end

"""
    get_nested(params::Dict, keys...; default=nothing)

Safely retrieve a nested value from a parameter dictionary.

# Example
```julia
params = Dict(:solver => Dict(:type => :newton, :rtol => 1e-8))
get_nested(params, :solver, :type)  # => :newton
get_nested(params, :solver, :atol, default=1e-10)  # => 1e-10
```
"""
function get_nested(params::Dict, keys...; default=nothing)
    d = params
    for k in keys
        if d isa Dict && haskey(d, k)
            d = d[k]
        else
            return default
        end
    end
    return d
end

"""
    merge_nested(base::Dict, override::Dict)

Deep-merge two nested dictionaries. Values in `override` take precedence.
"""
function merge_nested(base::Dict, override::Dict)
    result = copy(base)
    for (k, v) in override
        if haskey(result, k) && result[k] isa Dict && v isa Dict
            result[k] = merge_nested(result[k], v)
        else
            result[k] = v
        end
    end
    return result
end

"""
    ensure_directory(path::AbstractString)

Create directory (and parents) if it does not exist. Returns the path.
"""
function ensure_directory(path::AbstractString)
    mkpath(path)
    return path
end
