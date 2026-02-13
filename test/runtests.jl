using GridapHTS
using Test

@testset "GridapHTS.jl" begin
    include("test_materials.jl")
    include("test_formulations.jl")
    include("test_solvers.jl")
end
