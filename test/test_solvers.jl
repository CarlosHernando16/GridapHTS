using GridapHTS
using GridapHTS: add_default_params, validate_params, get_nested, merge_nested,
    _symbolize_keys, norm_safe
using Gridap
using Test

@testset "Solvers & Utilities" begin

    @testset "Parameter defaults" begin
        # Empty params get all defaults
        params = add_default_params(Dict())
        @test params[:formulation] == :A
        @test params[:fe_order] == 1
        @test params[:mesh][:type] == :cartesian
        @test params[:material][:jc] == GridapHTS.JC_DEFAULT
        @test params[:solver][:type] == :newton
    end

    @testset "Parameter override" begin
        user_params = Dict(
            :formulation => :TA,
            :material => Dict(:jc => 5e10, :n_exponent => 30),
        )
        params = add_default_params(user_params)
        @test params[:formulation] == :TA
        @test params[:material][:jc] == 5e10
        @test params[:material][:n_exponent] == 30
        # Defaults preserved for unspecified keys
        @test params[:material][:ec] == GridapHTS.E_C_DEFAULT
        @test params[:solver][:max_iter] == 50
    end

    @testset "Parameter validation" begin
        # Valid params
        params = add_default_params(Dict())
        @test validate_params(params) === nothing

        # Invalid formulation
        bad_params = deepcopy(params)
        bad_params[:formulation] = :invalid
        @test_throws ArgumentError validate_params(bad_params)

        # Invalid J_c
        bad_mat = deepcopy(params)
        bad_mat[:material][:jc] = -1.0
        @test_throws ArgumentError validate_params(bad_mat)

        # Gmsh without file
        bad_mesh = deepcopy(params)
        bad_mesh[:mesh][:type] = :gmsh
        bad_mesh[:mesh][:file] = nothing
        @test_throws ArgumentError validate_params(bad_mesh)

        # Axisymmetric cartesian mesh with negative r is invalid
        bad_axi = deepcopy(params)
        bad_axi[:mesh][:coordinate_system] = :axisymmetric2d
        bad_axi[:mesh][:domain] = (-1.0, 1.0, 0.0, 1.0)
        @test_throws ArgumentError validate_params(bad_axi)
    end

    @testset "String key symbolization" begin
        d = Dict("a" => 1, "b" => Dict("c" => 2))
        sd = _symbolize_keys(d)
        @test haskey(sd, :a)
        @test haskey(sd, :b)
        @test sd[:b][:c] == 2
    end

    @testset "Nested get helper" begin
        params = Dict(:a => Dict(:b => Dict(:c => 42)))
        @test get_nested(params, :a, :b, :c) == 42
        @test get_nested(params, :a, :b, :d; default=99) == 99
        @test get_nested(params, :x; default=nothing) === nothing
    end

    @testset "Nested merge" begin
        base = Dict(:a => 1, :b => Dict(:c => 2, :d => 3))
        override = Dict(:b => Dict(:c => 10, :e => 5))
        result = merge_nested(base, override)
        @test result[:a] == 1
        @test result[:b][:c] == 10
        @test result[:b][:d] == 3
        @test result[:b][:e] == 5
    end

    @testset "norm_safe" begin
        v = VectorValue(3.0, 4.0)
        @test isapprox(norm_safe(v; eps=0.0), 5.0, atol=1e-10)

        # Zero vector with regularization
        v_zero = VectorValue(0.0, 0.0)
        @test norm_safe(v_zero) > 0.0
        @test isfinite(norm_safe(v_zero))

        # Scalar version
        @test isapprox(norm_safe(3.0; eps=0.0), 3.0, atol=1e-10)
        @test norm_safe(0.0) > 0.0
    end

    @testset "main() driver - A-formulation" begin
        params = Dict(
            :formulation => :A,
            :fe_order => 1,
            :mesh => Dict(
                :type => :cartesian,
                :domain => (0, 1, 0, 1),
                :partition => (5, 5),
            ),
        )

        xh, fullparams, output = GridapHTS.main(params)
        @test fullparams[:formulation] == :A
        @test haskey(output, :elapsed_time)
        @test output[:elapsed_time] > 0.0
    end

end
