using GridapHTS
using GridapHTS: a_bilinear_form, a_linear_form, grad_grad_bilinear,
    setup_a_formulation, solve_a_formulation, _build_model
using Gridap
using Test

@testset "Formulations" begin

    @testset "Poisson equation validation" begin
        # Solve -Δu = f on [0,1]² with u=0 on boundary
        # Analytical solution: u = sin(πx)sin(πy), f = 2π²sin(πx)sin(πy)
        domain = (0, 1, 0, 1)
        partition = (20, 20)
        model = CartesianDiscreteModel(domain, partition)

        order = 1
        reffe = ReferenceFE(lagrangian, Float64, order)
        V = TestFESpace(model, reffe; conformity=:H1, dirichlet_tags="boundary")
        U = TrialFESpace(V, 0.0)

        Ω = Triangulation(model)
        dΩ = Measure(Ω, 2 * order + 1)

        # Weak forms
        a(u, v) = ∫(∇(u) ⋅ ∇(v))dΩ
        f(x) = 2π^2 * sin(π * x[1]) * sin(π * x[2])
        l(v) = ∫(f * v)dΩ

        op = AffineFEOperator(a, l, U, V)
        uh = solve(op)

        # Compute L² error
        u_exact(x) = sin(π * x[1]) * sin(π * x[2])
        e = u_exact - uh
        l2_error = sqrt(sum(∫(e * e)dΩ))

        # Error should be small (order-dependent)
        @test l2_error < 0.05
    end

    @testset "A-formulation setup (2D Cartesian)" begin
        params = Dict(
            :formulation => :A,
            :fe_order => 1,
            :mesh => Dict(
                :type => :cartesian,
                :domain => (0, 1, 0, 1),
                :partition => (10, 10),
            ),
            :material => Dict(
                :mu => GridapHTS.MU_0,
            ),
            :bcs => Dict(
                :dirichlet_tags => ["boundary"],
            ),
        )

        fullparams = add_default_params(params)
        setup = setup_a_formulation(fullparams)

        @test setup.model isa Gridap.Geometry.DiscreteModel
        @test setup.Ω isa Gridap.Geometry.Triangulation
        @test setup.μ_inv ≈ 1.0 / GridapHTS.MU_0
    end

    @testset "A-formulation solve (zero source)" begin
        params = Dict(
            :formulation => :A,
            :fe_order => 1,
            :mesh => Dict(
                :type => :cartesian,
                :domain => (0, 1, 0, 1),
                :partition => (10, 10),
            ),
            :material => Dict(:mu => GridapHTS.MU_0),
            :bcs => Dict(:dirichlet_tags => ["boundary"]),
        )

        fullparams = add_default_params(params)
        uh, setup = solve_a_formulation(fullparams)

        # With zero source and zero Dirichlet BC, solution should be ≈ 0
        Ω = setup.Ω
        dΩ = setup.dΩ
        l2_norm = sqrt(sum(∫(uh * uh)dΩ))
        @test l2_norm < 1e-10
    end

    @testset "A-formulation solve (axisymmetric 2D, zero source)" begin
        params = Dict(
            :formulation => :A,
            :fe_order => 1,
            :mesh => Dict(
                :type => :cartesian,
                :coordinate_system => :axisymmetric2d,
                :domain => (0.0, 1.0, 0.0, 1.0),
                :partition => (10, 10),
            ),
            :material => Dict(:mu => GridapHTS.MU_0),
            :bcs => Dict(:dirichlet_tags => ["boundary"]),
        )

        fullparams = add_default_params(params)
        uh, setup = solve_a_formulation(fullparams)

        @test setup.coordinate_system == :axisymmetric2d

        # With zero source and homogeneous Dirichlet BC, solution should remain ≈ 0.
        Ω = setup.Ω
        dΩ = setup.dΩ
        l2_norm = sqrt(sum(∫(uh * uh)dΩ))
        @test l2_norm < 1e-10
    end

    @testset "Model builder" begin
        params_cart = Dict(
            :mesh => Dict(
                :type => :cartesian,
                :domain => (0, 1, 0, 1),
                :partition => (5, 5),
            ),
        )
        model = _build_model(params_cart)
        @test model isa CartesianDiscreteModel

        # Gmsh without file should error
        params_gmsh = Dict(
            :mesh => Dict(:type => :gmsh, :file => nothing),
        )
        @test_throws Exception _build_model(params_gmsh)
    end

end
