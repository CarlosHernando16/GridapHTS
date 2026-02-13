using GridapHTS
using GridapHTS: resistivity, resistivity_derivative, electric_field,
    critical_current_density, material_from_params
using Gridap
using Test

@testset "Materials" begin

    @testset "PowerLawMaterial construction" begin
        # Default constructor
        mat = PowerLawMaterial()
        @test mat.ec == GridapHTS.E_C_DEFAULT
        @test mat.jc == GridapHTS.JC_DEFAULT
        @test mat.n == GridapHTS.N_DEFAULT

        # Custom constructor
        mat2 = PowerLawMaterial(; ec=1e-4, jc=1e10, n=30)
        @test mat2.jc == 1e10
        @test mat2.n == 30

        # Convenience constructor
        mat3 = PowerLawMaterial(3e10, 25)
        @test mat3.jc == 3e10
        @test mat3.n == 25

        # Invalid parameters
        @test_throws ArgumentError PowerLawMaterial(; jc=-1.0)
        @test_throws ArgumentError PowerLawMaterial(; n=0)
        @test_throws ArgumentError PowerLawMaterial(; ec=0.0)
    end

    @testset "Power-law resistivity" begin
        mat = PowerLawMaterial(; jc=1e9, n=25, ec=1e-4)

        # At J = J_c, ρ should equal E_c / J_c
        rho_jc = resistivity(mat, 1e9)
        @test isapprox(rho_jc, 1e-4 / 1e9, rtol=1e-10)

        # Resistivity should increase with |J|
        rho_low = resistivity(mat, 0.5e9)
        rho_high = resistivity(mat, 2.0e9)
        @test rho_low < rho_jc < rho_high

        # For very small J (regularized), should not be NaN/Inf
        rho_zero = resistivity(mat, 0.0)
        @test isfinite(rho_zero)
        @test rho_zero >= 0.0

        # Higher n gives sharper transition
        mat_low_n = PowerLawMaterial(; jc=1e9, n=5, ec=1e-4)
        mat_high_n = PowerLawMaterial(; jc=1e9, n=50, ec=1e-4)
        # Below J_c: higher n => lower resistivity
        @test resistivity(mat_low_n, 0.5e9) > resistivity(mat_high_n, 0.5e9)
        # Above J_c: higher n => higher resistivity
        @test resistivity(mat_low_n, 2.0e9) < resistivity(mat_high_n, 2.0e9)
    end

    @testset "Power-law resistivity derivative" begin
        mat = PowerLawMaterial(; jc=1e9, n=25, ec=1e-4)

        # Derivative should be positive for positive J
        drho = resistivity_derivative(mat, 1e9)
        @test drho > 0.0

        # Finite difference check
        J0 = 1e9
        δ = 1e3
        rho_plus = resistivity(mat, J0 + δ)
        rho_minus = resistivity(mat, J0 - δ)
        drho_fd = (rho_plus - rho_minus) / (2δ)
        @test isapprox(drho, drho_fd, rtol=1e-3)
    end

    @testset "Electric field vector" begin
        mat = PowerLawMaterial(; jc=1e9, n=25, ec=1e-4)

        J = VectorValue(1e9, 0.0)
        E = electric_field(mat, J)

        # E should be parallel to J
        @test abs(E[2]) < abs(E[1]) * 1e-10

        # |E| ≈ E_c when |J| = J_c
        E_norm = sqrt(E ⋅ E)
        @test isapprox(E_norm, 1e-4, rtol=1e-2)
    end

    @testset "KimModel" begin
        model = KimModel(; jc0=1e10, b0=0.1)

        # Self-field: J_c = J_c0
        B_zero = VectorValue(0.0, 0.0)
        jc_sf = critical_current_density(model, B_zero)
        @test isapprox(jc_sf, 1e10, rtol=1e-5)

        # At B = B_0: J_c = J_c0 / 2
        B_b0 = VectorValue(0.0, 0.1)
        jc_b0 = critical_current_density(model, B_b0)
        @test isapprox(jc_b0, 0.5e10, rtol=1e-2)

        # J_c decreases with increasing B
        B_high = VectorValue(0.0, 1.0)
        jc_high = critical_current_density(model, B_high)
        @test jc_high < jc_b0 < jc_sf

        # Scalar version
        jc_scalar = critical_current_density(model, 0.1)
        @test isapprox(jc_scalar, 0.5e10, rtol=1e-2)
    end

    @testset "MaterialLibrary" begin
        rebco = rebco_default()
        @test rebco.jc == 3e10
        @test rebco.n == 25

        bscco = bscco_default()
        @test bscco.jc == 1e9
        @test bscco.n == 15

        # Custom overrides
        rebco2 = rebco_default(; jc=5e10, n=30)
        @test rebco2.jc == 5e10
        @test rebco2.n == 30
    end

    @testset "material_from_params" begin
        # Power-law
        params_pl = Dict(:material => Dict(:type => :power_law, :jc => 2e9, :n_exponent => 20))
        mat_pl = material_from_params(params_pl)
        @test mat_pl isa PowerLawMaterial
        @test mat_pl.jc == 2e9
        @test mat_pl.n == 20

        # Field-dependent
        params_fd = Dict(:material => Dict(
            :type => :field_dependent, :jc => 1e10, :b0 => 0.2, :n_exponent => 25
        ))
        mat_fd = material_from_params(params_fd)
        @test mat_fd isa FieldDependentMaterial

        # Unknown type
        params_bad = Dict(:material => Dict(:type => :unknown))
        @test_throws ArgumentError material_from_params(params_bad)
    end

end
