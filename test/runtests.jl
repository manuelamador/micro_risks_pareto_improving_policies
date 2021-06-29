# -*- coding: utf-8 -*-
using Test
using Aiyagari 
using Roots

Ay = Aiyagari

const atol = 1e-3

@testset verbose=true "households structs" begin

    @testset "CRRA and Log" begin 
        u = CRRA(2.0)
        @test u(1.0) ≈ -1 

        u = Log()
        @test u(2.0) ≈ log(2.0) 

        @test Ay.outside(u, 1.0, 2.0, 0.8) ≈ log(1.0) + 0.8 * 2.0
        @test Ay.inside(u, 0.8) ≈ 0.8
        @test Ay.inverse_inside(u, 0.8) ≈ 0.8
        @test Ay.ss_value(u, 0.7, 0.8) ≈ log(0.7) / (1 - 0.8)
    end 

    @testset "EZ" begin 
        u = EZ(ies = 1/2, ra = 2.0)
        @test Ay.inside(u, 1.0) ≈ 1.0
        @test Ay.inverse_inside(u, Ay.inside(u, 2.0)) ≈ 2.0 

        # with EZ, v = 0 is like minus infinite 
        # -> should delivered zero value no matter c
        @test Ay.outside(u, 100.0, 0.0, 0.9) ≈ 0.0

        u = EZ(ies = 1, ra = 1)
        @test Ay.inside(u, 1.0) ≈ 0.0
        @test Ay.inverse_inside(u, Ay.inside(u, 2.0)) ≈ 2.0 

        # with EZ, v = 0 is like minus infinite 
        # -> should delivered zero value no matter c
        @test Ay.outside(u, 100.0, 0.0, 0.9) ≈ 0.0
    end 

    @testset "GHH" begin
        # GHH Preferences
        g = GHH(θ = 1.0, ν = 1.0)
        @test Ay.n_of_w(g, 1.0) ≈ 1.0
        @test Ay.disutility(g, 1.0) ≈ 1/2
        @test Ay.disutility_given_w(g, 1.0) ≈ 0.5
        @test Ay.labor_income_given_w(g, 1.0) ≈ 1.0 
        @test Ay.disutility(g, Ay.n_of_w(g, 10.0)) ≈ Ay.disutility_given_w(g, 10.0) 
        @test Ay.labor_income_given_w(g, 10.0) ≈ 10.0 * Ay.n_of_w(g, 10.0)
    end 

    @testset "Fixed Labor" begin 
        # FixedLabor Preferences
        g = FixedLabor(n = 2.0)
        @test Ay.n_of_w(g, 1.0) ≈ 2.0
        @test Ay.disutility(g, 1.0) ≈ 0.0
        @test Ay.disutility_given_w(g, 1.0) ≈ 0.0
        @test Ay.labor_income_given_w(g, 2.0) ≈ 4.0 
        @test Ay.disutility(g, Ay.n_of_w(g, 10.0)) ≈ 0.0 
        @test Ay.labor_income_given_w(g, 10.0) ≈ 10.0 * Ay.n_of_w(g, 10.0)
    end

end 

@testset verbose=true "household solutions" begin

    @testset "CRRA" begin
        # Standard household 
        h = Household(u=CRRA(2.0), v = FixedLabor(n = 1.0))

        @test Ay.get_c(h, 1.0, 1.0) ≈ 2.0 
        @test Ay.get_x(h, 1.0, 1.0) ≈ 2.0 

        sol = solve_stationary_household(h, 0.02, 1.0, tol = 1e-6)

        @test isapprox(sol.v[1], -25.54764599280174; atol)
        @test isapprox(sol.v[end], -15.109944737737447; atol)
        @test isapprox(stationary_pdf(h, sol)[1], 0.02333572741322308; atol)

    end 

    @testset "log" begin
        h = Household(u=Log(), v = FixedLabor())
        sol = solve_stationary_household(h, 0.02, 1.0, tol = 1e-6)

        @test isapprox(sol.v[1], -4.375219445045956; atol)
        @test isapprox(sol.v[end], 6.460514883157743; atol)
        @test isapprox(stationary_pdf(h, sol)[1], 0.05929518440648782; atol)
    end 

    @testset "EZ" begin
        # EZ preferences
        h = Household(u=EZ(ra = 2.0), v = FixedLabor())
        sol = solve_stationary_household(h, 0.02, 1.0, tol = 1e-6)

        @test isapprox(sol.v[1], 0.7828451950204168; atol)
        @test isapprox(sol.v[end], 1.3236150071700645; atol)
        @test isapprox(stationary_pdf(h, sol)[1], 0.023335727413377416; atol)

        h = Household(u = EZ(ies = 1/2, ra = 6), v = FixedLabor())
        sol = solve_stationary_household(h, 0.02, 1.0, tol = 1e-6)

        @test isapprox(sol.v[1], 0.769928033365381; atol)
        @test isapprox(sol.v[end], 1.315461685135013; atol)
        @test isapprox(stationary_pdf(h, sol)[1], 0.015128445715577437; atol)

    end 

    @testset "CRRA + GHH" begin
        # GHH versions
        h = Household(u=CRRA(2), v = GHH(θ = 1.0, ν = 0.2))

        @test Ay.get_c(h, 1.0, 1.0) ≈ 2.0 
        @test Ay.get_x(h, 1.0, 1.0) ≈ 2.0 - Ay.disutility(h.v, 1.0)


        sol = solve_stationary_household(h, 0.02, 1.0, tol = 1e-6)

        @test isapprox(sol.v[1], -32.17875754326171; atol)
        @test isapprox(sol.v[end], -17.320584417232183; atol)

        pdf = stationary_pdf(h, sol) 

        # policy iterator 
        @test all(isapprox.(Ay.iterate_pdf(pdf, h, sol.pol), pdf; atol))

        # testing the functions that check validity of household solution

        @test !Ay.amax_binding(sol)

        @test is_valid(sol)
    end 

    @testset "labor supply" begin
        # labor supply aggregator
        h = Household(u=CRRA(2), v = GHH(θ = 1.0))
        @test labor_supply(h, 1.0) ≈ 0.8386500803145994
    end 
end 


@testset "calibration" begin
    cal = calibration(10, 10, 0.9, 0.03, 0.01)

    @test cal[1][1] ≈ 0.0012309558783683773
    @test cal[2][1] ≈ 0.7874996334707172
end 


@testset verbose=true "Technology     " begin 

    @testset "CRS - CD" begin
        t = Technology()

        rK = rK_from_r(;t, r = 0.01)
        @test rK ≈ 0.11
        mpk = rK
        @test k_from_mpk(t; mpk, n = 1) ≈ 5.153725410455877
        @test Ay.rL_from_mpk(t, mpk) ≈ 1.1509986750018124
        @test get_y(t, n = 1, k = 1) ≈ 1.0 

        t = Technology(A = 2.0^(1/0.67))
        @test get_y(t, n = 1, k = 1) ≈ 2.0
    end 

    @testset "CRS - CES" begin 
        t = Technology(f = CES())
        rK = rK_from_r(;t, r = 0.01)
        mpk = rK
        @test get_y(t, n = 1, k = 1) ≈ 1.0
        @test k_from_mpk(t; mpk, n = 1) ≈ 2.09261314562519
        @test Ay.rL_from_mpk(t, mpk) ≈ 0.977983316917682

        # checking A: 
        let A = 2.0^(1/0.67) 
            t = Technology(A = 1.0, f=CES())
            t2 = Technology(A = A, f=CES())
            mpk = 0.11
            @test get_y(t, n = A * 3, k = 2) ≈ get_y(t2, n = 3, k = 2)
            @test A * Ay.rL_from_mpk(t, mpk) ≈ Ay.rL_from_mpk(t2, mpk)
            @test A * k_from_mpk(t; mpk, n = 1) ≈ k_from_mpk(t2; mpk, n = 1)
        end 
    end 

    @testset verbose=true "Euler's formula" begin 
        # #### Checking Euler's formula 

        # CobbDouglas
        @testset "CobbDouglas" begin 
            t = Technology(A = 2.0)
            mpk = 0.2 
            n = 3
            k = k_from_mpk(t; mpk, n)
            @test k * mpk + Ay.rL_from_mpk(t, mpk) * n ≈ get_y(t; k, n) 
        end 

        # CES
        @testset "CES" begin 
            t = Technology(A = 2.0, f = CES())
            mpk = 0.2 
            n = 2
            k = k_from_mpk(t; mpk, n)
            @test mpk * k  + Ay.rL_from_mpk(t, mpk) * n ≈ get_y(t; k, n) 
        end 
    end 

    @testset verbose=true "Markups" begin
        # ### Checking Technology with Markups 

        @testset "CobbDouglas" begin
            # CobbDouglas
            t = MarkupTechnology(μ = 1.2, A = 2.0, X = 0.0, m = 0.0)
            rK = 0.2
            n = 3
            mpk = mpk_from_after_tax_rK(t, rK) 
            k = k_from_mpk(t, mpk = mpk, n = n)
            w = Ay.rL_from_mpk(t, mpk)
            mpl = Ay.mpn_from_mpk(t, mpk)

            @test get_mpk(t; k, n) ≈ mpk
            @test mpl ≈ Ay.get_μ(t) * w
            @test k * mpk + mpl * n ≈ get_y(t; k, n)
            @test k * rK + w * n  + Ay.get_Π(t, (n = n, k = k, rK = rK, w =  w)) ≈ get_y(t; k, n)

        end 

        @testset "CES" begin
            # CES
            t = MarkupTechnology(μ = 1.2, A = 2.0, f = CES())
            rK = 0.2
            n = 3
            mpk = mpk_from_after_tax_rK(t, rK) 
            k = k_from_mpk(t;  mpk = mpk, n = n)
            w = Ay.rL_from_mpk(t, mpk)
            mpl = Ay.mpn_from_mpk(t, mpk)

            @test get_mpk(t; k, n) ≈ mpk
            @test mpl ≈ Ay.get_μ(t) * w
            @test k * mpk + mpl * n ≈ get_y(t; k, n)
            @test k * rK + w * n  + Ay.get_Π(t, (n = n, k = k, rK = rK, w =  w)) ≈ get_y(t; k, n)
        end 

        @testset "CES + X" begin 
            # CES + X (fixed cost)
            t = MarkupTechnology(μ = 1.2, A = 2.0, f = CES(), X = 1.0)
            rK = 0.2
            n = 3
            mpk = mpk_from_after_tax_rK(t, rK) 
            k = k_from_mpk(t;  mpk = mpk, n = n)
            w = Ay.rL_from_mpk(t, mpk)
            mpl = Ay.mpn_from_mpk(t, mpk)

            @test get_mpk(t; k, n) ≈ mpk
            @test mpl ≈ t.μ * w
            @test k * mpk + mpl * n ≈ get_y(t; k, n) + t.X
            @test k * rK + w * n  + Ay.get_Π(t, (n = n, k = k, rK = rK, w =  w)) ≈ get_y(t; k, n)
        end 

        @testset "CES + X + m (intermediates)" begin
            # CES + X (fixed cost) + m (intermediates)
            t = MarkupTechnology(μ = 1.2, A = 2.0, f = CES(), X = 1.0, m = 0.05)
            rK = 0.2
            n = 3
            mpk = mpk_from_after_tax_rK(t, rK) 
            k = k_from_mpk(t;  mpk = mpk, n = n)
            w = Ay.rL_from_mpk(t, mpk)
            mpl = Ay.mpn_from_mpk(t, mpk)

            @test get_mpk(t; k, n) ≈ mpk
            @test mpl ≈ t.μ * w
            @test k * mpk + mpl * n ≈ get_y(t; k, n) + t.X
            @test k * rK + w * n  + Ay.get_Π(t, (n = n, k = k, rK = rK, w =  w)) ≈ get_y(t; k, n)
        end
    end 
end

