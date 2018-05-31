## GLOBAL OPTIMIZATION WITH POLYNOMIALS AND THE PROBLEM OF MOMENTS
# JEAN B. LASSERRE, 2001
# http://www.ii.uib.no/~lennart/drgrad/Lasserre2001.pdf

@testset "Lasserre2001 small real problems" begin

    @testset "Example 5 (p. 811)" begin
        problem = Problem()
        x1 = Variable("x1", Real); add_variable!(problem, x1)
        x2 = Variable("x2", Real); add_variable!(problem, x2)

        set_objective!(problem, -(x1-1)^2 -(x1-x2)^2 -(x2-3)^2)
        add_constraint!(problem, "crt1", (1-(x1-1)^2) >> 0)
        add_constraint!(problem, "crt2", (1-(x1-x2)^2) >> 0)
        add_constraint!(problem, "crt3", (1-(x2-3)^2) >> 0)

        logpath = joinpath("Mosek_runs", "testfolder")
        ispath(logpath) && rm(logpath, recursive=true); mkpath(logpath)
        ## Order 1
        relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                            d = 1)

        primobj, dualobj =run_hierarchy(problem, relax_ctx, logpath)
        @test primobj ≈ -3 atol=1e-6
        @test dualobj ≈ -3 atol=1e-6

        ## Order 2
        relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                            d = 2)

        primobj, dualobj =run_hierarchy(problem, relax_ctx, logpath)
        @test primobj ≈ -2 atol=1e-6
        @test dualobj ≈ -2 atol=1e-6

        rm(logpath, recursive = true)
    end
end

