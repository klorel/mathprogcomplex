using Mosek

ROOT = pwd()
include(joinpath(ROOT, "src_SOShierarchy", "SOShierarchy.jl"))

## Lasserre2001 problem 1
problem, relax_ctx = lasserre_ex1()

max_cliques = get_maxcliques(relax_ctx, problem)
moments_params = build_sparsity(relax_ctx, problem, max_cliques)

mmtrel_pb = MomentRelaxationPb(relax_ctx, problem, moments_params, max_cliques)

sdpinstance = build_SDPInstance(relax_ctx, mmtrel_pb)
sdp_instance = read_SDPInstance(pwd())

sdp = SDP_Problem()

set_constraints!(sdp, sdp_instance)
set_blocks!(sdp, sdp_instance)
set_matrices!(sdp, sdp_instance)
set_linear!(sdp, sdp_instance)
set_const!(sdp, sdp_instance)

primal=SortedDict{Tuple{String,String,String}, Float64}()
dual=SortedDict{String, Float64}()

solve_mosek(sdp::SDP_Problem, primal::SortedDict{Tuple{String,String,String}, Float64}, dual::SortedDict{String, Float64}, debug=true)


@test primal[("moment_cstr", "1", "1")] ≈ 0.492635 atol=1e-4
# @test primal[("moment_cstr", "x1", "1")] ≈ 1.000000 atol=1e-4
# @test primal[("moment_cstr", "x1", "x1")] ≈ 3.073741 atol=1e-4
# @test primal[("moment_cstr", "x1^2", "1")] ≈-0.036871 atol=1e-4
# @test primal[("moment_cstr", "x1^2", "x1")] ≈ 0.000000 atol=1e-4
# @test primal[("moment_cstr", "x1^2", "x1^2")] ≈ 1.000000 atol=1e-4
# @test primal[("moment_cstr", "x1^2", "x2")] ≈-0.002228 atol=1e-4
# @test primal[("moment_cstr", "x1_*_x2", "1")] ≈-0.044236 atol=1e-4
# @test primal[("moment_cstr", "x1_*_x2", "x1")] ≈ 0.002228 atol=1e-4
# @test primal[("moment_cstr", "x1_*_x2", "x1^2")] ≈-0.000000 atol=1e-4
# @test primal[("moment_cstr", "x1_*_x2", "x1_*_x2")] ≈ 0.768593 atol=1e-4
# @test primal[("moment_cstr", "x1_*_x2", "x2")] ≈ 0.002228 atol=1e-4
# @test primal[("moment_cstr", "x2", "1")] ≈ 1.000000 atol=1e-4
# @test primal[("moment_cstr", "x2", "x1")] ≈ 1.044236 atol=1e-4
# @test primal[("moment_cstr", "x2", "x2")] ≈ 3.073741 atol=1e-4
# @test primal[("moment_cstr", "x2^2", "1")] ≈-0.036871 atol=1e-4
# @test primal[("moment_cstr", "x2^2", "x1")] ≈-0.002228 atol=1e-4
# @test primal[("moment_cstr", "x2^2", "x1^2")] ≈-0.384296 atol=1e-4
# @test primal[("moment_cstr", "x2^2 x", "1_*_x2")] ≈ 0.000000 atol=1e-4
# @test primal[("moment_cstr", "x2^2", "x2")] ≈ 0.000000 atol=1e-4
# @test primal[("moment_cstr", "x2^2", "x2^2")] ≈ 1.000000 atol=1e-4



## Lasserre2001 problem 2
problem, relax_ctx = lasserre_ex2()

max_cliques = get_maxcliques(relax_ctx, problem)
moments_params = build_sparsity(relax_ctx, problem, max_cliques)

mmtrel_pb = MomentRelaxationPb(relax_ctx, problem, moments_params, max_cliques)

sdpinstance = build_SDPInstance(relax_ctx, mmtrel_pb)
sdp_instance = read_SDPInstance(pwd())

sdp = SDP_Problem()

set_constraints!(sdp, sdp_instance)
set_blocks!(sdp, sdp_instance)
set_matrices!(sdp, sdp_instance)
set_linear!(sdp, sdp_instance)
set_const!(sdp, sdp_instance)

primal=SortedDict{Tuple{String,String,String}, Float64}()
dual=SortedDict{String, Float64}()

solve_mosek(sdp::SDP_Problem, primal::SortedDict{Tuple{String,String,String}, Float64}, dual::SortedDict{String, Float64}, debug=true)


@test primal[("moment_cstr", "1", "1")] ≈ 11.4580 atol=1e-4