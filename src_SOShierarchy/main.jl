ROOT = pwd()
include(joinpath(ROOT, "src_PowSysMod", "PowSysMod_body.jl"))
include(joinpath("func_definitions.jl"))

# function main()

########################################
# Construction du problème type
OPFpbs = load_OPFproblems(MatpowerInput, joinpath("data_Matpower", "matpower", "WB2.m"))
OPF_problem = build_globalpb!(OPFpbs)

println("WB2 problem built")

########################################
# Normalizing pb and setting relaxation order by constraint
problem = normalize_problem(OPF_problem)
relax_ctx = set_relaxation(problem, issparse = false, ismultiordered = false, d = 2)

########################################
# Construction du sparsity pattern
sparsity_pattern = compute_sparsitypattern(problem, relax_ctx)

# Extension chordale et détection des cliques maximales
compute_chordalextension!(sparsity_pattern)
max_cliques = compute_maxcliques(sparsity_pattern)

########################################
# Relaxation degree par clique and variables par constrainte
varsbycstr = compute_varsbycstr(problem)
cliquevarsbycstr = compute_varsbycstr(sparsity_pattern, max_cliques, varsbycstr)
orderbyclique = compute_cliqueorders(sparsity_pattern, varsbycstr, max_cliques, relax_ctx)

########################################
# Calcul des matrices B_i et pose du probleme
B_i = compute_Bibycstr(problem, max_cliques, cliquevarsbycstr, orderbyclique, relax_ctx)
SDP_SOS = build_SDP_SOS(problem, max_cliques, B_i, cliquevarsbycstr, orderbyclique, relax_ctx)

########################################
# Calcul d'une solution par un solveur
m = make_JuMPproblem(SDP_SOS)

# end

# main()
