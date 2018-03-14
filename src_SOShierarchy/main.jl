ROOT = pwd()
include(joinpath(ROOT, "src_PowSysMod", "PowSysMod_body.jl"))
include(joinpath("func_definitions.jl"))

function main()

    ########################################
    # Construction du problème type
    OPFpbs = load_OPFproblems(MatpowerInput, joinpath("data_Matpower", "matpower", "WB2.m"))
    problem = build_globalpb!(OPFpbs)

    println("WB2 problem built")
    # print(problem)

    ########################################
    # Setting relaxation order by constraint
    d_i = Dict("cstr1"=>1, "cstr2"=>1)
    k_i = get_ki(problem)
    check_di_ki!(d_i, k_i)

    ########################################
    # Construction du sparsity pattern
    sparsity_pattern = compute_sparsitypattern(problem, d_i, k_i)

    # Extension chordale et détection des cliques maximales
    compute_chordalextension!(sparsity_pattern)
    max_cliques = compute_maxcliques(sparsity_pattern)

    ########################################
    # Relaxation degree par clique and variables par constrainte
    varsbycstr = compute_varsbycstr(problem)
    cliquevarsbycstr = compute_varsbycstr(sparsity_pattern, max_cliques, varsbycstr)
    orderbyclique = compute_cliqueorders(sparsity_pattern, d_i, varsbycstr, max_cliques)

    ########################################
    # Calcul des matrices B_i et pose du probleme
    B_i = compute_Bibycstr(problem, max_cliques, cliquevarsbycstr, orderbyclique)
    SDP_SOS = build_SDP_SOS(problem, d_i, max_cliques, B_i, cliquevarsbycstr, orderbyclique)

    ########################################
    # Calcul d'une solution par un solveur
    m = make_JuMPproblem(SDP_SOS)

end

main()
