## Testing global convergtence with right multi ordered hierarchy (https://arxiv.org/pdf/1508.02068.pdf)
    di = SortedDict{String, Int}("BaseCase_1_BALANCE-UNIT_Im_hi"  	=> 2,
                                 "BaseCase_1_BALANCE-UNIT_Re_hi"  	=> 2,
                                 "BaseCase_1_BALANCE-UNIT_Im_lo"  	=> 2,
                                 "BaseCase_1_BALANCE-UNIT_Re_lo"  	=> 2,
                                 "BaseCase_1_Volt_VOLTM_Re_hi"  	=> 2,
                                 "BaseCase_1_Volt_VOLTM_Re_lo"  	=> 2,
                                 "BaseCase_2_BALANCE-LOAD_Im_eq"  	=> 2,
                                 "BaseCase_2_BALANCE-LOAD_Re_eq"  	=> 2,
                                 "BaseCase_2_Volt_VOLTM_Re_hi"  	=> 2,
                                 "BaseCase_2_Volt_VOLTM_Re_lo"  	=> 2,
                                 "BaseCase_3_BALANCE-LOAD_Im_eq"  	=> 2,
                                 "BaseCase_3_BALANCE-LOAD_Re_eq"  	=> 2,
                                 "BaseCase_3_Volt_VOLTM_Re_hi"  	=> 2,
                                 "BaseCase_3_Volt_VOLTM_Re_lo"  	=> 2,
                                 "BaseCase_4_BALANCE-LOAD_Im_eq"  	=> 4,
                                 "BaseCase_4_BALANCE-LOAD_Re_eq"  	=> 4,
                                 "BaseCase_4_Volt_VOLTM_Re_hi"  	=> 4,
                                 "BaseCase_4_Volt_VOLTM_Re_lo"  	=> 4,
                                 "BaseCase_5_BALANCE-UNIT_Im_hi"  	=> 2,
                                 "BaseCase_5_BALANCE-UNIT_Im_lo"  	=> 2,
                                 "BaseCase_5_BALANCE-UNIT_Re_hi"  	=> 4,
                                 "BaseCase_5_BALANCE-UNIT_Re_lo"  	=> 4,
                                 "BaseCase_5_Volt_VOLTM_Re_hi"  	=> 4,
                                 "BaseCase_5_Volt_VOLTM_Re_lo"  	=> 4)
    logpath = joinpath(pwd(), "Mosek_runs", "worksdp_WB5sparsemulti")
    ispath(logpath) && rm(logpath, recursive = true)
    mkpath(logpath)

    OPFpbs = load_OPFproblems(MatpowerInput, joinpath("..", "data", "data_Matpower", "matpower", "WB5.m"))
    WB5 = build_globalpb!(OPFpbs)
    problem = pb_cplx2real(WB5)

    relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                               # symmetries=[PhaseInvariance],
                               issparse=true,
                               di = di)

    max_cliques = get_WB5cliques(relax_ctx, problem)
    primobj, dualobj = run_hierarchy(problem, relax_ctx, logpath, max_cliques=max_cliques)

    # Saving max_cliques
    open(joinpath(logpath, "maxcliques_relaxctx.txt"), "w") do fcliques
        println(fcliques, "max_cliques are:")
        println(fcliques, max_cliques)
        println(fcliques, "relaxation_ctx is:")
        println(fcliques, relax_ctx)
    end