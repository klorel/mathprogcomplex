ROOT = pwd()
include(joinpath(ROOT, "src_SOShierarchy", "SOShierarchy.jl"))

using Base.Profile
using ProfileView
using BenchmarkTools

function build_pb_from_dat(dat_path)
    (problem_c, point) = import_from_dat(dat_path)
    problem = pb_cplx2real(problem_c)

    return problem
end


function main()
    repo = LibGit2.GitRepo(pwd()); branch = LibGit2.shortname(LibGit2.head(repo))
    date = String(Dates.format(now(), "mm_dd-HHhMM"))
    workdir = joinpath("Mosek_runs", "benchmark", branch, date)
    ispath(workdir) && rm(workdir, recursive=true)
    mkpath(workdir)


    sols = OrderedDict("WB2"        => (2, 885.71, 905.73, true),
                    #    "WB3"        => (1, 417.25, 417.25, false),
                       "LMBM3"      => (1, 386.42, 386.42, false),
                       "WB5"        => (2, 954.82, 1146.4, true),
                       "case6ww"    => (1, 2986, 2986, false),
                       "case9"      => (2, 373.8, 1458.8, true))
                    #   "case9mod"   => (2, 234.6, 1320.4, true),
                    #   "case14"     => (2, 721.5, 5371.5, true),
                    #    "case22loop" => (1, 4538.8, 4538.8, false), ## Absent in data repo...
                    #   "case30"     => (2, 268.915, 316.49, true))


    suite = BenchmarkGroup()

    suite["order_1"] = BenchmarkGroup() #["POP_build", "pb_construction", "pb_build_SDPMosekstruct", "pb_mosek_solve"])
    suite["order_2"] = BenchmarkGroup() #["POP_build", "pb_construction", "pb_build_SDPMosekstruct", "pb_mosek_solve"])

    primal = SortedDict{Tuple{String,String,String}, Float64}()
    dual = SortedDict{Tuple{String, String, String}, Float64}()

    for d=1:2

        suite["order_$d"]["1.POP_build"] = BenchmarkGroup()
        suite["order_$d"]["2.pb_construction"] = BenchmarkGroup()
        suite["order_$d"]["3.pb_build_SDPMosekstruct"] = BenchmarkGroup()
        suite["order_$d"]["4.pb_mosek_solve"] = BenchmarkGroup()


        for (instance, (dcv, obj_rankrel, obj_opt, lackconstant)) in sols

            logpath = joinpath(workdir, "order_$d", instance)
            mkpath(logpath)

            dat_path = joinpath("..", "data", "data_Matpower", "matpower_QCQP", instance*".dat")
            suite["order_$d"]["1.POP_build"][instance] = @benchmarkable (problem = build_pb_from_dat($dat_path))
            problem = build_pb_from_dat(dat_path)

            relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                                d = d)

            suite["order_$d"]["2.pb_construction"][instance] = @benchmarkable (SOS_pb = build_relaxation($problem, $relax_ctx))
            SOS_pb = build_relaxation(problem, relax_ctx)
            export_SDP(SOS_pb, logpath)

            suite["order_$d"]["3.pb_build_SDPMosekstruct"][instance] = @benchmarkable (sdpinstance = build_mosekpb($logpath))
            sdpinstance = build_mosekpb(logpath)

            suite["order_$d"]["4.pb_mosek_solve"][instance] = @benchmarkable ((primobj, dualobj) = solve_mosek($sdpinstance, $primal, $dual; logname = joinpath($logpath, "Mosek_run.log")))
            primobj, dualobj = solve_mosek(sdpinstance, primal, dual; logname = joinpath(logpath, "Mosek_run.log"))
        end
    end

    # If a cache of tuned parameters already exists, use it, otherwise, tune and cache
    # the benchmark parameters. Reusing cached parameters is faster and more reliable
    # than re-tuning `suite` every time the file is included.
    paramspath = joinpath(dirname(@__FILE__), "params_hierarchybenchmark.json")

    if isfile(paramspath)
        loadparams!(suite, BenchmarkTools.load(paramspath)[1], :evals);
    else
        #tune!(suite)
        BenchmarkTools.save(paramspath, params(suite));
    end

    return suite
end

function showresults(results; io=STDOUT)
    data = OrderedDict()

    orders = Set()
    instances = Set()
    measurekinds = SortedSet()

    for (order, val1) in results
        for (measurekind, val2) in val1
            for (instance, t) in val2
                !haskey(data, instance) && (data[instance] = Dict())
                data[instance][(order, measurekind)] = t

                push!(orders, order)
                push!(instances, instance)
                push!(measurekinds, measurekind)
            end
        end
    end

    orders = OrderedSet(sort(collect(orders), by=x->parse(matchall(r"\d+", x)[1])))
    instances = OrderedSet(sort(collect(instances), by=x->parse(matchall(r"\d+", x)[1])))

    colofset = 10
    instanceslen = maximum(map(x->length(x), instances))
    instanceslen = max(instanceslen, length("instances"))
    ordercollen = sum([length(mkind)+colofset+3 for mkind in measurekinds])-3

    for order in orders
        print(io, order, "\n")

        print(io, "| ")
        print_string(io, "instances", instanceslen)
        print(io, "| ")
        for measurekind in measurekinds
            print_string(io, measurekind, length(measurekind)+colofset)
            print(io, "| ")
        end
        print(io, "\n")

        print(io, "| ")
        print(io, "-"^instanceslen)
        print(io, " | ")
        for measurekind in measurekinds
            print(io, "-"^(length(measurekind)+colofset), ":")
            print(io, "| ")
        end
        print(io, "\n")


        for instance in instances
            print(io, "| ")
            print_string(io, instance, instanceslen)
            print(io, "| ")
            for measurekind in measurekinds
                t = mean(results[order][measurekind][instance].times)
                print_string(io, BenchmarkTools.prettytime(t), length(measurekind)+colofset)
                print(io, "| ")
            end
            print(io, "\n")
        end
        print(io, "\n\n")
    end
end
