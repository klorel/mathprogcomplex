using DataStructures

function main()

    #############################################################
    # Parameters
    nmax_order1 = 3
    nmax_order2 = 3

    date = String(Dates.format(now(), "mm_dd-HHhMM"))
    workdir = joinpath(pwd(), "Mosek_runs", "pararuns", date)

    ispath(workdir) && rm(workdir, recursive=true)
    mkpath(workdir)


    params = OrderedSet()
    tempfolders = Set()

    matpower_path = joinpath("..", "data", "data_Matpower", "matpower")
    instances = OrderedSet(sort(readdir(matpower_path), by=x->parse(first(matchall(r"\d+", x)))))

    instances_order1 = filter(x->parse(first(matchall(r"\d+", x))) ≤ nmax_order1, instances)
    instances_order2 = filter(x->parse(first(matchall(r"\d+", x))) ≤ nmax_order2, instances)

    for instance in instances_order1
        tempfolder = randstring(4)
        while tempfolder in tempfolders
            tempfolder = randstring(4)
        end
        push!(tempfolders, tempfolder)
        push!(params, (instance[1:end-2], 1, joinpath(workdir, tempfolder)))
    end

    for instance in instances_order2
        tempfolder = randstring(4)
        while tempfolder in tempfolders
            tempfolder = randstring(4)
        end
        push!(tempfolders, tempfolder)
        push!(params, (instance[1:end-2], 2, joinpath(workdir, tempfolder)))
    end

    for param in params
        println(param)
    end

    return
    for param in params
        instancename, d, logpath = param
        run(`qsub -e $logpath/hierarchy_$(instancename).e -o $logpath/hierarchy_$(instancename).o julia dev/script_dat_to_hierarchysimple.jl $instancename $d $logpath`)
    end

end

main()