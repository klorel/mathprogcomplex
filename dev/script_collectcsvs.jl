using CSV, DataFrames, DataArrays

function parse_commandline()
  s = ArgParseSettings()
  s.description = "Merge all atomic csv into one large csv."
  @add_arg_table s begin
    "run_dir"
    help = "Directory which contains each run folder."
    arg_type = String
    required = true
  end
  return parse_args(s)
end

function main(args)
    println("--- Merge_csv.jl ---")

    ## Setting working directory
    input_params = parse_commandline()
    workdir = input_params["run_dir"]

    info("Working in $workdir")

    run_folders = Set(filter(x->isdir(joinpath(workdir, x)), readdir(workdir)))
    println("Found $(length(run_folders)) folders")

    global_csv = DataFrames.DataFrame()

    n_addedinstances = 0
    for folder in run_folders
        curfolder = joinpath(workdir, folder)
        println("-- working on $curfolder")


        ## Collecting information
        csvfile = joinpath(curfolder, "knitro_log.csv")


        # Has timed out ?
        if isfile(csvfile) != 0
            hastimedout = false
        end

        ## Reading Knitro csv
        if !hastimedout
            cur_csv = CSV.read(csvfile, delim=";")
        end
        cur_csv[:slv_hasTimedOut] = hastimedout
        !hastimedout || info("hastimedout...")


        if size(global_csv) == (0,0)
            global_csv = deepcopy(cur_csv)
        elseif names(global_csv) != names(cur_csv)
            @show names(cur_csv)
            @show length(names(cur_csv))
            error("Ignoring  $(params["instance"])  $(params["knname"]) :   Names")
        elseif DataFrames.eltypes(global_csv) != DataFrames.eltypes(cur_csv)
            @show DataFrames.eltypes(global_csv)
            @show DataFrames.eltypes(cur_csv)
            @show names cur_csv
    	    error("Ignoring  $(params["instance"])  $(params["knname"]) :   Elttypes")
        else
            n_addedinstances += 1
            append!(global_csv, cur_csv)
        end

        println("   Size global CSV : $(size(global_csv))")
    end

    println(global_csv)

    csvname = "report_log.csv"
    println("Writing full csv file $(joinpath(workdir, csvname))")

    CSV.write(joinpath(workdir, csvname), global_csv, delim=';')

    nbinter = size(global_csv[global_csv[:slv_hasTimedOut].==true, :], 1)
    info("$nbinter / $(size(global_csv, 1)) runs were interrupted.")

    info("$n_addedinstances / $(length(run_folders)) added to final csv.")

    global_csv
end

main(ARGS)