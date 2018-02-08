"""
  Profiling for all matpower instances
    To be adapted...
"""
include(joinpath(ROOT,"..","src_PowSysMod", "PowSysMod_body.jl"))

struct Measurement
  nb_buses
  time
  bytes
  gctime
  nbpoolallocs
  nbmalloc
  summarysize
end

## Input
files_m = [split(name, ".")[1] for name in readdir(joinpath("instances", "matpower"))]
files_dat = [split(name, ".")[1] for name in readdir(joinpath("instances", "matpower_QCQP"))]

instances = intersect(files_m, files_dat)

nb_repeat = 3

instancename = sort(collect(instances))[1]

data_OPFpbs = Dict{String, Measurement}()
data_Problems = Dict{String, Measurement}()
data_errors = Dict{String, Float64}()

instnb = 1
for instancename in sort(collect(instances))
  println("\n-- Working on $instancename ($instnb/$(length(instances)))")
  nb_buses = parse(matchall(r"\d+", instancename)[1])

  ## Building OPF_problems
  t = bytes = gctime = nbpoolallocs = nbmalloc = sumsize = 0
  OPFpbs = build_OPFproblems(MatpowerSimpleInput, [joinpath("instances", "matpower", "$instancename.m")])
  for i=1:nb_repeat
    print("$i ")
    OPFpbs, t_, bytes_, gctime_, gc_ = @timed build_OPFproblems(MatpowerSimpleInput, [joinpath("instances", "matpower", "$instancename.m")])
    t += t_/nb_repeat
    bytes += bytes_/nb_repeat
    gctime += gctime_/nb_repeat
    nbpoolallocs += gc_.poolalloc/nb_repeat
    nbmalloc += gc_.malloc/nb_repeat
    sumsize += Base.summarysize(OPFpbs)/nb_repeat
  end
  data_OPFpbs[instancename] = Measurement(nb_buses, t, bytes, gctime, nbpoolallocs, nbmalloc, sumsize)

  for (busname, busdata) in OPFpbs["BaseCase"].ds.bus
    if haskey(busdata, "Gen_reference")
      OPFpbs["BaseCase"].mc.node_formulations[busname]["Gen_reference"] = :None
    end
  end

  ## Building OPF problem
  t = bytes = gctime = nbpoolallocs = nbmalloc = sumsize = 0
  pb = build_Problem!(OPFpbs, "BaseCase")
  for i=1:nb_repeat
    print("$i ")
    pb, t_, bytes_, gctime_, gc_ = @timed build_Problem!(OPFpbs, "BaseCase")
    t += t_/nb_repeat
    bytes += bytes_/nb_repeat
    gctime += gctime_/nb_repeat
    nbpoolallocs += gc_.poolalloc/nb_repeat
    nbmalloc += gc_.malloc/nb_repeat
    sumsize += Base.summarysize(pb)/nb_repeat
  end
  data_Problems[instancename] = Measurement(nb_buses, t, bytes, gctime, nbpoolallocs, nbmalloc, sumsize)

  export_to_dat(pb, "my$instancename.dat")
  filename_dat = joinpath("instances", "matpower_QCQP", "$instancename.dat")
  error, _ = compare_dat("my$instancename.dat", filename_dat)
  rm("my$instancename.dat")
  data_errors[instancename] = error

  instnb+=1
end

# output
filename = "buildOPFpbs_$(Dates.format(now(), "yy_u_dd_HH_MM")).csv"
touch(filename)

open(filename, "w") do outfile
  print(outfile, "instancename;nb_buses;time;bytes;gctime;nbpoolallocs;nbmalloc;summarysize;dat_error\n")

  insancenames_ord = sort(collect(keys(data_OPFpbs)), by=x->data_OPFpbs[x].nb_buses)

  for instancename in insancenames_ord
    measure = data_OPFpbs[instancename]
    print(outfile, "$instancename;$(measure.nb_buses);$(measure.time);$(measure.bytes);$(measure.gctime);$(measure.nbpoolallocs);$(measure.nbmalloc);$(measure.summarysize);$(data_errors[instancename])\n")
  end
end

filename = "build_Problem_$(Dates.format(now(), "yy_u_dd_HH_MM")).csv"
touch(filename)

open(filename, "w") do outfile
  print(outfile, "instancename;nb_buses;time;bytes;gctime;nbpoolallocs;nbmalloc;summarysize;dat_error\n")

  insancenames_ord = sort(collect(keys(data_Problems)), by=x->data_Problems[x].nb_buses)

  for instancename in insancenames_ord
    measure = data_Problems[instancename]
    print(outfile, "$instancename;$(measure.nb_buses);$(measure.time);$(measure.bytes);$(measure.gctime);$(measure.nbpoolallocs);$(measure.nbmalloc);$(measure.summarysize);$(data_errors[instancename])\n")
  end
end
