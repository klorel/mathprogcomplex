"""
    pb, init_point = import_from_dat(instancepath::String, precondcstrspath::String)

Build the polynomial optimization problem described by the `instancepath` file,
with the dispensory preconditionning descritpion `precondcstrspath` along with
the iniitial point possibly provided in the file (defaults value is null).
"""
function import_from_dat(instancepath::String, precondcstrspath::String)
    init_point = Point()
    variables = Dict{String, Variable}()
    exponents = Dict{String, Exponent}()
    pb = Problem()

    instance_str = open(instancepath)
    l = jump_comments!(instance_str)


    ## Collect and define variables
    line = matchall(r"\S+", l)
    while line[1] == "VAR_TYPE" && !eof(instance_str)
        if line[2] == "REAL"
            variables[line[3]] = Variable(line[3], Real)
        elseif line[2] == "BOOL"
            variables[line[3]] = Variable(line[3], Bool)
        elseif line[2] == "CPLX"
            variables[line[3]] = Variable(line[3], Complex)
        else
            error("import_to_dat(): Unknown variable type $(line[2]) for variable $(line[3]).")
        end

        init_point[variables[line[3]]] = parse(line[5]) + im*parse(line[6])

        l = readline(instance_str)
        line = matchall(r"\S+", l)
    end

    ## Mark where objective defintion begins
    mark(instance_str)


    ## Move forward to the monomial definition section
    while line[1] != "MONO_DEF" && !eof(instance_str)
        l = readline(instance_str)
        line = matchall(r"\S+", l)
    end

    while line[1] == "MONO_DEF" && !eof(instance_str)
        exponame = line[2]
        var = variables[line[3]]
        if !haskey(exponents, exponame)
            exponents[exponame] = Exponent()
        end
        exponents[exponame] = product(exponents[exponame], Exponent(Dict(var=>Degree(parse(line[5]), parse(line[6])))))
        l = readline(instance_str)
        line = matchall(r"\S+", l)
    end

    ## Reset stream to the objective defintion
    reset(instance_str)
    l = readline(instance_str)
    line = matchall(r"\S+", l)

    ## Build polynomial objective
    p = Polynomial()
    while line[2] == "OBJ"
        λ = parse(line[5]) + im*parse(line[6])
        var1, var2 = line[3:4]
        if line[1] == "MONO"
            p += λ * exponents[var1]
        else
            p += λ * (var1!="NONE" ? variables[var1] : 1) * (var2!="NONE" ? variables[var2] : 1)
        end
        l = readline(instance_str)
        line = matchall(r"\S+", l)
    end
    set_objective!(pb, p)

    while line[1] != "MONO_DEF" && !eof(instance_str)
        cstrname = line[2]
        lb = -Inf
        ub = +Inf
        p = Polynomial()
        while line[2] == cstrname
            λ = parse(line[5]) + im*parse(line[6])
            var1, var2 = line[3:4]
            if line[1] == "MONO"
                p += λ * exponents[var1]
            elseif line[1] ∈ Set(["CONST", "LIN", "QUAD"])
                p += λ * (var1!="NONE" ? variables[var1] : 1) * (var2!="NONE" ? variables[var2] : 1)
            elseif line[1] == "UB"
                ub = real(λ)
            elseif line[1] == "LB"
                lb = real(λ)
            else
                error("import_to_dat(): Unknown variable type $(line[1]) for constraint $(line[2]).")
            end
            l = readline(instance_str)
            line = matchall(r"\S+", l)
        end
        add_constraint!(pb, String(cstrname), lb << p << ub)
        l = readline(instance_str)
        line = matchall(r"\S+", l)
    end


    precond_str = open(precondcstrspath)

    l = jump_comments!(precond_str)
    line = matchall(r"\S+", l)
    while !eof(precond_str)
        if line[2] == "SQRT"
            pb.constraints[line[1]].precond = :sqrt
        else
            error("import_from_dat(): Unknown preconditioning $(line[2]) for constraint $(line[1]).")
        end
        l = readline(precond_str)
        line = matchall(r"\S+", l)
    end
    return pb, init_point
end


"""
    l = jump_comments!(io::IOStream)

Jump comments, i.e. lines starting by '#', in the `io` stream, and return the
first non commented line.
"""
function jump_comments!(io::IOStream)
    l = ""
    while !eof(io)
        l = readline(io)
        ismatch(r"\s*#", l) || break
    end
    return l
end
