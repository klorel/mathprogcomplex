push!(LOAD_PATH,"D:\\repo\\complex-modeler\\src_PowSysMod")
using Documenter, PowSysMod

makedocs(
    #options
    format = :html,
    sitename = "Documentation PowSysMod",
    pages = [
        "Main functions" => "main_functions.md"
        "General structures" => "general_structures.md"
        "MatpowerInput" => "indexmatpower.md"
        "GOCInput" => "indexGOC.md"
        "PolynomialOptim" => "polynomialoptim_structures.md"
    ]
)
