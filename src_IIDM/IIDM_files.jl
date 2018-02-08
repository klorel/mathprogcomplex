##file to include in ComplexModeler.jl which contains all files associated to GOC

include("read_IIDM.jl")
include("read_IIDMeod.jl")

include("IIDMGenerator.jl")
include("IIDMLine_π.jl")
include("IIDMLoad.jl")
include("IIDMVolt.jl")

include("IIDMVoltMeasure.jl")
include("IIDMLine_πMeasure.jl")
