## Element type
struct IIDMVoltMeasure <: AbstractNodeLabel
  V::Complex    # p.u.
  S::Complex       # MW + im*MVar
end
