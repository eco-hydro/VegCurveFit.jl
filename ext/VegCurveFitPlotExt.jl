export VegCurveFitPlotExt
module VegCurveFitPlotExt


import VegCurveFit
import Plots

# :vline!, :scatter!, 
funcs = [:stroke]

for f in funcs
  @eval VegCurveFit.$f(args...; kw...) = Plots.$f(args...; kw...)
end

end
