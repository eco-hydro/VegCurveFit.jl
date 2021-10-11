# Write your own tests here.
using BenchmarkTools





begin
    using doubleLogistics
end


using doubleLogistics

## ex1
sumsq(x, y) = sum((x -y).^2)
y = repeat([1], 5)
x0 = rand(5)
nlminb(x0, sumsq, y)


## ex2
f(x) = x[1] - cos(x[1]) 
r = nlminb([ 1.0 ], f; lower = [-10], upper = [12])


## wsl version -----------------------------------------------------------------
# Dict{String,Any} with 5 entries:
#   "par"         => [-1.5703]
#   "convergence" => 0
#   "iterations"  => 16
#   "evaluations" => Dict{String,Int32}("function"=>17,"gradient"=>23)
#   "objective"   => -1.5708

## windows version -------------------------------------------------------------
# Dict{String,Any} with 5 entries:
#   "par"         => [-1.5703]
#   "convergence" => 0
#   "iterations"  => 16
#   "evaluations" => Dict{String,Int32}("function"=>17,"gradient"=>23)
#   "objective"   => -1.5708
# upper = [ 4.0];

# function f2(x::Array{Float64, 1})::Float64
#     x[1] - cos(x[1]) 
# end
