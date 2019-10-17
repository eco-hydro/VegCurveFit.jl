# Write your own tests here.
using nlminb

begin
    start = [ 1.0 ]; # [1, 2, 3, 4]
    f(x) = x[1] - cos(x[1]) 
    r = optim_nlminb(start, f, verbose = false)
    r
    # par  = start
    # npar = length(par);
    # iv   = zeros(Int32, 78 + 3 * npar);
    # v    = zeros(Float64, 130 + Int((npar * (npar + 27)) / 2));

    # libpath2 = "/mnt/e/Research/julia/nlminb.jl/deps/nlminb.so"
    # # # Init iv and v, .Call("port_ivset", 2, iv, v)    
    # ccall((:Rf_divset, libpath2), Cvoid, 
    #     (Int32, Ptr{Int32}, Int32, Int32, Ptr{Float64}), 
    #     2, iv, length(iv), length(v), v);
end

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
upper = [ 4.0];

# function f2(x::Array{Float64, 1})::Float64
#     x[1] - cos(x[1]) 
# end
