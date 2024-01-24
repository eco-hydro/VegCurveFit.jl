# function improved(b_new, b_old)
#     @show b_new
#     println("optimized improved:", mean(b_old.times)/mean(b_new.times), " times")
# end

using CategoricalArrays: CategoricalArray, CategoricalValue, levels, compress, cut

factor(args...) = CategoricalArray(args...) |> compress


null_default(x, default) = x === nothing ? default : x

day2num(x::AbstractArray{Dates.Day,1}) = map(x -> x.value, x)

function aggregate(x::AbstractVector{T}, by::AbstractArray{T2,1}, FUN::Function=mean) where {
  T<:Real,T2<:Real}

  grps = unique(by)
  # vals = ones(length(grps));
  vals = []
  @views for grp in grps
    val = FUN(x[by.==grp])
    push!(vals, val)
  end
  Dict(grps .=> vals)
end


function pmap(vec::Vector{Vector}, prop::Int)
  map(x -> x[prop], vec)
end

function pmap(vec::Vector{Dict}, prop)
  map(x -> x[prop], vec)
end

export day2num, aggregate, pmap2
