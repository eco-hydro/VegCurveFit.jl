# function improved(b_new, b_old)
#     @show b_new
#     println("optimized improved:", mean(b_old.times)/mean(b_new.times), " times")
# end

first(x::AbstractArray{T, 1}) where {T<:Real} = x[1] 
last(x::AbstractArray{T, 1}) where {T<:Real} = x[end] 

day2num(x::AbstractArray{Dates.Day, 1}) = map(x -> x.value, x)

function aggregate(x::AbstractArray{T, 1}, by::AbstractArray{T2, 1}, FUN::Function = mean) where {
    T<:Real, T2<:Real }

    grps = unique(by)
    # vals = ones(length(grps));
    vals = []
    @views for grp in grps
        val = FUN(x[by .== grp])
        push!(vals, val)
    end
    Dict(grps .=> vals)
end


export first, last, day2num, aggregate
