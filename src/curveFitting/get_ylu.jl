get_range(y::AbstractArray{T, 1}) where T <: Real = [minimum(y), maximum(y)]

get_ylu(y) = get_range(y)

function get_ylu(y, years, w0, width, I; wmin = 0.2)
    n = length(w0)
    I_beg = first(I)
    I_end = last(I)
    I_all = max(1, I_beg - width) : min(n, I_end + width)
    
    w0_win = @view w0[I_all]
    I_all  = @view I_all[w0_win .> wmin]
    
    if (isempty(I_all)) 
        [NaN, NaN] # return ylu
    else
        y_win = @view y[I_all]
        ylu_min = minimum(y_win)
        ylu_max = maximum(y_win)
        
        year_win = @view years[I_all]
        if (width > length(I)*2/12) 
            ylu_min = aggregate(y_win, year_win, minimum) |> values |> median
            # ylu_min <- median(aggregate.data.frame(y_win, list(year = year_win), min)$x)
        end
        if (width > length(I)*7/12) 
            ylu_max <- aggregate(y_win, year_win, maximum) |> values |> median
        end
    end
end

function aggregate(x::AbstractArray{T, 1}, by::AbstractArray{T2, 1}, FUN::Function = mean) where {
    T<:Real, T2<:Real }

    grps = unique(by)
    vals = [];
    @views for grp in grps
        val = FUN(x[by .== grp])
        push!(vals, val)
    end
    Dict(grps .=> vals)
end

export get_range, get_ylu, aggregate
