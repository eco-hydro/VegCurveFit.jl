get_range(y::AbstractArray{T, 1}) where T <: Real = [minimum(y), maximum(y)]

get_ylu(y) = get_range(y)

function get_ylu(y, years, w0, width::Integer, I; wmin = 0.2)
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
            ylu_max = aggregate(y_win, year_win, maximum) |> values |> median
        end
        [ylu_min, ylu_max]
    end
end

# update ylu_new directly
function merge_ylu!(ylu, ylu_new)
    # ylu = copy(ylu_org)
    if !isnan(ylu_new[1]); ylu_new[1] = max(ylu[1], ylu_new[1]); end
    if !isnan(ylu_new[2]); ylu_new[2] = min(ylu[2], ylu_new[2]); end
    # ylu
end


export get_range, get_ylu
