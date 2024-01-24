"""
    get_extent(w0, I_beg, I_end, nptperyear)

# Parameters
- `w0`: original initialized weights
- `opt`: 
    + maxExtendMonth
    + minExtendMonth
    + nextend
    + wmin

# Return
- `I_beg : I_end`
"""
function get_extentI(w0::AbstractVector, I_beg, I_end, nptperyear;
  maxExtendMonth=2, minExtendMonth=1, nextend=2, wmin=0.2)

  MaxExtendWidth = ceil(Int, nptperyear / 12 * maxExtendMonth)
  MinExtendWidth = ceil(Int, nptperyear / 12 * minExtendMonth)

  n = length(w0)
  I_beg2 = I_end2 = NaN
  # period <- floor(nptperyear/12*extend_month)
  if nextend === nothing
    I_beg_raw = max(1, I_beg - 1):max(1, I_beg - MaxExtendWidth)
    I_end_raw = min(n, I_end + 1):min(n, I_end + MaxExtendWidth)

    # at least `nextend` marginal point in previous and following season
    # I_beg2 and I_end2 have been constrained in the range of [1, n]
    I_beg2 = I_beg_raw[findall(w0[I_beg_raw] .> wmin)[nextend]]
    I_end2 = I_end_raw[findall(w0[I_end_raw] .> wmin)[nextend]]
  end

  # in case of previous and subsequent season good values too closed
  max_Beg = max(1, I_beg - MinExtendWidth)
  min_End = min(n, I_end + MinExtendWidth)

  I_beg2 = isnan(I_beg2) ? max_Beg : min(I_beg2, max_Beg)
  I_end2 = isnan(I_end2) ? min_End : min(I_end2, min_End)

  if isempty(I_beg2) || isempty(I_end2)
    return (I_beg:I_end)
  else
    return (I_beg2:I_end2)
  end
end


function getDateId(dates, t, direction="before")
  function sub(date)
    if direction == "before"
      # before
      ans = findall(t .<= date) |> last
      if isempty(ans)
        ans = findall(t .>= date) |> first
      end
    elseif direction == "after"
      # after
      ans = findall(t .>= date) |> first
      if isempty(ans)
        ans = findall(t .<= date) |> last
      end
    else
      ans = indexin(dates, t)
    end
    ans
  end
  map(sub, dates)
end


export get_extentI, getDateId
