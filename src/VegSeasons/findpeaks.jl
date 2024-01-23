using DataFrames: sort!, nrow


"""
    findpeaks(x::AbstractVector{T};
        nups::Int=1,
        ndowns::Int=nups,
        zerostr::Char='0',
        peakpat=nothing, 
        verbose=false,
        options...)

# Examples

# References
- Gerhard Aigner, <https://github.com/halleysfifthinc/Peaks.jl/issues/11#issuecomment-689998279>
"""
function findpeaks(x::AbstractVector{T};
  nups::Int=1,
  ndowns::Int=nups,
  zerostr::Char='0',
  peakpat=nothing, 
  verbose=false,
  options...) where {T<:Real}

  zerostr ∉ ('0', '+', '-') && error("zero must be one of `0`, `-` or `+`")

  # generate the peak pattern with no of ups and downs or use provided one
  peakpat = Regex(peakpat === nothing ? "[+]{$nups,}[0]*[-]{$ndowns,}" : peakpat)

  # transform x into a "+-+...-+-" character string
  xs = String(map(diff(x)) do e
    e < 0 && return '-'
    e > 0 && return '+'
    return zerostr
  end)
  verbose && @show(xs)

  grps = filter(x -> length(x) > 0, findall(peakpat, xs))
  # find index positions and maximum values
  peaks = map(grps) do m
    v, i = findmax(@view x[m])
    # fix extreme value positions on the plateau
    i = floor(Int, median(findall(@view(x[m]) .== v)))
    pos_beg = first(m)
    pos_end = last(m) + 1
    pos_peak = first(m) + i - 1
    h_left = v - x[pos_beg]
    h_right = v - x[pos_end]

    (; pos_beg, pos_peak, pos_end,
      val_beg=x[pos_beg], val_peak=v, val_end=x[pos_end],
      h_left, h_right,
      h_min=min(h_left, h_right),
      h_max=max(h_left, h_right))
  end

  df_peaks = DataFrame(peaks)
  length(options) > 0 && (df_peaks = filter_peaks(df_peaks; options...))

  df_peaks
end

"""
  filter_peaks(peaks, 
    minpeakheight=typemin(T), minpeakdistance::Int=1,
    A_max=zero(T), A_min=zero(T), 
    npeaks::Int=0,
    sortstr=false, 
    history=false,
    ignored...)

# Parameters

- `peaks`: Vector{NamedTuple{}}
- `history`: if true, removed reason will be returned

# Future Improvements
- 修改pos_beg, pos_end（向左还是向右合并?）

"""
function filter_peaks(df_peaks::DataFrame;
  minpeakheight=-Inf, minpeakdistance::Int=1,
  A_max=0, A_min=0,
  npeaks::Int=0,
  verbose=true,
  ignored...)

  if (verbose)
    println("minpeakheight = $minpeakheight, A_max = $A_max, A_min = $A_min")
  end
  n = nrow(df_peaks)

  removal = falses(n) # if true, then will be removed
  # 记录每个点被剔除的原因
  status = fill("", n)
  @inbounds for i = 1:n
    p = df_peaks[i, :]
    if p.val_peak < minpeakheight
      status[i] = "minpeakheight"
    elseif abs(p.h_min) < A_min
      status[i] = "A_min"
    elseif abs(p.h_max) < A_max
      status[i] = "A_max"
    end
  end
  inds_good = findall(.!removal)

  # TODO: 剔除peaks之后，可能存在beg, end变更的情况，如何考虑这种情况？
  # 2. 剔除距离较近的peaks
  if minpeakdistance > 1
    @inbounds for i in inds_good
      removal[i] && continue
      for j in inds_good
        (removal[j] || i == j) && continue
        peak_i = df_peaks.pos_peak[i]
        peak_j = df_peaks.pos_peak[j]
        dist = abs(peak_i - peak_j)
        if dist < minpeakdistance
          # 剔除一个比较小的
          ind = ifelse(peak_i < peak_j, i, j)
          removal[ind] = true
          status[ind] = "minpeakdistance"
        end
      end
    end
  end
  
  df_peaks[:, :status] = status

  if (npeaks > 0)
    df = df_peaks[.!removal, :]
    sort!(df, :val_peak, rev=true)
    npeaks = min(npeaks, nrow(df))
    df[1:npeaks, :]
  else
    df_peaks[status .!== "minpeakheight", :]
    # deleteat!(df_peaks, df_peaks.status .== "minpeakheight")
    # df_peaks
  end
end


export findpeaks, filter_peaks
