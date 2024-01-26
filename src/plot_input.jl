"""
    get_date_lims(dates; by=nothing)

# Examples
```
x_lims, x_ticks = date_ticks(dates)
```
"""
function date_ticks(date_min, date_max; by=nothing, expand=Dates.Month(2))
  date_lims = [date_min, date_max]
  if isa(expand, Dates.Month)
    expand = expand * [1, 1]
  end
  date_min, date_max = date_lims
  year_min, year_max = year.(date_lims)
  if year_max - year_min + 1 <= 3
    by = null_default(by, Month(3))

    x_ticks = date_min:by:date_max
    x_lims = (date_min, date_max)
  else
    by = null_default(by, Year(1))

    date_min = Date(year_min, 1, 1)
    date_max = Date(year_max + 1, 1, 1)
    x_ticks = date_min:by:date_max
    x_lims = (date_min, date_max)
  end
  xticks = (x_ticks, Dates.format.(x_ticks, "yyyy-mm-dd"))

  x_lims = x_lims[1] - expand[1], x_lims[2] + expand[2]
  x_lims, xticks
end

date_ticks(dates::AbstractVector; kw...) = date_ticks(extrema(dates)...; kw...)


function plot_input(t, y, QC_flag; date_lims=nothing, base_size=4.5, kw...)
  # level_names = ["snow", "cloud", "shadow", "aerosol", "marginal", "good"]
  level_names_r = ["good", "marginal", "snow", "cloud", "aerosol", "shadow"]
  # I_x, I_y = match2(level_names, level_names_r)
  # I_x = [1, 2, 3, 4, 5, 6]
  # flgs = I_x
  flgs = factor(level_names_r)
  # flgs = [6, 5, 1, 2, 4, 3]
  qc_shape = [:circle, :rect, :xcross, :dtriangle, :dtriangle, :utriangle]
  qc_colors = ["grey60", "#00BFC4", "#F8766D", "#C77CFF", "#B79F00", "#C77CFF"]
  qc_size = [0.5, 0.5, 0.5, 0, 0, 0] .+ base_size

  x_lims, x_ticks = date_lims === nothing ? date_ticks(t) : date_ticks(date_lims...)

  p = plot(t, y,
    xticks=x_ticks, xlims=x_lims,
    gridlinewidth=1, grid=:x,
    label="",
    color="black",
    framestyle=:box, kw...)

  for i = 1:6
    ind = findall(QC_flag .== flgs[i])
    # print(ind)
    # println(i, " ", length(ind))
    scatter!(p, t[ind], y[ind],
      markersize=qc_size[i],
      markerstrokewidth=1,
      markerstrokecolor=qc_colors[i],
      label=level_names_r[i],
      markercolor=qc_colors[i],
      markershape=qc_shape[i])
  end
  # p = spike_rm!(val, 0.5, x = date, QC_flag = QC_flag, p = p)
  p
end


export plot_input, date_ticks
