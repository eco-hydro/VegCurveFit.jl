DATE_ORIGIN = Date("2000-01-01")

FUNCS_FITDL = Dict(
  "AG"     => FitDL_AG,
  "Zhang"  => FitDL_Zhang,
  "Beck"   => FitDL_Beck,
  "Elmore" => FitDL_Elmore,
  "Gu"     => FitDL_Gu,
  "Klos"   => FitDL_Klos
)

FUNCS_doubleLog = Dict(
  "AG"     => doubleLog_AG,
  "Zhang"  => doubleLog_Zhang,
  "Beck"   => doubleLog_Beck,
  "Elmore" => doubleLog_Elmore,
  "Gu"     => doubleLog_Gu,
  "Klos"   => doubleLog_Klos
)
