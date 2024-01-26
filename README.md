# VegCurveFit.jl: Curve fitting for Vegetation Index 

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://kongdd.github.io/nlminb.jl/stable) -->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://kongdd.github.io/nlminb.jl/dev)
[![CI](https://github.com/eco-hydro/VegCurveFit.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/eco-hydro/VegCurveFit.jl/actions/workflows/CI.yml)
[![Codecov](https://codecov.io/gh/eco-hydro/VegCurveFit.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/eco-hydro/VegCurveFit.jl/tree/master)

> Dongdong Kong

# Functions

- [x] weighted Whittaker-Henderson smoothing and interpolation
- [x] weighted Savitzky Golay filter
- [x] HANTS (untested)
- [x] Double Logistics
- [ ] The GPU version
- [ ] Growing season dividing


# References

> [1] Kong, D., McVicar, T. R., Xiao, M., Zhang, Y., Peña-Arancibia, J. L., Filippa, G., Xie, Y., Gu, X. (2022). phenofit: An R package for extracting vegetation phenology from time series remote sensing. __*Methods in Ecology and Evolution*__, 13, 1508-1527. <https://doi.org/10.1111/2041-210X.13870>
> 
> [2] Kong, D., Zhang, Y.\*, Wang, D., Chen, J., & Gu, X\*. (2020).
> Photoperiod Explains the Asynchronization Between Vegetation Carbon
> Phenology and Vegetation Greenness Phenology. 
> *Journal of Geophysical Research: Biogeosciences*, 125(8), e2020JG005636.
> <https://doi.org/10.1029/2020JG005636>
>
> [3] Kong, D., Zhang, Y.\*, Gu, X., & Wang, D. (2019). A robust method
> for reconstructing global MODIS EVI time series on the Google Earth
> Engine. __*ISPRS Journal of Photogrammetry and Remote Sensing*__, 155,
> 13–24.
>
> [4] Zhang, Q.\*, Kong, D.\*, Shi, P., Singh, V.P., Sun, P., 2018.
> Vegetation phenology on the Qinghai-Tibetan Plateau and its response
> to climate change (1982–2013). __*Agricultural and Forest Meteorology*__. 248, 408–417.
> <https://doi.org/10.1016/j.agrformet.2017.10.026>
