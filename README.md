# Nvib.jl
[![DOI](https://zenodo.org/badge/72933503.svg)](https://zenodo.org/badge/latestdoi/72933503)

Julia source code for the chapter titled “Nonlocal continuum modeling of nanostructures” by Ekrem TUFEKCI and Serhan Aydın AYA for the upcoming book "Experimental characterization, predictive mechanical and thermal modeling of nanostructures and their polymer composite”  to be published by Elsevier in early 2017.

__Example__

This function is created for a parametric analysis in mind, however we can calculate nondimensional
frequencies for specific conditions:

If we want to determine the nondimensional frequency for the first and second mode of a nonlocal
beam (in-plane vibration analysis) having slenderness ratio of Λ=150, opening angle θT=120 and small
scale parameter R0/γ=1 where γ=1.56 nm we will follow the given steps:

```julia
julia> Pkg.clone("git@github.com:serhanaya/Nvib.jl")
julia> analysis1 = IPBeam()  # In-plane analysis
julia> nondimfrq = freqTab(analysis1, min_lambda=150, max_lambda=150, no_lambda=1, min_rad0=1.56,
    max_rad0=1.56, no_rad0=1, min_thetat=2pi/3, max_thetat=2pi/3, no_thetat=1, min_step=0,
    max_step=10e10, size_step=2, tolerance=1e-6, iter_max=10000, no_of_roots=2)
julia> nondimfrq[1, 1, 1, :]  # gives an array consist of the nondimensional frequencies for the
                              # first and second mode.
```
