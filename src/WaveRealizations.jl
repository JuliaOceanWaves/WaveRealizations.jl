"""
Package for generating realizations of surface elevation from a wave spectrum.

Can generate surfaces and time series.
Includes intermediate representation as discrete complex amplitudes.

This is part of the `JuliaOceanWaves` ecosystem.
The input spectra types are from `WaveSpectra.jl`.
"""
module WaveRealizations

using WaveSpectra: AbstractOmnidirectionalSpectrum, AbstractSpectrum, Spectrum, isdirection,
                   isevenlyspaced, m, periodic, rad, s, uconvert, °
using Unitful: eltype, ustrip, NoUnits
using Random: AbstractRNG, default_rng
using AxisArrays: axisvalues
using Dates: DateTime, Nanosecond, TimeType, value
using TimeSeries: TimeArray

import AxisArrays # axes # in the future, do `import AxisArrays: axes as AAaxes`
const axes = Base.axes # name conflict will be fixed by AxisArrays in the future
import WaveSpectra: axesnames, axestypes, coordinates, iscartesian, ispolar, isspatial,
                    istemporal
import Unitful: unit

function surface_gif end

export ComplexAmplitudes, equal_energy_bins, fft_surface, isspatial, istemporal,
       surface_function, surface_gif, WaveSurface

include("complex_amplitudes.jl")
include("equal_energy_bins.jl")
include("surfaces.jl")
include("surfaces_fft.jl")

end
