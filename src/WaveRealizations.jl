"""
Package for generating realizations of surface elevation from a wave spectrum.

Can generate surfaces and time series.
Includes intermediate representation as discrete complex amplitudes.

This is part of the `JuliaOceanWaves` ecosystem.
The input spectra types are from `WaveSpectra.jl`.
"""
module WaveRealizations

using SpectralSuperpositions: AbstractSuperposition, coordinates, iscartesian, isdirection,
                              isevenlyspaced, ispolar, isspatial, istemporal,
                              validate_superposition, _show_axis, (..)
using WaveSpectra: AbstractOmnidirectionalSpectrum, AbstractSpectrum, m, periodic, rad, s,
                   uconvert, °
using Unitful: eltype, ustrip, NoUnits
using Random: AbstractRNG, default_rng
using Dates: DateTime, Nanosecond, TimeType, value
using FFTW: ifft, irfft

import AxisArrays: AxisArrays, AxisArray, axisvalues # axes # in the future, do `import AxisArrays: axes as AAaxes`
const axes = Base.axes # name conflict will be fixed by AxisArrays in the future
import Unitful: unit
import SpectralSuperpositions: superposition_unit_aliases
import TimeSeries: TimeArray

function surface_gif end

export ComplexAmplitudes, WaveSurface, equal_energy_bins, fft_surface, isspatial,
       istemporal, surface_function, surface_gif, (..)

include("complex_amplitudes.jl")
include("equal_energy_bins.jl")
include("surfaces.jl")
include("surfaces_fft.jl")

end
