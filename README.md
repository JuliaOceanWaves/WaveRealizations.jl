# WaveRealizations.jl

[![Test](https://github.com/JuliaOceanWaves/WaveRealizations.jl/actions/workflows/Test.yml/badge.svg)](https://github.com/JuliaOceanWaves/WaveRealizations.jl/actions/workflows/Test.yml)
[![Docs](https://github.com/JuliaOceanWaves/WaveRealizations.jl/actions/workflows/Documentation.yml/badge.svg)](https://github.com/JuliaOceanWaves/WaveRealizations.jl/actions/workflows/Documentation.yml)
[![Coverage](https://codecov.io/gh/JuliaOceanWaves/WaveRealizations.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaOceanWaves/WaveRealizations.jl)

Generate wave realizations from wave elevation variance spectra.

## Installation

```julia
using Pkg
Pkg.add("WaveRealizations")
```

## Usage

```julia
using WaveRealizations
using WaveSpectra

spectrum = OmnidirectionalSpectrum([1.0, 0.25] .* (m^2 / Hz), [0.1, 0.2] .* Hz)
amplitudes = ComplexAmplitudes(spectrum)
```

## Development

```julia
using Pkg
Pkg.test("WaveRealizations")
```

## License

MIT. See `LICENSE`.
