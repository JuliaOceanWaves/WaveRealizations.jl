# discrete surfaces using FFT
using FFTW: ifft, irfft
using Unitful: ustrip

"""
    fft_surface(amplitudes::ComplexAmplitudes)

Generate a realization on the natural FFT grid.

Polar amplitudes must have one direction and an evenly spaced one-sided frequency or
wavenumber axis starting at the first positive FFT frequency. Spatial polar amplitudes
must point along positive x or positive y.

Cartesian amplitudes must be spatial, use full signed FFT wavenumber grids, and already
be Hermitian. A positive-quadrant Cartesian grid is insufficient because its missing
directional quadrants cannot be inferred.
"""
function fft_surface(amplitudes::ComplexAmplitudes)
    return ispolar(amplitudes) ? _fft_surface_polar(amplitudes) :
           _fft_surface_cartesian(amplitudes)
end

function _fft_surface_polar(amplitudes::ComplexAmplitudes)
    size(amplitudes, 2) == 1 ||
        throw(ArgumentError("Polar FFT realization requires exactly one direction."))
    angular_axis = istemporal(amplitudes) ?
                   uconvert.(rad / s, amplitudes.axis1, periodic) :
                   uconvert.(rad / m, amplitudes.axis1, periodic)
    output_axis, data = _irfft_axis(amplitudes.data[:, 1], angular_axis)

    if istemporal(amplitudes)
        return WaveSurface(reshape(data, 1, 1, :), [0m], [0m], output_axis)
    end

    θ = amplitudes.axis2[1]
    if isapprox(cos(θ), one(cos(θ))) && isapprox(sin(θ), zero(sin(θ)))
        return WaveSurface(reshape(data, :, 1, 1), output_axis, [0m], [0s])
    elseif isapprox(sin(θ), one(sin(θ))) && isapprox(cos(θ), zero(cos(θ)))
        return WaveSurface(reshape(data, 1, :, 1), [0m], output_axis, [0s])
    end
    throw(ArgumentError("Spatial polar FFT direction must be positive x or positive y."))
end

function _irfft_axis(amplitudes, angular_axis)
    length(angular_axis) >= 2 ||
        throw(ArgumentError("FFT realization requires at least two amplitudes."))
    isevenlyspaced(angular_axis) ||
        throw(ArgumentError("FFT spectral axis must be evenly spaced."))
    dκ = angular_axis[2] - angular_axis[1]
    isapprox(first(angular_axis), dκ) ||
        throw(ArgumentError("One-sided FFT spectral axis must start at its spacing."))

    n_amplitudes = length(amplitudes)
    n = 2 * (n_amplitudes + 1)
    amplitude_unit = unit(eltype(amplitudes))
    coefficients = zeros(ComplexF64, n_amplitudes + 2)
    coefficients[2:(n_amplitudes + 1)] .= (n / 2) .* ustrip.(amplitude_unit, amplitudes)
    data = irfft(coefficients, n) .* amplitude_unit
    output_axis = collect((0:(n - 1)) .* (2π * rad / (n * dκ)))
    return output_axis, data
end

function _fft_surface_cartesian(amplitudes::ComplexAmplitudes)
    isspatial(amplitudes) ||
        throw(ArgumentError("Cartesian FFT realization requires spatial amplitudes."))
    kx = uconvert.(rad / m, amplitudes.axis1, periodic)
    ky = uconvert.(rad / m, amplitudes.axis2, periodic)
    px, dx = _fft_permutation(kx)
    py, dy = _fft_permutation(ky)
    ordered = amplitudes.data[px, py]
    _check_hermitian(ordered)

    nx, ny = size(ordered)
    amplitude_unit = unit(eltype(ordered))
    data = real.(ifft(ustrip.(amplitude_unit, ordered) .* (nx * ny))) .* amplitude_unit
    x = collect((0:(nx - 1)) .* (2π * rad / (nx * dx)))
    y = collect((0:(ny - 1)) .* (2π * rad / (ny * dy)))
    return WaveSurface(reshape(data, nx, ny, 1), x, y, [0s])
end

function _fft_permutation(axis)
    length(axis) >= 2 || throw(ArgumentError("FFT axes require at least two values."))
    isevenlyspaced(axis) || throw(ArgumentError("FFT axes must be evenly spaced."))
    dk = axis[2] - axis[1]
    n = length(axis)
    indices = round.(Int, ustrip.(axis ./ dk))
    required = iseven(n) ? collect((-n ÷ 2):(n ÷ 2 - 1)) : collect((-(n ÷ 2)):(n ÷ 2))
    sort(indices) == required ||
        throw(ArgumentError("Cartesian FFT axes must be full signed FFT grids."))
    return sortperm(mod.(indices, n)), dk
end

function _check_hermitian(coefficients)
    nx, ny = size(coefficients)
    @inbounds for i in 1:nx, j in 1:ny

        ii = mod1(2 - i, nx)
        jj = mod1(2 - j, ny)
        isapprox(coefficients[i, j], conj(coefficients[ii, jj])) ||
            throw(ArgumentError("Cartesian FFT amplitudes must be Hermitian."))
    end
end
