using Test
using AxisArrays
using Random
using WaveSpectra
using WaveRealizations

@test isdefined(WaveRealizations, :_complex_amplitudes)
@test isdefined(WaveRealizations, :ComplexAmplitudes)
@test !isdefined(WaveRealizations, :complex_amplitudes)

let
    amplitudes = WaveRealizations._complex_amplitudes(MersenneTwister(1), [1.0, 4.0], :phase)

    @test size(amplitudes) == (2,)
    @test abs.(amplitudes) ≈ sqrt.([2.0, 8.0])
end

let
    amplitudes = ComplexAmplitudes(ones(2, 1), [1, 2] .* Hz, [90] .* °)

    @test WaveRealizations.istemporal(amplitudes)
    @test amplitudes.axis1 == [1, 2] .* Hz
    @test amplitudes.axis2 == [90] .* °
    @test amplitudes.coordinates == :polar
    @test amplitudes.axestypes == (:frequency, :direction)
    @test amplitudes.axesnames == (:frequency, :direction)
    @test WaveSpectra.coordinates(amplitudes) == :polar
    @test WaveSpectra.ispolar(amplitudes)
    @test !WaveSpectra.iscartesian(amplitudes)

    compact = sprint(show, amplitudes)
    @test compact == "2×1 ComplexAmplitudes{1}{Hz}{°}"

    plain = sprint(show, MIME"text/plain"(), amplitudes)
    @test occursin("with polar coordinates and axes:", plain)
    @test occursin(":frequency: [1 Hz, 2 Hz]", plain)
    @test occursin(":direction: [90°]", plain)
    @test occursin("and data:", plain)

    selection = amplitudes[2, 1]
    @test selection isa ComplexAmplitudes
    @test size(selection) == (1, 1)
    @test axisvalues(selection) == ([2Hz], [90°])

    selection = amplitudes[frequency = 1Hz .. 2Hz, direction = 90°]
    @test selection isa ComplexAmplitudes
    @test size(selection) == (2, 1)
    @test selection.data == ones(2, 1)

    @test amplitudes[:] isa ComplexAmplitudes
    @test size(amplitudes[:]) == size(amplitudes)
    @test size(amplitudes[frequency = 2]) == (1, 1)
    @test axisvalues(amplitudes[frequency = (1 + 1e-9) * Hz])[1] == [1Hz]
    @test_throws BoundsError amplitudes[frequency = 1.2Hz]
end

let
    amplitudes = ComplexAmplitudes(ones(1, 2), [0] .* rad, [2π, 4π] .* (rad / m))

    @test WaveRealizations.isspatial(amplitudes)
    @test amplitudes.axis1 == [2π, 4π] .* (rad / m)
    @test amplitudes.axis2 == [0] .* rad
    @test size(amplitudes) == (2, 1)
end

let
    kx = [-1, 1] .* (rad / m)
    ky = [-2, 0, 2] .* (rad / m)
    amplitudes = ComplexAmplitudes(ones(2, 3), kx, ky)

    @test amplitudes.coordinates == :cartesian
    @test amplitudes.axestypes == (:angular_wavenumber, :angular_wavenumber)
    @test amplitudes.axesnames == (:angular_wavenumber_1, :angular_wavenumber_2)
    @test WaveSpectra.iscartesian(amplitudes)

    selection = amplitudes[angular_wavenumber_1 = -1rad / m,
        angular_wavenumber_2 = 0rad / m .. 2rad / m]
    @test selection isa ComplexAmplitudes
    @test size(selection) == (1, 2)
    @test axisvalues(selection) == ([-1] .* (rad / m), [0, 2] .* (rad / m))

    plain = sprint(show, MIME"text/plain"(), amplitudes)
    @test occursin(":angular_wavenumber_1 (angular_wavenumber):", plain)
    @test occursin(":angular_wavenumber_2 (angular_wavenumber):", plain)
end

let
    amplitudes = ComplexAmplitudes(ones(7, 1), collect(1:7) .* Hz, [0] .* °)
    plain = sprint(show, MIME"text/plain"(), amplitudes)
    @test occursin(":frequency: [1 Hz, 2 Hz, 3 Hz, …, 5 Hz, 6 Hz, 7 Hz]", plain)
end

let
    amplitudes = ComplexAmplitudes(ones(1, 36), [1] .* Hz, 0°:10°:350°)
    plain = sprint(show, MIME"text/plain"(), amplitudes)
    @test occursin(":direction: [0°, 10°, 20°, …, 330°, 340°, 350°]", plain)
end

@test_throws ArgumentError ComplexAmplitudes(ones(2, 1), [1, 2] .* °, [1] .* °)
@test_throws ArgumentError ComplexAmplitudes(
    ones(2, 2), [1, 2] .* Hz, [1, 2] .* (rad / m))
@test_throws ArgumentError ComplexAmplitudes(
    ones(2, 2), [1, 2] .* Hz, [1, 2] .* (rad / s))

let
    f = [1, 2] .* Hz
    spectrum = OmnidirectionalSpectrum([1, 2] .* (m^2 / Hz), f)

    amplitudes = ComplexAmplitudes(MersenneTwister(1), spectrum; direction = 30°)

    @test amplitudes isa ComplexAmplitudes
    @test size(amplitudes) == (2, 1)
    @test WaveRealizations.istemporal(amplitudes)
    @test !WaveRealizations.isspatial(amplitudes)
    @test amplitudes.axis1 == f
    @test amplitudes.axis2 == [30] .* °
    @test amplitudes.axestypes == (:frequency, :direction)
end

let
    f = [1, 2] .* Hz
    spectrum = OmnidirectionalSpectrum([1, 2] .* (m^2 / Hz), f)
    amplitudes = ComplexAmplitudes(MersenneTwister(1), spectrum, :phase, 45°)

    @test amplitudes.axis2 == [45] .* °
end

let
    k = [1, 2] .* (rad / m)
    θ = [0, 90] .* °
    spectrum = Spectrum(ones(2, 2) .* (m^2 / ((rad / m) * °)), k, θ)

    amplitudes = ComplexAmplitudes(MersenneTwister(1), spectrum)

    @test amplitudes isa ComplexAmplitudes
    @test size(amplitudes) == (2, 2)
    @test WaveRealizations.isspatial(amplitudes)
    @test !WaveRealizations.istemporal(amplitudes)
    @test amplitudes.axis1 == k
    @test amplitudes.axis2 == θ
    @test amplitudes.axestypes == (:angular_wavenumber, :direction)
end

let
    kx = [-1, 1] .* (rad / m)
    ky = [-2, 0, 2] .* (rad / m)
    spectrum = Spectrum(ones(2, 3) .* m^4 / rad^2, kx, ky)

    amplitudes = ComplexAmplitudes(MersenneTwister(1), spectrum)

    @test amplitudes.coordinates == :cartesian
    @test amplitudes.axis1 == kx
    @test amplitudes.axis2 == ky
end
