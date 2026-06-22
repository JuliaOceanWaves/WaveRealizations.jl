using Test
using AxisArrays
using Dates
import Plots
using PlotUtils
using RecipesBase
using TimeSeries
using WaveSpectra
using WaveSpectra.DispersionRelations: gravitywaves_deepwater
using WaveRealizations

@test !isdefined(WaveRealizations, :surface)
@test isdefined(WaveRealizations, :surface_function)
@test isdefined(WaveRealizations, :WaveSurface)
@test !isdefined(WaveRealizations, :surface_slice)
@test isdefined(WaveRealizations, :fft_surface)
@test isdefined(WaveRealizations, :surface_gif)
@test !isdefined(WaveRealizations, :timeseries)

let
    amplitudes = ComplexAmplitudes(
        reshape([1 + 0im, 0.5 + 0.25im] .* m, :, 1),
        [1, 2] .* Hz,
        [0] .* °
    )
    realization = surface_function(amplitudes)

    @test realization isa Function
    @test realization(0m, 0m, 0s) ≈ 1.5m
    @test realization(0m, 0m, 1 / 4 * s) ≈ -0.5m
    @test_throws ArgumentError realization(1m, 0m, 0s)

    result = WaveSurface(amplitudes; time = [0, 1 / 4] .* s)
    @test result isa WaveSurface
    @test size(result) == (1, 1, 2)
    @test axisvalues(result) == ([0m], [0m], [0, 1 / 4] .* s)
    time_array = TimeArray(result)
    @test timestamp(time_array) == [DateTime(0), DateTime(0) + Millisecond(250)]
    @test values(time_array) ≈ [1.5, -0.5] .* m
    @test colnames(time_array) == [:elevation]
    start = DateTime(2020, 1, 2)
    @test timestamp(TimeArray(result; start)) ==
          [start, start + Millisecond(250)]
    @test timestamp(TimeArray(amplitudes, [0, 1 / 4] .* s; start)) ==
          [start, start + Millisecond(250)]
    @test values(TimeArray(amplitudes, [0, 1 / 4] .* s; start)) ≈
          [1.5, -0.5] .* m

    timestamps = [start, start + Millisecond(250)]
    absolute_time_array = TimeArray(amplitudes, timestamps)
    @test timestamp(absolute_time_array) == timestamps
    @test values(absolute_time_array) ≈ [1.5, -0.5] .* m

    dates = [Date(2020, 1, 1), Date(2020, 1, 2)]
    @test timestamp(TimeArray(amplitudes, dates)) == dates
    @test_throws ArgumentError TimeArray(amplitudes, DateTime[])
    @test result[1, 1, :].data ≈ reshape([1.5, -0.5] .* m, 1, 1, :)
    @test only(WaveSurface(amplitudes)[1, 1, 1].data) ≈ 1.5m
    @test_throws ArgumentError WaveSurface(amplitudes; x = 1m)
    @test WaveSurface(amplitudes; time = Any[0s, 1 / 4 * s])[1, 1, :].data ≈
          reshape([1.5, -0.5] .* m, 1, 1, :)

    compact = sprint(show, result)
    @test compact == "1×1×2 WaveSurface{m}"
    plain = sprint(show, MIME"text/plain"(), result)
    @test occursin(":x: [0 m]", plain)
    @test occursin(":time: [0.0 s, 0.25 s]", plain)

    recipe = only(RecipesBase.apply_recipe(Dict{Symbol, Any}(), result))
    @test recipe.plotattributes[:seriestype] == :path
    @test recipe.plotattributes[:xlabel] == "time (s)"
    @test recipe.plotattributes[:ylabel] == "surface elevation (m)"
    @test recipe.plotattributes[:title] == "x = 0 m, y = 0 m"
    @test recipe.args[1] == [0.0, 0.25]
    @test recipe.args[2] ≈ [1.5, -0.5]

    selection = result[x = 0m, y = 0m, time = 0s .. 1 / 4 * s]
    @test selection isa WaveSurface
    @test size(selection) == (1, 1, 2)
    @test axisvalues(selection) == ([0m], [0m], [0, 1 / 4] .* s)
    @test result[:] isa WaveSurface
    @test size(result[:]) == size(result)
    @test size(result[x = 1, y = 1, time = 2]) == (1, 1, 1)
    @test size(result[0m .. 0m, :, 0s .. 1 / 4 * s]) == (1, 1, 2)

    approximate = WaveSurface(
        reshape(1:6, 3, 1, 2) .* m,
        [1, 2, 3] .* m,
        [0] .* m,
        [0, 1] .* s
    )[x = (1 + 1e-9) * m]
    @test approximate.x == [1m]
    @test size(approximate) == (1, 1, 2)
    @test_throws BoundsError WaveSurface(
        reshape(1:6, 3, 1, 2) .* m,
        [1, 2, 3] .* m,
        [0] .* m,
        [0, 1] .* s
    )[x = 1.2m]

    scalar = WaveSurface(amplitudes)
    @test_throws ArgumentError RecipesBase.apply_recipe(Dict{Symbol, Any}(), scalar)
end

let
    result = WaveSurface(zeros(2, 2, 1) .* m, [0, 1] .* m, [0, 1] .* m, [0] .* s)

    @test_throws ArgumentError TimeArray(result)
    recipe = only(RecipesBase.apply_recipe(Dict{Symbol, Any}(), result))
    @test recipe.plotattributes[:seriescolor] == cgrad(:blues; rev = true)
    @test recipe.plotattributes[:xlabel] == "x (m)"
    @test recipe.plotattributes[:ylabel] == "y (m)"
    @test recipe.plotattributes[:zlabel] == "surface elevation (m)"
    @test recipe.plotattributes[:title] == "time = 0 s"

    recipe = only(RecipesBase.apply_recipe(Dict{Symbol, Any}(:seriescolor => :reds), result))
    @test recipe.plotattributes[:seriescolor] == :reds
end

let
    extension = Base.get_extension(WaveRealizations, :WaveRealizationsPlotsExt)
    result = WaveSurface(
        reshape([-2, -1, 1, 2, -1, 0, 0, 1] .* m, 2, 2, 2),
        [0, 1] .* m,
        [0, 1] .* m,
        [0, 1 / 4] .* s
    )

    @test extension._infer_fps(result) == 4
    @test extension._symmetric_limits(result) == (-2, 2)
    @test_throws ArgumentError surface_gif(WaveSurface(
        zeros(2, 2, 1) .* m, [0, 1] .* m, [0, 1] .* m, [0] .* s
    ))
    @test_throws ArgumentError extension._infer_fps(WaveSurface(
        zeros(2, 2, 3) .* m, [0, 1] .* m, [0, 1] .* m, [0, 1, 3] .* s
    ))
end

let
    amplitudes = ComplexAmplitudes(
        reshape([1 + 0im, 0.5im] .* m, :, 1),
        [1, 2] .* (rad / m),
        [0] .* °
    )
    realization = surface_function(amplitudes)

    @test realization(0m, 0m, 0s) ≈ 1m
    @test_throws ArgumentError realization(0m, 0m, 1s)

    result = WaveSurface(amplitudes, [0, π / 2] .* m, [0m], [0s])
    @test size(result) == (2, 1, 1)
    @test vec(result[:, 1, 1].data) ≈ [1, 0] .* m atol = 1e-12m
    recipe = only(RecipesBase.apply_recipe(Dict{Symbol, Any}(), result))
    @test recipe.plotattributes[:xlabel] == "x (m)"
    @test recipe.plotattributes[:ylabel] == "surface elevation (m)"
    @test recipe.plotattributes[:title] == "y = 0 m, time = 0 s"

    multiple_times = WaveSurface(
        cat(result.data, result.data; dims = 3),
        result.x,
        result.y,
        [0, 1 / 4] .* s
    )
    @test_throws ArgumentError RecipesBase.apply_recipe(Dict{Symbol, Any}(), multiple_times)
    @test_throws ArgumentError RecipesBase.apply_recipe(
        Dict{Symbol, Any}(:time => 1 / 4 * s), multiple_times)
    recipe = only(RecipesBase.apply_recipe(
        Dict{Symbol, Any}(), multiple_times[time = nextfloat(1 / 4) * s]))
    @test recipe.plotattributes[:title] == "y = 0 m, time = 0.25 s"
    @test_throws ArgumentError WaveSurface(amplitudes; time = 1s)
end

let
    kx = [-1, 0, 1] .* (rad / m)
    ky = [-1, 0, 1] .* (rad / m)
    data = zeros(ComplexF64, 3, 3) .* m
    data[3, 2] = 1m
    amplitudes = ComplexAmplitudes(data, kx, ky)

    result = WaveSurface(amplitudes; x = [0, π] .* m)
    @test vec(result[:, 1, 1].data) ≈ [1, -1] .* m
end

let
    amplitudes = ComplexAmplitudes(
        reshape([1 + 0im, 1 + 0im] .* m, 1, :),
        [1] .* (rad / m),
        [0, 90] .* °
    )
    realization = surface_function(amplitudes)

    @test realization(π * m, 0m, 0s) ≈ 0m atol = 1e-12m
    @test realization(π * m, π * m, 0s) ≈ -2m
end

let
    dispersion = gravitywaves_deepwater()
    amplitudes = ComplexAmplitudes(reshape([1 + 0im] .* m, 1, 1), [1] .* Hz, [0] .* °)
    realization = surface_function(amplitudes; dispersion)

    @test realization(0m, 0m, 0s) ≈ 1m
    @test realization(0m, 0m, 1 / 4 * s) ≈ 0m atol = 1e-12m
    quarter_wavelength = (π / 2 * rad) / uconvert(rad / m, 1Hz, dispersion)
    @test realization(quarter_wavelength, 0m, 0s) ≈ 0m atol = 1e-12m
end

let
    dispersion = gravitywaves_deepwater()
    ω = [-1, 0, 1] .* (rad / s)
    data = zeros(ComplexF64, 3, 3) .* m
    data[3, 2] = 1m
    amplitudes = ComplexAmplitudes(data, ω, ω)
    realization = surface_function(amplitudes; dispersion)

    @test realization(0m, 0m, π / 2 * s) ≈ 0m atol = 1e-12m
    @test realization((π / 2 * rad) / uconvert(rad / m, 1rad / s, dispersion), 0m, 0s) ≈
          0m atol = 1e-12m
end

let
    amplitudes = ComplexAmplitudes(
        reshape([1 + 0im, 0.5 + 0.25im] .* m, :, 1),
        [1, 2] .* Hz,
        [0] .* °
    )
    result = fft_surface(amplitudes)

    @test result isa WaveSurface
    @test size(result) == (1, 1, 6)
    @test axisvalues(result)[3] ≈ collect((0:5) .* (s / 6))
    @test only(result[1, 1, 1].data) ≈ 1.5m
end

@test_throws ArgumentError fft_surface(ComplexAmplitudes(
    ones(ComplexF64, 2, 2) .* m,
    [1, 2] .* (rad / m),
    [1, 2] .* (rad / m)
))

let
    k = [-1, 0, 1] .* (rad / m)
    data = zeros(ComplexF64, 3, 3) .* m
    data[1, 2] = 0.5m
    data[3, 2] = 0.5m
    amplitudes = ComplexAmplitudes(data, k, k)
    result = fft_surface(amplitudes)

    @test size(result) == (3, 3, 1)
    @test vec(result[:, 1, 1].data) ≈ cos.(2π .* (0:2) ./ 3) .* m
end
