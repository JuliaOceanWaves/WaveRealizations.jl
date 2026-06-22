using Test
using WaveRealizations

@test isdefined(WaveRealizations, :equal_energy_bins)

let
    edges, widths, centers, values = equal_energy_bins(x -> x, 2.0, 2)

    @test edges ≈ [0, √2, 2] rtol = 1e-6
    @test widths ≈ [√2, 2 - √2] rtol = 1e-6
    @test centers ≈ [√2 / 2, (2 + √2) / 2] rtol = 1e-6
    @test values ≈ centers
end

let
    edges, widths, centers, values = equal_energy_bins([0, 1, 2], [0, 1, 2], 2)

    @test eltype(values) == Float64
    @test values ≈ centers
end
