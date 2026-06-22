using Test, SafeTestsets

@time @testset verbose=true "WaveRealizations.jl" begin
    @time @safetestset "Complex Amplitudes" begin
        include("test_complex_amplitudes.jl")
    end
    @time @safetestset "Equal Energy Bins" begin
        include("test_equal_energy_bins.jl")
    end
    @time @safetestset "Surfaces" begin
        include("test_surfaces.jl")
    end
    # documentation
    @time @safetestset "Doc Tests" begin
        include("test_doctest.jl")
    end
end
