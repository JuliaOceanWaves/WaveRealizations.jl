# Create random complex amplitudes realizations from a wave spectrum

"""
    ComplexAmplitudes(data, axis1, axis2)
    ComplexAmplitudes(
        [rng::AbstractRNG,] spectrum::AbstractSpectrum,
        [axis1_spacing = nothing, axis2_spacing = nothing, randommethod = :complexamplitude]
    )
    ComplexAmplitudes(
        [rng::AbstractRNG,] spectrum::AbstractOmnidirectionalSpectrum,
        [spectral_spacing = nothing, randommethod = :complexamplitude];
        direction = 0°
    )

Complex amplitude matrix with physical-unit axes.

The first axis is a spectral variable.
The second axis is either another spectral variable (cartesian) or direction (polar).

Indexing returns a `ComplexAmplitudes` with all dimensions retained. Select by integer
indices, by axis-name keywords, or by AxisArrays `..` intervals.

The complex amplitudes can be constructed from a wave spectrum, which is assumed to be a
variance spectrum.
An omnidirectional spectrum produces a polar `ComplexAmplitudes` with one direction,
specified by `direction`.
Axis spacing may be a vector, number, or `nothing`.
When `nothing`, spacing is inferred from the corresponding evenly spaced spectrum axis.

Two randomization methods are available:

- `:complexamplitude` (default) independently samples the real and imaginary components:
  ``A = \\sqrt{V}(X + iY)``, where ``X, Y \\sim \\mathcal{N}(0, 1)``.
- `:phase` uses the spectrum-derived amplitude with a uniformly random phase:
  ``A = \\sqrt{2V}e^{i\\phi}``, where ``\\phi \\sim \\mathcal{U}[0, 2\\pi]``.

Here, ``V`` is the variance in the corresponding discrete spectral bin.
For example, for a given bin ``j`` of variance spectrum ``S(f)``, ``V=S_j Δf_j``

A random number generator can be passed as the first argument, which can be useful for
reproducible results.
"""
struct ComplexAmplitudes{
    TDAT,
    TAX1 <: AbstractVector,
    TAX2 <: AbstractVector
} <: AbstractMatrix{TDAT}
    data::Matrix{TDAT}
    axis1::TAX1
    axis2::TAX2
    coordinates::Symbol
    axestypes::Tuple{Symbol, Symbol}
    axesnames::Tuple{Symbol, Symbol}

    function ComplexAmplitudes(
            data::AbstractMatrix,
            axis1::AbstractVector,
            axis2::AbstractVector
    )
        size(data) == (length(axis1), length(axis2)) ||
            throw(DimensionMismatch(
                "Complex amplitude data size must match axis lengths."))

        if isdirection(axis1) && !isdirection(axis2)
            data, axis1, axis2 = permutedims(data), axis2, axis1
        end

        validated = Spectrum(data, axis1, axis2)
        if iscartesian(validated) && (validated.axestypes[1] != validated.axestypes[2])
            throw(ArgumentError(
                "Cartesian complex-amplitude axes must use the same spectral variable type."
            ))
        end
        return new{eltype(validated.data), typeof(validated.axis1), typeof(validated.axis2)}(
            validated.data,
            validated.axis1,
            validated.axis2,
            validated.coordinates,
            validated.axestypes,
            validated.axesnames
        )
    end
end

# default behavior / array interface
Base.size(x::ComplexAmplitudes) = size(x.data)
Base.IndexStyle(::Type{<:ComplexAmplitudes}) = IndexLinear()
Base.eltype(::Type{<:ComplexAmplitudes{TDAT}}) where {TDAT} = TDAT
function Base.copy(x::ComplexAmplitudes)
    ComplexAmplitudes(copy(x.data), copy(x.axis1), copy(x.axis2))
end

AxisArrays.axes(x::ComplexAmplitudes) = (x.axis1, x.axis2)
AxisArrays.axisvalues(x::ComplexAmplitudes) = (x.axis1, x.axis2)
coordinates(x::ComplexAmplitudes) = x.coordinates
axestypes(x::ComplexAmplitudes) = x.axestypes
axesnames(x::ComplexAmplitudes) = x.axesnames
ispolar(x::ComplexAmplitudes) = x.coordinates == :polar
iscartesian(x::ComplexAmplitudes) = x.coordinates == :cartesian
istemporal(x::ComplexAmplitudes) = istemporal(x.axis1)
isspatial(x::ComplexAmplitudes) = isspatial(x.axis1)

# indexing
function AxisArrays.AxisArray(x::ComplexAmplitudes)
    axis1 = AxisArrays.Axis{x.axesnames[1]}(x.axis1)
    axis2 = AxisArrays.Axis{x.axesnames[2]}(x.axis2)
    return AxisArrays.AxisArray(x.data, axis1, axis2)
end

function ComplexAmplitudes(x::AxisArrays.AxisArray)
    axis1, axis2 = AxisArrays.axisvalues(x)
    return ComplexAmplitudes(x.data, axis1, axis2)
end

Base.getindex(x::ComplexAmplitudes, i::Int) =
    getindex(x, Tuple(CartesianIndices(size(x))[i])...)
Base.getindex(x::ComplexAmplitudes, i::CartesianIndex) = getindex(x, Tuple(i)...)
Base.getindex(x::ComplexAmplitudes, ::Colon) = getindex(x, :, :)

function Base.getindex(x::ComplexAmplitudes, I...)
    selection = getindex(AxisArrays.AxisArray(x), _axis_selectors(AxisArrays.axes(x), I...)...)
    return ComplexAmplitudes(selection)
end

function Base.getindex(x::ComplexAmplitudes; kwargs...)
    selection = getindex(AxisArrays.AxisArray(x);
        _axis_selector_kwargs(x.axesnames, AxisArrays.axes(x), kwargs)...)
    return ComplexAmplitudes(selection)
end

function unit(x::ComplexAmplitudes, quantity::Symbol)
    ux, u1, u2 = unit(eltype(x)), unit(eltype(x.axis1)), unit(eltype(x.axis2))
    (quantity == :axis1) && return u1
    (quantity == :axis2) && return u2
    (quantity == x.axesnames[1]) && return u1
    (quantity == x.axesnames[2]) && return u2
    (quantity == :amplitude) && return ux
    throw(ArgumentError("Unknown `quantity`."))
end

@inline _axis_selector(axis, x::Integer) = x:x
@inline function _axis_selector(axis, x::Number)
    indices = findall(axis_value -> isapprox(axis_value, x), axis)
    isempty(indices) && throw(BoundsError(axis, x))
    length(indices) == 1 ||
        throw(ArgumentError("multiple axis coordinates are approximately equal to $x"))
    return only(indices):only(indices)
end
@inline _axis_selector(axis, x) = x
@inline _axis_selectors(axes, I...) =
    ntuple(index -> index <= length(I) ? _axis_selector(axes[index], I[index]) : Colon(),
        length(axes))

@inline function _axis_selector_kwargs(names, axes, kwargs)
    axis_by_name = Dict(zip(names, axes))
    return (; (key => _axis_selector(axis_by_name[key], value) for
               (key, value) in kwargs)...)
end

unit(x::ComplexAmplitudes) = unit(x, :amplitude)

# show methods
function Base.show(io::IO, x::ComplexAmplitudes)
    shape = size(x)
    amplitude_unit = unit(x)
    amplitude_unit == NoUnits && (amplitude_unit = 1)
    io_fancy = IOContext(io, :fancy_exponent => true)
    print(io_fancy, shape[1], "×", shape[2], " ComplexAmplitudes{", amplitude_unit, "}{",
        unit(x, :axis1), "}{", unit(x, :axis2), "}")
end

function Base.show(io::IO, ::MIME"text/plain", x::ComplexAmplitudes)
    show(io, x)
    println(io, " with ", coordinates(x), " coordinates and axes:")
    _show_axis(io, x.axesnames[1], x.axestypes[1], x.axis1)
    _show_axis(io, x.axesnames[2], x.axestypes[2], x.axis2)
    println(io, "and data:")
    Base.print_matrix(io, x.data)
end

function _show_axis(io, name, type, axis)
    print(io, "  :", name)
    name == type || print(io, " (", type, ")")
    print(io, ": [")
    if length(axis) <= 6
        _show_axis_values(io, axis)
    else
        _show_axis_values(io, (axis[1], axis[2], axis[3]))
        print(io, ", …, ")
        n = length(axis)
        _show_axis_values(io, (axis[n - 2], axis[n - 1], axis[n]))
    end
    println(io, "]")
end

function _show_axis_values(io, values)
    for (i, value) in enumerate(values)
        i > 1 && print(io, ", ")
        show(IOContext(io, :fancy_exponent => true), value)
    end
end

# constructors from Omnidirectional Spectrum
function ComplexAmplitudes(
        rng::AbstractRNG,
        x::AbstractOmnidirectionalSpectrum,
        spectral_spacing::Union{AbstractVector, Number, Nothing} = nothing,
        randommethod::Symbol = :complexamplitude;
        direction = 0°
)
    spectral_spacing = _axis_spacing(x.axis, spectral_spacing, "Spectral axis")
    variance = x.data .* spectral_spacing
    amplitudes = reshape(_complex_amplitudes(rng, variance, randommethod), :, 1)

    # direction vector
    direction = direction isa AbstractVector ? direction : [direction]
    length(direction) == 1 ||
        throw(DimensionMismatch("Omnidirectional spectra require exactly one direction."))
    isdirection(direction) ||
        throw(ArgumentError("Direction must have DimensionfulAngles angular units."))

    return ComplexAmplitudes(amplitudes, x.axis, direction)
end

function ComplexAmplitudes(
        x::AbstractOmnidirectionalSpectrum,
        spectral_spacing::Union{AbstractVector, Number, Nothing} = nothing,
        randommethod::Symbol = :complexamplitude;
        direction = 0°
)
    return ComplexAmplitudes(default_rng(), x, spectral_spacing, randommethod; direction)
end

function ComplexAmplitudes(
        rng::AbstractRNG,
        x::AbstractOmnidirectionalSpectrum,
        randommethod::Symbol,
        direction = 0°
)
    return ComplexAmplitudes(rng, x, nothing, randommethod; direction)
end

function ComplexAmplitudes(
        x::AbstractOmnidirectionalSpectrum,
        randommethod::Symbol,
        direction = 0°
)
    return ComplexAmplitudes(default_rng(), x, nothing, randommethod; direction)
end

# constructors from Spectrum
function ComplexAmplitudes(
        rng::AbstractRNG,
        x::AbstractSpectrum,
        axis1_spacing::Union{AbstractVector, Number, Nothing} = nothing,
        axis2_spacing::Union{AbstractVector, Number, Nothing} = nothing,
        randommethod::Symbol = :complexamplitude
)
    axis1_spacing = _axis_spacing(x.axis1, axis1_spacing, "Axis 1")
    axis2_spacing = _axis_spacing(x.axis2, axis2_spacing, "Axis 2")
    variance = x.data .* axis1_spacing .* axis2_spacing'
    amplitudes = _complex_amplitudes(rng, variance, randommethod)
    return ComplexAmplitudes(amplitudes, x.axis1, x.axis2)
end

function ComplexAmplitudes(
        x::AbstractSpectrum,
        axis1_spacing::Union{AbstractVector, Number, Nothing} = nothing,
        axis2_spacing::Union{AbstractVector, Number, Nothing} = nothing,
        randommethod::Symbol = :complexamplitude
)
    return ComplexAmplitudes(default_rng(), x, axis1_spacing, axis2_spacing, randommethod)
end

function ComplexAmplitudes(
        rng::AbstractRNG,
        x::AbstractSpectrum,
        randommethod::Symbol
)
    return ComplexAmplitudes(rng, x, nothing, nothing, randommethod)
end

function ComplexAmplitudes(
        x::AbstractSpectrum,
        randommethod::Symbol
)
    return ComplexAmplitudes(default_rng(), x, nothing, nothing, randommethod)
end

# complex amplitudes methods
function _complex_amplitudes(
        rng::AbstractRNG,
        variance::AbstractVecOrMat,
        randommethod::Symbol = :complexamplitude
)
    dims = size(variance)
    if randommethod == :complexamplitude
        a, b = randn(rng, dims...), randn(rng, dims...)
        return @. complex(a, b) * √(variance)
    elseif randommethod == :phase
        ϕ = rand(rng, dims...) * 2π
        return @. √(2variance) * cis(ϕ)
    end
    throw(ArgumentError("Unknown random method `$randommethod`."))
end

function _complex_amplitudes(
        variance::AbstractVecOrMat,
        randommethod::Symbol = :complexamplitude
)
    return _complex_amplitudes(default_rng(), variance, randommethod)
end

# axis spacing options and checks
function _axis_spacing(axis::AbstractVector, spacing::AbstractVector, axis_name::String)
    if length(spacing) != length(axis)
        msg = axis_name * " spacing length must match spectrum " * axis_name * " length."
        throw(DimensionMismatch(msg))
    end

    if unit(eltype(spacing)) != unit(eltype(axis))
        msg = axis_name * " spacing units must match spectrum " * axis_name * " units."
        throw(DimensionMismatch(msg))
    end

    return spacing
end

function _axis_spacing(axis::AbstractVector, spacing::Number, axis_name::String)
    if unit(typeof(spacing)) != unit(eltype(axis))
        msg = axis_name * " spacing units must match spectrum " * axis_name * " units."
        throw(DimensionMismatch(msg))
    end

    return fill(spacing, length(axis))
end

function _axis_spacing(axis::AbstractVector, ::Nothing, axis_name::String)
    if length(axis) < 2
        msg = "axis must have at least two elements."
        throw(ArgumentError(msg))
    end

    if !isevenlyspaced(axis)
        msg = axis_name * " spacing cannot be inferred from unevenly spaced spectrum " *
              axis_name * ". Provide spacing directly."
        throw(ArgumentError(msg))
    end

    return fill(axis[2] - axis[1], length(axis))
end
