# wave elevation using superposition of sinusoidal waves

"""
    WaveSurface(data, x, y, time)
    WaveSurface(amplitudes::ComplexAmplitudes, x, y, time; dispersion=periodic)
    WaveSurface(amplitudes::ComplexAmplitudes; x=nothing, y=nothing, time=nothing,
                dispersion=periodic)

Discrete wave-surface values with spatial axes `x` and `y` and temporal axis `time`.

`WaveSurface` behaves as a three-dimensional array indexed in `(x, y, time)`.
Indexing returns a `WaveSurface` with all dimensions retained. Select by integer indices,
by `x`, `y`, and `time` keywords, or by AxisArrays `..` intervals.

Convert a surface with singleton `x` and `y` axes to a TimeSeries `TimeArray` with
`TimeArray(surface)`.
"""
struct WaveSurface{
    T,
    TDAT <: AbstractArray{T, 3},
    TX <: AbstractVector,
    TY <: AbstractVector,
    TT <: AbstractVector
} <: AbstractArray{T, 3}
    data::TDAT
    x::TX
    y::TY
    time::TT

    function WaveSurface(data::TDAT,
            x::TX,
            y::TY,
            time::TT) where {
            T,
            TDAT <: AbstractArray{T, 3},
            TX <: AbstractVector,
            TY <: AbstractVector,
            TT <: AbstractVector
    }
        size(data) == (length(x), length(y), length(time)) ||
            throw(DimensionMismatch("Surface data size must match axis lengths."))
        foreach(xi -> _check_coordinate(xi, :space), x)
        foreach(yi -> _check_coordinate(yi, :space), y)
        foreach(ti -> _check_coordinate(ti, :time), time)
        return new{T, TDAT, TX, TY, TT}(data, x, y, time)
    end
end

Base.size(x::WaveSurface) = size(x.data)
Base.IndexStyle(::Type{<:WaveSurface}) = IndexCartesian()
Base.parent(x::WaveSurface) = x.data
Base.copy(x::WaveSurface) = WaveSurface(copy(x.data), copy(x.x), copy(x.y), copy(x.time))

"""
    TimeArray(surface::WaveSurface; start=DateTime(0))

Convert a `WaveSurface` at a single spatial point to a `TimeArray` of surface elevation.
The relative surface time axis is added to `start`.
"""
function TimeArray(surface::WaveSurface; start = DateTime(0))
    (length(surface.x) == 1 && length(surface.y) == 1) ||
        throw(ArgumentError("TimeArray conversion requires singleton x and y axes."))
    return _time_array(surface, start .+ float.(surface.time))
end

"""
    TimeArray(amplitudes::ComplexAmplitudes, time; start=DateTime(0), dispersion=periodic)

Evaluate `amplitudes` at a single spatial point and return the resulting surface
elevation as a `TimeArray`.

Unitful relative times are added to `start`. A vector of `Dates.TimeType` values, such
as `Date`, `DateTime`, or `Time`, is preserved as the `TimeArray` timestamp axis and
evaluated using elapsed time from its first value.
"""
function TimeArray(amplitudes::ComplexAmplitudes, time::AbstractVector;
        start = DateTime(0),
        dispersion = periodic)
    return TimeArray(WaveSurface(amplitudes; time, dispersion); start)
end

function TimeArray(amplitudes::ComplexAmplitudes, time::AbstractVector{<:TimeType};
        dispersion = periodic)
    isempty(time) && throw(ArgumentError("The time axis cannot be empty."))
    surface = WaveSurface(amplitudes; time = _elapsed_time(time), dispersion)
    return _time_array(surface, time)
end

_time_array(surface::WaveSurface, timestamps) =
    TimeArray(timestamps, vec(surface.data), [:elevation])

_elapsed_time(time::AbstractVector{<:TimeType}) =
    [value(Nanosecond(t - first(time))) * 1e-9s for t in time]

AxisArrays.axes(x::WaveSurface) = (x.x, x.y, x.time)
AxisArrays.axisvalues(x::WaveSurface) = (x.x, x.y, x.time)

# indexing
function AxisArrays.AxisArray(x::WaveSurface)
    return AxisArrays.AxisArray(
        x.data,
        AxisArrays.Axis{:x}(x.x),
        AxisArrays.Axis{:y}(x.y),
        AxisArrays.Axis{:time}(x.time)
    )
end

function WaveSurface(x::AxisArrays.AxisArray)
    x_axis, y_axis, time_axis = AxisArrays.axisvalues(x)
    return WaveSurface(x.data, x_axis, y_axis, time_axis)
end

Base.getindex(x::WaveSurface, i::Int) =
    getindex(x, Tuple(CartesianIndices(size(x))[i])...)
Base.getindex(x::WaveSurface, i::CartesianIndex) = getindex(x, Tuple(i)...)
Base.getindex(x::WaveSurface, ::Colon) = getindex(x, :, :, :)

function Base.getindex(x::WaveSurface, I...)
    selection = getindex(AxisArrays.AxisArray(x), _axis_selectors(AxisArrays.axes(x), I...)...)
    return WaveSurface(selection)
end

function Base.getindex(x::WaveSurface; kwargs...)
    selection = getindex(AxisArrays.AxisArray(x);
        _axis_selector_kwargs((:x, :y, :time), AxisArrays.axes(x), kwargs)...)
    return WaveSurface(selection)
end

function Base.show(io::IO, x::WaveSurface)
    shape = size(x)
    surface_unit = unit(eltype(x))
    surface_unit == NoUnits && (surface_unit = 1)
    io_fancy = IOContext(io, :fancy_exponent => true)
    print(io_fancy, shape[1], "×", shape[2], "×", shape[3], " WaveSurface{",
        surface_unit, "}")
end

function Base.show(io::IO, ::MIME"text/plain", x::WaveSurface)
    show(io, x)
    println(io, " with axes:")
    _show_axis(io, :x, :x, x.x)
    _show_axis(io, :y, :y, x.y)
    _show_axis(io, :time, :time, x.time)
    println(io, "and data:")
    show(IOContext(io, :limit => true), MIME"text/plain"(), x.data)
end

# discrete surfaces ηᵢ(xᵢ,yᵢ,timeᵢ)
"""
Evaluate the real surface represented by `amplitudes`.

The surface is the sum of regular sinusoidal components
`real(A * cis(kx * x + ky * y - ω * t))`. The negative sign on `ω * t` is the
package's temporal sign convention.

Evaluate the surface on every combination of `x`, `y`, and `time`, with omitted
coordinates set to zero. Return a `WaveSurface` with axes `x`, `y`, and `time`.
"""
function WaveSurface(
        amplitudes::ComplexAmplitudes;
        x = nothing,
        y = nothing,
        time = nothing,
        dispersion = periodic
)
    x, y, time = _coordinate(x, m), _coordinate(y, m), _coordinate(time, s)
    foreach(xi -> _check_coordinate(xi, :space), x)
    foreach(yi -> _check_coordinate(yi, :space), y)
    foreach(ti -> _check_coordinate(ti, :time), time)
    spatial = any(!iszero, x) || any(!iszero, y)
    temporal = any(!iszero, time)
    kx, ky, ω = _surface_components(amplitudes, dispersion, spatial, temporal)
    realization = _surface_function(amplitudes.data, kx, ky, ω)
    data = [realization(xi, yi, ti) for xi in x, yi in y, ti in time]
    return WaveSurface(data, x, y, time)
end

function WaveSurface(amplitudes::ComplexAmplitudes, x, y, time; dispersion = periodic)
    return WaveSurface(amplitudes; x, y, time, dispersion)
end

# continuous surfaces η(x,y,time)
"""
    surface_function(amplitudes::ComplexAmplitudes; dispersion=periodic)

Return a function accepting `(x, y, time)` that evaluates the real surface represented
by `amplitudes`.

The surface is the sum of regular sinusoidal components
`real(A * cis(kx * x + ky * y - ω * t))`. The negative sign on `ω * t` is the
package's temporal sign convention.
"""
function surface_function(amplitudes::ComplexAmplitudes; dispersion = periodic)
    space = isspatial(amplitudes) || dispersion !== periodic
    temporal = istemporal(amplitudes) || dispersion !== periodic
    kx, ky, ω = _surface_components(amplitudes, dispersion, space, temporal)
    realization = _surface_function(amplitudes.data, kx, ky, ω)

    return function (x, y, time)
        _check_coordinate(x, :space)
        _check_coordinate(y, :space)
        _check_coordinate(time, :time)
        # throw error if dispersion cannot convert between temporal<->spatial and
        #   value is not zero
        if !space && (!iszero(x) || !iszero(y))
            uconvert(rad / m, first(amplitudes.axis1), dispersion)
        end
        if !temporal && !iszero(time)
            uconvert(rad / s, first(amplitudes.axis1), dispersion)
        end
        return realization(x, y, time)
    end
end

# common functions
function _surface_components(amplitudes::ComplexAmplitudes, dispersion, space, temporal)
    if ispolar(amplitudes)
        return _polar_components(amplitudes, dispersion, space, temporal)
    end
    return _cartesian_components(amplitudes, dispersion, space, temporal)
end

function _polar_components(amplitudes::ComplexAmplitudes, dispersion, space, temporal)
    spectral_axis = amplitudes.axis1
    n = length(spectral_axis)
    θ = amplitudes.axis2
    k = space ? uconvert.(rad / m, spectral_axis, dispersion) : fill(0 * rad / m, n)
    ω = temporal ? uconvert.(rad / s, spectral_axis, dispersion) : fill(0 * rad / s, n)
    kx = k .* permutedims(cos.(θ))
    ky = k .* permutedims(sin.(θ))
    return kx, ky, repeat(ω, 1, length(θ))
end

function _cartesian_components(amplitudes::ComplexAmplitudes, dispersion, space, temporal)
    if isspatial(amplitudes)
        kx_axis = uconvert.(rad / m, amplitudes.axis1, dispersion)
        ky_axis = uconvert.(rad / m, amplitudes.axis2, dispersion)
        kx = repeat(kx_axis, 1, length(ky_axis))
        ky = repeat(permutedims(ky_axis), length(kx_axis), 1)
        k = sqrt.(kx .^ 2 .+ ky .^ 2)
        ω = temporal ? uconvert.(rad / s, k, dispersion) : fill(0 * rad / s, size(k))
        return kx, ky, ω
    end

    ωx_axis = uconvert.(rad / s, amplitudes.axis1, dispersion)
    ωy_axis = uconvert.(rad / s, amplitudes.axis2, dispersion)
    ωx = repeat(ωx_axis, 1, length(ωy_axis))
    ωy = repeat(permutedims(ωy_axis), length(ωx_axis), 1)
    ω = sqrt.(ωx .^ 2 .+ ωy .^ 2)
    if !space
        return fill(0 * rad / m, size(ω)), fill(0 * rad / m, size(ω)), ω
    end

    k = uconvert.(rad / m, ω, dispersion)
    scale = map((ki, ωi) -> iszero(ωi) ? zero(ki / oneunit(ωi)) : ki / ωi, k, ω)
    return scale .* ωx, scale .* ωy, ω
end

function _surface_function(amplitudes, kx, ky, ω)
    return function (x, y, time)
        value = zero(real(amplitudes[1]))
        @inbounds for i in eachindex(amplitudes)
            phase = kx[i] * x + ky[i] * y - ω[i] * time
            value += real(amplitudes[i] * cis(phase))
        end
        return value
    end
end

function _coordinate(x, unit)
    isnothing(x) && return [0 * unit]
    return x isa AbstractVector ? x : [x]
end

function _check_coordinate(x, domain::Symbol)
    valid = domain == :space ? isspatial(x) : istemporal(x)
    valid || throw(DimensionMismatch(
        domain == :space ? "x and y must have length dimensions." :
        "time must have time dimensions."))
end
