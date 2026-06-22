module WaveRealizationsPlotsExt

using Plots
using Unitful: NoUnits, s, unit, ustrip
using WaveRealizations: WaveSurface, isevenlyspaced

import WaveRealizations: surface_gif

"""
    plot(surface::WaveSurface)

Plot a static frame of a `WaveSurface` using Plots.jl.

At a single spatial point, plot all times as a time series. For spatial surfaces with
multiple times, select a frame by indexing the surface:

```julia
plot(surface[time = 0.2s])
```
"""
@recipe function f(surface::WaveSurface)
    x, y = surface.x, surface.y

    if length(x) == 1 && length(y) == 1
        length(surface.time) > 1 ||
            throw(ArgumentError("Cannot plot a surface with only one x, y, and time value."))
        seriestype --> :path
        xlabel --> _axis_label("time", surface.time)
        ylabel --> _axis_label("surface elevation", surface.data)
        title --> "x = $(only(x)), y = $(only(y))"
        ustrip.(surface.time), vec(ustrip.(surface.data))
    elseif length(x) == 1
        _check_single_time(surface)
        z = surface.data[:, :, 1]
        seriestype --> :path
        xlabel --> _axis_label("y", y)
        ylabel --> _axis_label("surface elevation", surface.data)
        title --> "x = $(only(x)), time = $(only(surface.time))"
        ustrip.(y), vec(ustrip.(z))
    elseif length(y) == 1
        _check_single_time(surface)
        z = surface.data[:, :, 1]
        seriestype --> :path
        xlabel --> _axis_label("x", x)
        ylabel --> _axis_label("surface elevation", surface.data)
        title --> "y = $(only(y)), time = $(only(surface.time))"
        ustrip.(x), vec(ustrip.(z))
    else
        _check_single_time(surface)
        z = surface.data[:, :, 1]
        seriestype --> :surface
        seriescolor --> cgrad(:blues; rev = true)
        xlabel --> _axis_label("x", x)
        ylabel --> _axis_label("y", y)
        zlabel --> _axis_label("surface elevation", surface.data)
        title --> "time = $(only(surface.time))"
        ustrip.(x), ustrip.(y), ustrip.(z)
    end
end

"""
    surface_gif(surface::WaveSurface, filename="surface.gif"; fps=nothing, limits=nothing,
                plot_kwargs...)

Create a GIF animation of all frames in `surface` using Plots.jl.

By default, `fps` is inferred from an evenly spaced time axis and symmetric color and
vertical limits are inferred from the largest absolute surface elevation. Pass `limits`
to override both limits. Additional keywords are forwarded to `Plots.plot`.
"""
function surface_gif(surface::WaveSurface, filename = "surface.gif";
        fps = nothing,
        limits = nothing,
        plot_kwargs...)
    length(surface.time) > 1 ||
        throw(ArgumentError("A surface animation requires at least two time frames."))
    (length(surface.x) > 1 || length(surface.y) > 1) ||
        throw(ArgumentError("A surface animation requires at least two spatial points."))

    frame_rate = isnothing(fps) ? _infer_fps(surface) : fps
    frame_limits = isnothing(limits) ? _symmetric_limits(surface) : limits

    animation = Plots.Animation()
    for time in surface.time
        frame = Plots.plot(surface[time = time];
            clims = frame_limits, zlims = frame_limits, plot_kwargs...)
        Plots.frame(animation, frame)
    end

    return Plots.gif(animation, filename; fps = frame_rate)
end

_check_single_time(surface::WaveSurface) =
    length(surface.time) == 1 ||
    throw(ArgumentError("Select one frame by indexing the surface before plotting."))

function _axis_label(name, values)
    value_unit = unit(eltype(values))
    value_unit == NoUnits && return name
    return "$name ($value_unit)"
end

function _infer_fps(surface::WaveSurface)
    isevenlyspaced(surface.time) ||
        throw(ArgumentError("Specify `fps` when the surface time axis is not evenly spaced."))
    Δt = ustrip(s, surface.time[2] - surface.time[1])
    Δt > 0 || throw(ArgumentError("Surface times must be increasing to infer `fps`."))
    return max(1, round(Int, inv(Δt)))
end

function _symmetric_limits(surface::WaveSurface)
    limit = maximum(abs, ustrip.(surface.data))
    return (-limit, limit)
end

end
