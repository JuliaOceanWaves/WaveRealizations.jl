"""
    equal_energy_bins(f::AbstractVector, x::AbstractVector, n::Integer)
    equal_energy_bins(f, xf::Number, n::Integer; n_points::Integer = max(4n + 1, 1025))

Divide `f(x)` into `n` bins with equal integral/area.

The function returns `(edges, widths, centers, values)`, where `edges` are the bin edges,
`widths` are the bin widths, `centers` are the midpoint locations of each bin, and `values`
are the spectrum values at the bin centers.
If `f` is a vector, the sample locations `x` must be provided.
If `f` is callable, it is discretized over `[0, xf]` with `n_points` samples.
"""
function equal_energy_bins(
        f::AbstractVector,
        x::AbstractVector,
        n::Integer
)
    (n < 1) && throw(ArgumentError("Number of bins must be positive."))

    n_points = length(f)
    (n_points < 2) &&
        throw(ArgumentError("At least 2 function values are required."))
    (length(x) != n_points) &&
        throw(DimensionMismatch(
            "Function values and sample locations must have the same length."))
    any(y -> y < zero(y), f) &&
        throw(ArgumentError("Function values must be nonnegative."))

    x0 = x[1]
    y_prev = f[1]
    dx = x[2] - x[1]
    cumulative = Vector{typeof((y_prev + y_prev) * dx / 2)}(undef, n_points)
    cumulative[1] = zero(eltype(cumulative))

    @inbounds for i in 2:n_points
        dx = x[i] - x[i - 1]
        y = f[i]
        cumulative[i] = cumulative[i - 1] + (y_prev + y) * dx / 2
        y_prev = y
    end

    total_area = cumulative[end]
    iszero(total_area) && throw(ArgumentError("Total integral must be nonzero."))

    edges = Vector{typeof((x[1] + x[2]) / 2)}(undef, n + 1)
    edges[1] = x0
    edges[end] = x[end]

    target_step = total_area / n
    j = 2
    @inbounds for k in 1:(n - 1)
        target = k * target_step
        while (j < n_points) && (cumulative[j] < target)
            j += 1
        end

        c_left = cumulative[j - 1]
        c_right = cumulative[j]
        x_left = x[j - 1]
        dx = x[j] - x_left
        edges[k + 1] = iszero(c_right - c_left) ? x_left :
                       x_left + (target - c_left) / (c_right - c_left) * dx
    end

    widths = Vector{typeof(edges[2] - edges[1])}(undef, n)
    centers = Vector{typeof((edges[1] + edges[2]) / 2)}(undef, n)
    @inbounds for i in 1:n
        widths[i] = edges[i + 1] - edges[i]
        centers[i] = (edges[i] + edges[i + 1]) / 2
    end

    values = Vector{typeof((f[1] + f[2]) / 2)}(undef, n)
    j = 2
    @inbounds for i in eachindex(centers)
        center = centers[i]
        while (j < n_points) && (x[j] < center)
            j += 1
        end

        x_left = x[j - 1]
        x_right = x[j]
        y_left = f[j - 1]
        y_right = f[j]
        values[i] = y_left + (y_right - y_left) * (center - x_left) / (x_right - x_left)
    end

    return edges, widths, centers, values
end

function equal_energy_bins(
        f,
        xf::Number,
        n::Integer;
        n_points::Integer = max(4n + 1, 1025)
)
    (n_points < 2) &&
        throw(ArgumentError("Number of discretization points must be at least 2."))
    x = range(zero(xf); stop = xf, length = n_points)
    values = f.(x)
    return equal_energy_bins(values, x, n)
end
