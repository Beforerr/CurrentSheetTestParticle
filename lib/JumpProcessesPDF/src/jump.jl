using StatsBase, Distributions
using DimensionalData
using DataFrames, DataFramesMeta

"""
    EmpiricalDistribution(data::Vector{T} where T <: Real)

Create a discrete empirical distribution based on observations.
"""
function EmpiricalDistribution(data)
    data = sort(data) #sort the observations
    empirical_cdf = ecdf(data) #create empirical cdf
    data_clean = unique(data) #remove duplicates to avoid allunique error
    cdf_data = empirical_cdf.(data_clean) #apply ecdf to data
    pmf_data = vcat(cdf_data[1], diff(cdf_data)) #create pmf from the cdf
    return DiscreteNonParametric(data_clean, pmf_data) #define distribution
end

function rand_jump(x0, dists; incremental = false, bounds = nothing)
    dist = local_dist(x0, dists)
    x_new = incremental ? rand(dist) + x0 : rand(dist)
    return isnothing(bounds) ? x_new : reflect_into(x_new, bounds...)
end

function local_dist(x, dists)
    grid = dists.dims[1].val
    n = length(grid)
    return if x <= grid[1]
        dists[1]
    elseif x >= grid[end]
        dists[n]
    else
        # find upper index
        j = searchsortedfirst(grid, x)
        i = j - 1
        # linear position
        t = (x - grid[i]) / (grid[j] - grid[i])
        rand() < t ? dists[j] : dists[i]
    end
end

function reflect_into(x, a, b)
    return if x < a
        a + (a - x)
    elseif x > b
        b - (x - b)
    else
        x
    end
end

function rand_jumps!(x::AbstractVector, dists; n = 1, kw...)
    sizehint!(x, length(x) + n)
    for _ in 1:n
        push!(x, rand_jump(x[end], dists; kw...))
    end
    return x
end

rand_jumps(x0, args...; kw...) = rand_jumps!([x0], args...; kw...)

function ecdf_da(df, x, y)
    ecdfs = combine(
        groupby(df, x),
        y => EmpiricalDistribution;
        renamecols = false
    )
    sort!(ecdfs, x)
    return DimArray(ecdfs[!, y], (x = ecdfs[!, x],))
end

function df_rand_jumps(results, x, y; n = 100, counts = 1024, kw...)
    da = ecdf_da(results, x, y)
    # min_x, max_x = extrema(da.dims[1])
    # min_x, max_x = extrema(results[!, x])
    # xs = range(min_x, max_x, length = counts)
    xs = results[!, x]
    return stack(rand_jumps.(xs, (da,); n, kw...))
end
