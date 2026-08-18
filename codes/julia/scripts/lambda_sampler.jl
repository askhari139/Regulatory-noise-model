using StatsBase  # for sample(::Weights)

# ============================================================
# CORE: simulate one stationary trajectory
# ============================================================

"""
    simulate_lambda_trajectory(noise_mode, noise_level, lambda_init,
                               min_lambda, max_lambda;
                               dt, t_burnin, t_sample)

Run the clamped lambda random walk that `simulate_with_noise`'s `affect!`
callback performs, and return the stationary samples (post-burnin).

The update rule per `dt` matches `affect!` exactly.
"""
function simulate_lambda_trajectory(noise_mode::String,
                                    noise_level::Float64,
                                    lambda_init::Float64,
                                    min_lambda::Float64,
                                    max_lambda::Float64;
                                    dt::Float64       = 0.01,
                                    t_burnin::Float64 = 50.0,
                                    t_sample::Float64 = 50.0)

    n_burnin = ceil(Int, t_burnin / dt)
    n_sample = ceil(Int, t_sample / dt)
    samples  = Vector{Float64}(undef, n_sample)

    λ = clamp(lambda_init, min_lambda, max_lambda)

    for step in 1:(n_burnin + n_sample)
        ξ = randn() * noise_level

        if     noise_mode == "Additive"
            λ += ξ
        elseif noise_mode == "Multiplicative"
            λ *= (1 + ξ)
        elseif noise_mode == "MultiplicativeInvLambda"
            # Mirrors simulate_with_noise's affect! branch: noise_level here is
            # already sigma*max_lambda (via effective_noise_level below) for
            # activation ranges (max_lambda > 1) -- dividing by the CURRENT λ
            # instead of using it flat makes the kick shrink once λ has already
            # strengthened, instead of eroding it back toward the boundary.
            # Inhibition ranges (max_lambda <= 1) fall back to plain ξ.
            λ *= (1 + (max_lambda > 1.0 ? randn() * (noise_level / λ) : ξ))
        elseif noise_mode == "Lognormal"
            λ *= exp(ξ)
        elseif noise_mode == "Fluctuating"
            λ = lambda_init + ξ
        elseif noise_mode == "Jumping"
            λ += (ξ < 0 ? -noise_level : noise_level)
        elseif noise_mode == "Extreme"
            λ += (ξ < 0 ? -1.0 : 1.0)
        else
            error("Unknown noise_mode '$noise_mode'")
        end

        λ = clamp(λ, min_lambda, max_lambda)

        if step > n_burnin
            samples[step - n_burnin] = λ
        end
    end

    return samples
end

# ============================================================
# HISTOGRAM REPRESENTATION
# ============================================================

"""
    LambdaHist

Compact stationary distribution: bin edges (linear scale) and bin counts.
Sample with `rand_from_hist(h)`.

For 200 bins, total memory ≈ 200 × 16 bytes = 3.2 KB regardless of
how many trajectories were used to build it.
"""
struct LambdaHist
    edges::Vector{Float64}     # length n_bins + 1
    counts::Vector{Int64}      # length n_bins
end

function build_histogram(samples::Vector{Float64},
                         min_lambda::Float64,
                         max_lambda::Float64;
                         n_bins::Int = 200)

    edges  = collect(range(min_lambda, max_lambda; length=n_bins + 1))
    counts = zeros(Int64, n_bins)
    width  = (max_lambda - min_lambda) / n_bins

    for s in samples
        # Clamp into range then bin
        s_clamped = clamp(s, min_lambda, max_lambda)
        idx = min(n_bins, max(1, ceil(Int, (s_clamped - min_lambda) / width)))
        counts[idx] += 1
    end

    return LambdaHist(edges, counts)
end

"""
    rand_from_hist(h::LambdaHist)

Draw one sample from the histogram. Picks a bin weighted by counts,
then samples uniformly within that bin.
"""
function rand_from_hist(h::LambdaHist)
    bin = sample(1:length(h.counts), Weights(h.counts))
    return h.edges[bin] + rand() * (h.edges[bin + 1] - h.edges[bin])
end

# ============================================================
# BUILD THE FULL CACHE
# ============================================================

"""
    build_lambda_cache(noise_mode, noise_levels, min_lambda, max_lambda;
                       dt, t_burnin, t_sample, n_traj, n_bins, init_grid)

Build a cache mapping noise level → LambdaHist (or, for Fluctuating mode,
noise level → Dict{Float64, LambdaHist} keyed by init-λ).

The grid-over-init-λ is only used for Fluctuating mode, where the
stationary distribution genuinely depends on the centre. For all other
modes the random walk forgets its starting point, so one histogram per
(mode, noise_level, inh-vs-act) suffices.
"""
function build_lambda_cache(noise_mode::String,
                            noise_levels::Vector{Float64},
                            min_lambda::Float64,
                            max_lambda::Float64;
                            dt::Float64       = 0.01,
                            t_burnin::Float64 = 50.0,
                            t_sample::Float64 = 50.0,
                            n_traj::Int       = 20,
                            n_bins::Int       = 200,
                            init_grid::Union{Nothing, Vector{Float64}} = nothing)

    effective_noise_level(σ) = max_lambda > 1.0 ? σ * max_lambda : σ

    # All noise modes use an init-λ grid: the stationary distribution of a
    # clamped random walk depends on starting position whenever the boundaries
    # are asymmetric, which is always the case here.
    grid = init_grid === nothing ?
        (min_lambda < 1 ? collect(1:100) ./ 100 : collect(1.0:100.0)) :
        init_grid

    cache = Dict{Float64, Dict{Float64, LambdaHist}}()
    for σ in noise_levels
        inner = Dict{Float64, LambdaHist}()
        for λ0 in grid
            samples = Float64[]
            for _ in 1:n_traj
                append!(samples,
                        simulate_lambda_trajectory(noise_mode, effective_noise_level(σ), λ0,
                                                   min_lambda, max_lambda;
                                                   dt=dt, t_burnin=t_burnin,
                                                   t_sample=t_sample))
            end
            inner[λ0] = build_histogram(samples, min_lambda, max_lambda;
                                        n_bins=n_bins)
        end
        cache[σ] = inner
    end
    return cache
end

# ============================================================
# SAMPLING API (used by replace_lambdas)
# ============================================================

"""
    sample_lambda(cache, noise_mode, noise_level, current_lambda)

Draw one sample from the appropriate cached distribution.

For Fluctuating mode, looks up the histogram for the init-λ closest to
`current_lambda` on the cache's grid. For other modes, ignores
`current_lambda`.
"""
function sample_lambda(cache, noise_mode::String,
                       noise_level::Float64,
                       current_lambda::Float64)

    inner = cache[noise_level]
    # Snap to nearest grid point
    grid = sort(collect(keys(inner)))
    idx  = argmin(abs.(grid .- current_lambda))
    return rand_from_hist(inner[grid[idx]])
end