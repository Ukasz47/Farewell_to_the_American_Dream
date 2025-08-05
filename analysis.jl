
function analyze_wealth_distribution(h, W_stat)
    @assert length(h.W) == length(W_stat) "Grid and distribution must be of same length"

    cum_dist = cumsum(W_stat)
    total_wealth = dot(h.W, W_stat)

    mean_wealth = total_wealth
    median_index = findfirst(x -> x ≥ 0.5, cum_dist)
    median_wealth = h.W[median_index]

    function group_stats(p)
        cutoff = 1.0 - p
        idx = findfirst(x -> x ≥ cutoff, cum_dist)
        group_mass = sum(W_stat[idx:end])
        group_wealth = dot(h.W[idx:end], W_stat[idx:end])
        share = group_wealth / total_wealth
        avg_wealth = group_wealth / group_mass
        return (share, avg_wealth)
    end

    function gini_index(W, dist)
        cum_wealth = cumsum(W .* dist)
        B = sum(cum_wealth .* dist) / total_wealth
        G = 1 - 2 * B
        return G
    end

    return (
        mean_wealth = mean_wealth,
        median_wealth = median_wealth,
        gini = gini_index(h.W, W_stat),
        top50 = group_stats(0.50),
        top10 = group_stats(0.10),
        top1  = group_stats(0.01),
        top01 = group_stats(0.001)
    )
end

function portfolio_stats_by_wealth_decile(h, Wz_dist, φ)
    N, M = size(h.W, 1), size(h.z, 1)
    @assert size(Wz_dist) == (N, M)
    @assert size(φ) == (N, M)

    W_dist_marg = vec(sum(Wz_dist, dims=2))
    cum_dist = cumsum(W_dist_marg)

    avg_r_by_n = [sum(((1 .- φ[n, :]) .* (h.R - 1) .+ φ[n, :] .* h.E_d) .* Wz_dist[n, :]) / sum(Wz_dist[n, :]) for n in 1:N]
    avg_φ_by_n = [sum(φ[n, :] .* Wz_dist[n, :]) / sum(Wz_dist[n, :]) for n in 1:N]

    cutoffs = collect(0.0:0.1:1.0)
    results = Dict()

    for i in 1:10
        p_low = cutoffs[i]
        p_high = cutoffs[i+1]
        idx_low = findfirst(x -> x ≥ p_low, cum_dist)
        idx_high = findfirst(x -> x ≥ p_high, cum_dist)

        idx_low = isnothing(idx_low) ? N : idx_low
        idx_high = isnothing(idx_high) ? N : idx_high

        weights = W_dist_marg[idx_low:idx_high]
        avg_r = sum(avg_r_by_n[idx_low:idx_high] .* weights) / sum(weights)
        avg_φ = sum(avg_φ_by_n[idx_low:idx_high] .* weights) / sum(weights)

        key = Symbol("p$(Int(100 * p_low + 1))_to_p$(Int(100 * p_high))")
        results[key] = (avg_r, avg_φ)
    end

    idx_low = findfirst(x -> x ≥ 0.5, cum_dist)
    idx_high = findfirst(x -> x ≥ 0.9, cum_dist)

    idx_low = isnothing(idx_low) ? N : idx_low
    idx_high = isnothing(idx_high) ? N : idx_high

    weights = W_dist_marg[idx_low:idx_high]
    avg_r = sum(avg_r_by_n[idx_low:idx_high] .* weights) / sum(weights)
    avg_φ = sum(avg_φ_by_n[idx_low:idx_high] .* weights) / sum(weights)

    results[:p50_to_p90] = (avg_r, avg_φ)

    function stats_top_percentile(p)
        cutoff = 1.0 - p
        idx = findfirst(x -> x ≥ cutoff, cum_dist)
        idx = isnothing(idx) ? N : idx
        weights = W_dist_marg[idx:end]
        avg_r = sum(avg_r_by_n[idx:end] .* weights) / sum(weights)
        avg_φ = sum(avg_φ_by_n[idx:end] .* weights) / sum(weights)
        return (avg_r, avg_φ)
    end

    results[:top_1] = stats_top_percentile(0.01)
    results[:top_01] = stats_top_percentile(0.001)

    total_phi = sum(φ .* Wz_dist) / sum(Wz_dist)
    results[:average_phi] = total_phi

    return results
end

function simulate_agent(W_init, c_pol, φ_pol, h; T::Int = 250)
    N, M = size(φ)
    z_stat = stat_dist(h.π)
    z = Vector{Int}(undef, T)
    d = Vector{Int}(undef, T)
    W = Vector{Float64}(undef, T)

    z[1] = rand(Categorical(z_stat))
    d[1] = rand(Categorical(h.w_d))
    for t = 2:T
        z[t] = rand(Categorical(h.π[z[t-1],:]))
        d[t] = rand(Categorical(h.w_d))
    end

    c_itp = [extrapolate(interpolate((h.W,), c_pol[:, j], Gridded(Linear())), Line()) for j in 1:M]
    φ_itp = [extrapolate(interpolate((h.W,), φ_pol[:, j], Gridded(Linear())), Line()) for j in 1:M]

    W[1] = W_init
    for t = 1:T-1
        c = c_itp[z[t]](W[t])
        c = clamp(c, 1e-8, W[t] + h.z[z[t]] - 0.01)
        φ = φ_itp[z[t]](W[t])

        RET = (1 .- φ) .* h.R .+ φ .* (1 .+ h.d_nodes[d[t]])
        s = W[t] + h.w * h.z[z[t]] - c
        W[t+1] = s * RET
    end

    return W
end

function quantile_weighted(x::Vector{Float64}, w::Vector{Float64}, p::Float64)
    idx = sortperm(x)
    x_sorted = x[idx]
    w_sorted = w[idx]
    cum_w = cumsum(w_sorted)
    total = cum_w[end]
    threshold = p * total
    for i in eachindex(cum_w)
        if cum_w[i] ≥ threshold
            return x_sorted[i]
        end
    end
    return x_sorted[end]
end

function mobility_simulation(W_stat, c_pol, φ_pol, h; T::Int = 50, reps::Int = 100_000)

    start_labels = ["Median", "Mean", "P90", "P99", "P99_9"]
    starts = [
        quantile_weighted(h.W, W_stat, 0.5),
        dot(W_stat, h.W),
        quantile_weighted(h.W, W_stat, 0.90),
        quantile_weighted(h.W, W_stat, 0.99),
        quantile_weighted(h.W, W_stat, 0.999)
    ]

    top_0_1_cutoff = quantile_weighted(h.W, W_stat, 0.999)
    top_1_cutoff   = quantile_weighted(h.W, W_stat, 0.99)
    top_10_cutoff  = quantile_weighted(h.W, W_stat, 0.90)
    top_50_cutoff  = quantile_weighted(h.W, W_stat, 0.50)

    results = Dict()

    @showprogress 1 "Simulating agents" for (i, W_init) in enumerate(starts)
        final_wealth = [simulate_agent(W_init, c_pol, φ_pol, h, T = T)[end] for _ in 1:reps]

        results[Symbol(start_labels[i])] = (
            p_top01 = mean(w -> w ≥ top_0_1_cutoff, final_wealth),
            p_top1  = mean(w -> w ≥ top_1_cutoff, final_wealth),
            p_top10 = mean(w -> w ≥ top_10_cutoff, final_wealth),
            p_top50 = mean(w -> w ≥ top_50_cutoff, final_wealth)
        )
    end

    return results, starts
end


