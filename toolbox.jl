function power_grid(W_min::Float64, W_max::Float64, N::Int; ξ::Float64 = 2.0)
    base = LinRange(0.0, 1.0, N)
    return W_min .+ (W_max - W_min) .* base .^ ξ
end

function rouwenhorst(N::Int,ρ::Float64,σ::Float64)
    p = (1 + ρ) / 2
    q = p
    π = [p 1-p; 1-q q] 

    for n in 3:N
        π_new = zeros(n, n)
        π_new[1:n-1, 1:n-1] += p * π
        π_new[1:n-1, 2:n]   += (1-p) * π
        π_new[2:n, 1:n-1]   += (1-q) * π
        π_new[2:n, 2:n]     += q * π
        π_new[2:n-1, :] .= π_new[2:n-1, :] ./ 2

        π = π_new
    end

    σ_z = sqrt(σ^2 / (1 - ρ^2))
    z = collect(range(-σ_z * sqrt(N-1), stop=σ_z * sqrt(N-1), length=N))

    return π, z
end

function income_process(N::Int, ρ::Float64, σ::Float64, ρ_top::Float64, z_top::Float64)
    @assert 0 ≤ ρ_top ≤ 1 "ρ_top must be between 0 and 1"

    # 1. Rouwenhorst process
    π_r, z_r = rouwenhorst(N, ρ, σ)
    z_r = exp.(z_r)

    # 2. Extended transtition matrix
    π_ext = zeros(N + 1, N + 1)

    # 3. Scaled transitions for normal states
    μ_top = 0.05
    π_T = (1 - ρ_top) * μ_top / (1 - μ_top)

    π_ext[1:N, 1:N] .= (1 - π_T) * π_r
    π_ext[1:N, N+1] .= π_T

    # 4. Transitions from top 5%
    π_ext[N+1, N+1] = ρ_top
    π_ext[N+1, 1:N] .= (1 - ρ_top) / N

    # 5. Top income ratio
    z = vcat(z_r, z_top)

    # 6. Productivity normalization: E[z] = 1
    stat = stat_dist(π_ext)
    z_norm = z ./ dot(z, stat)

    return π_ext, z_norm
end

function risk_process(MD::Int, h)
    d_nodes = zeros(MD)
    w_d = zeros(MD)

    d_nodes[1] = -1
    w_d[1] = h.p_B

    mean_lognorm = (h.E_d + h.p_B) / (1 - h.p_B)
    x_d, w_d_n = gausshermite(MD - 1)
    d_nodes[2:end] = mean_lognorm * exp.(sqrt(2) * h.σ_d .* x_d .- 0.5 * h.σ_d^2)
    w_d[2:end] = (1 - h.p_B) * w_d_n ./ sqrt(Base.MathConstants.pi)
    return d_nodes, w_d
end

function stat_dist(P)
    n = size(P, 1)
    @assert size(P,1) == size(P,2) "P must be a square matrix"
    
    π = fill(1.0 / n, n)
    max_iter = 1_000_000
    tol = 1e-12

    for iter in 1:max_iter
        π_new = π' * P
        diff = norm(π_new - π', Inf)
        if diff < tol
            return vec(π_new)
        end
        π = vec(π_new)
    end

    @warn "Stationary distribution did not converge within $max_iter iterations. Max difference: $(norm(π' * P - π', Inf))"
    return π
end

function trans_mat(W::Vector, z::Vector, π::Matrix, W_prim::Matrix)
    N = length(W) 
    M = length(z)  
    MN = M * N  

    row_idx = Int[]
    col_idx = Int[]
    values  = Float64[]

    for i in 1:M   
        for j in 1:N   
            a′ = W_prim[j, i]  

            col = (i - 1) * N + j 

            if a′ <= W[1]
                k_low = 1
                w = 1.0
            elseif a′ >= W[end]
                k_low = N
                w = 1.0
            else
                k_low = searchsortedlast(W, a′)
                k_high = k_low + 1
                w = (W[k_high] - a′) / (W[k_high] - W[k_low])
            end

            for k in 1:M
                prob = π[i, k]

                row1 = (k - 1) * N + k_low
                push!(row_idx, row1)
                push!(col_idx, col)
                push!(values, w * prob)

                if k_low < N
                    row2 = (k - 1) * N + k_low + 1
                    push!(row_idx, row2)
                    push!(col_idx, col)
                    push!(values, (1 - w) * prob)
                end
            end
        end
    end

    return sparse(row_idx, col_idx, values, MN, MN)
end

function trans_mat_risky(W::Vector, z::Vector, π::Matrix, W_prim::Array{Float64, 3}, w_d::Vector)
    N = length(W)
    M = length(z)
    MD = length(w_d)
    MN = N * M

    row_idx = Int[]
    col_idx = Int[]
    values  = Float64[]

    for k in 1:MD
        for j in 1:M 
            for i in 1:N 

                a′ = W_prim[i, j, k]
                col = (j - 1) * N + i

                if a′ <= W[1]
                    k_low = 1
                    w = 1.0
                elseif a′ >= W[end]
                    k_low = N
                    w = 1.0
                else
                    k_low = searchsortedlast(W, a′)
                    k_high = k_low + 1
                    w = (W[k_high] - a′) / (W[k_high] - W[k_low])
                end

                for jp in 1:M
                    prob = π[j, jp] * w_d[k]

                    row1 = (jp - 1) * N + k_low
                    push!(row_idx, row1)
                    push!(col_idx, col)
                    push!(values, w * prob)

                    if k_low < N
                        row2 = (jp - 1) * N + k_low + 1
                        push!(row_idx, row2)
                        push!(col_idx, col)
                        push!(values, (1 - w) * prob)
                    end
                end
            end
        end
    end

    return sparse(row_idx, col_idx, values, MN, MN)
end

function add_fortune_destruction(T::SparseMatrixCSC{Float64, Int}, h, a_birth::Float64)
    N, M = length(h.W), length(h.z)
    MN = N * M
    z_stat = stat_dist(h.π)

    # Interpolation: a_birth between W[k_low] and W[k_low+1]
    if a_birth <= h.W[1]
        k_low = 1
        w = 1.0
    elseif a_birth >= h.W[end]
        k_low = N
        w = 1.0
    else
        k_low = searchsortedlast(h.W, a_birth)
        k_high = k_low + 1
        w = (h.W[k_high] - a_birth) / (h.W[k_high] - h.W[k_low])
    end

    # Target rows (where new agents are added)
    target_rows = Int[]
    target_weights = Float64[]

    for j in 1:M
        base = (j - 1) * N
        if 1 ≤ k_low ≤ N
            push!(target_rows, base + k_low)
            push!(target_weights, w * z_stat[j])
        end
        if k_low < N
            push!(target_rows, base + k_low + 1)
            push!(target_weights, (1 - w) * z_stat[j])
        end
    end

    birth_rows = repeat(target_rows, inner=MN)
    birth_cols = repeat(1:MN, outer=length(target_rows))
    birth_vals = repeat(target_weights, inner=MN) .* h.ξ

    T_scaled = (1 - h.ξ) * T
    B = sparse(birth_rows, birth_cols, birth_vals, MN, MN)
    return T_scaled + B
end

function get_distribution(mean_wealth, k_prim, S, h)

    W_prim = zeros(N, M, MD)
    for k in 1:MD
        W_prim[:, :, k] .= k_prim .* h.R .+ S .* (1 + h.d_nodes[k])
    end

    T = trans_mat_risky(h.W, h.z, h.π, W_prim, h.w_d)
    T = add_fortune_destruction(T, h, mean_wealth)
    Wz_stat = stat_dist(T')
    Wz_stat = reshape(Wz_stat, N, M)
    k_supply = sum(Wz_stat .* k_prim)
    S_demand = sum(Wz_stat .* S)
    W_stat = sum(Wz_stat, dims=2)[:]

    return W_stat, Wz_stat
end



