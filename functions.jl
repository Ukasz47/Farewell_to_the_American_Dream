include("toolbox.jl")

#|============================================================================================
#|--------------Functions for finding household optimum and general equilibrium---------------
#|============================================================================================

function HH_Ayiagari_two_assets_exo(M, N, MD, ε, h)
    π = h.π
    z = h.z * h.w
    d_nodes = h.d_nodes
    w_d = h.w_d
    W = h.W

    c = h.c_init
    φ = h.φ_init
    V = h.V_init
    c_hat = similar(c)
    W_hat = similar(c)
    c_ret = similar(c)
    V_ret = similar(c)
    V = value_function_two_assets(c, φ, h, V)

    @showprogress 1 "Solving EGM" for s in 1:1000
        c_prev = copy(c)

        RET_nodes = [(1 .- φ) .* h.R .+ φ .* (1 .+ d) for d in d_nodes]

        itp_c = [extrapolate(interpolate((W,), c[:, j], Gridded(Linear())), Line()) for j in 1:M]
        itp_V = [extrapolate(interpolate((W,), V[:, j], Gridded(Linear())), Line()) for j in 1:M]
        EU_prim = zero(c)
        EV_Vpow_1mθ = zero(c)

        for m in 1:MD
            RET = RET_nodes[m]
            W_ret = (W * ones(1, M)) .* RET

            for j in 1:M
                @inbounds for i in 1:N
                    c_ij = itp_c[j](W_ret[i, j])
                    c_ret[i, j] = clamp(c_ij, 1e-8, W_ret[i, j] + z[j])
                    V_ret[i, j] = itp_V[j](W_ret[i, j])
                end
            end

            EV_Vpow_1mθ .+= w_d[m] .* V_ret.^(1 - h.θ) * π'
            U_prim = h.A .* ((W_ret .+ h.W_bar).^(-h.η))

            u_prim_m = RET .* ((1 - h.γ) .* c_ret.^(-h.ρ) .+ h.γ * U_prim) .* V_ret.^(h.ρ - h.θ)
            EU_prim .+= w_d[m] .* u_prim_m * π'
        end

        c_hat .= (h.β .* EU_prim ./ EV_Vpow_1mθ.^((h.ρ - h.θ)/(1 - h.θ))) .^(-1 / h.ρ)

        W_hat .= W .+ c_hat .- z'

        itp_c = [extrapolate(interpolate((W_hat[:, j],), c_hat[:, j], Gridded(Linear())), Line()) for j in 1:M]
        for j in 1:M
            @inbounds for i in 1:N
                c_ij = itp_c[j](W[i])
                c[i, j] = clamp(c_ij, 1e-8, W[i] + z[j])
            end
        end

        V = value_function_two_assets(c, φ, h, V)

        if maximum(abs.(c - c_prev)) < ε
            println("Converged in $s iterations")
            break
        end
    end

    k_prim = (1 .- φ) .* ((W .+ z') .- c)
    S = φ .* ((W .+ z') .- c)

    return k_prim, S, c, φ, V
end

function value_function_two_assets(c, φ, h, V)
    N, M, MD = size(c, 1), size(c, 2), length(h.d_nodes)
    W, z, π = h.W, h.z * h.w, h.π
    d_nodes, w_d = h.d_nodes, h.w_d

    V_new = similar(V)
    Vp = similar(V)

    mask = φ .> 0.0

    for iter in 1:2000
        EVp_pow = zeros(N, M)
  
        for k in eachindex(d_nodes)
            RET = (1 .- φ) .* h.R .+ φ .* (1 .+ d_nodes[k])
            Wp = RET .* ((W .+ z') .- c .- mask .* h.κ)

            itp_V = [extrapolate(interpolate((W,), V[:, j], Gridded(Linear())), Line()) for j in 1:M]
            Vp = reshape([itp_V[j](Wp[i, j]) for i in 1:N, j in 1:M], N, M)

            EVp_pow .+= w_d[k] .* (Vp .^ (1 - h.θ)) * π'
        end

        Uc= (1 - h.γ) * (1 - h.β) * c.^(1 - h.ρ)
        UV = (1 - h.γ) * h.β * EVp_pow.^((1 - h.ρ) / (1 - h.θ))
        Ub = h.γ * (1 - h.ρ) / (1 - h.η) * (1 - h.β) * h.A * ((W .+ h.W_bar).^(1 - h.η))
        V_new .= (Uc .+ UV .+ Ub).^(1 / (1 - h.ρ))

        if maximum(abs.(V_new - V)) < 1e-6
            break
        elseif iter == 2000
            println("Value Function not converged in $iter iterations")
        end
        V .= V_new
    end

    return V
end

function value_function_one_update(c, φ, h, V, if_κ)
    N, M, MD = size(c, 1), size(c, 2), length(h.d_nodes)
    W, z, π = h.W, h.z * h.w, h.π
    d_nodes, w_d = h.d_nodes, h.w_d

    V_new = similar(V)
    Vp = similar(V)

    EVp_pow = zeros(N, M)

    if if_κ == 1
        for k in eachindex(d_nodes)
            RET = (1 .- φ) .* h.R .+ φ .* (1 .+ d_nodes[k])
            Wp = RET .* ((W .+ z') .- c .- h.κ)

            itp_V = [extrapolate(interpolate((W,), V[:, j], Gridded(Linear())), Line()) for j in 1:M]
            Vp = reshape([itp_V[j](Wp[i, j]) for i in 1:N, j in 1:M], N, M)

            EVp_pow .+= w_d[k] .* (Vp .^ (1 - h.θ)) * π'
        end
    else
        Wp = h.R .* ((W .+ z') .- c)

        itp_V = [extrapolate(interpolate((W,), V[:, j], Gridded(Linear())), Line()) for j in 1:M]
        Vp = reshape([itp_V[j](Wp[i, j]) for i in 1:N, j in 1:M], N, M)

        EVp_pow = (Vp .^ (1 - h.θ)) * π'
    end

    Uc= (1 - h.γ) * (1 - h.β) * c.^(1 - h.ρ)
    UV = (1 - h.γ) * h.β * EVp_pow.^((1 - h.ρ) / (1 - h.θ))
    Ub = h.γ * (1 - h.ρ) / (1 - h.η) * (1 - h.β) * h.A * ((W .+ h.W_bar).^(1 - h.η))
    V_new .= (Uc .+ UV .+ Ub).^(1 / (1 - h.ρ))

    return V_new
end

function phi_resid(
    φ::Float64,
    i::Int,
    j::Int,
    c_interp::AbstractVector,
    V_interp::AbstractVector,
    W::Vector{Float64},
    z::Vector{Float64},
    π::Matrix,
    d_nodes::Vector{Float64},
    w_d::Vector{Float64},
    h
)
    resid = 0.0

    for m in eachindex(d_nodes)
        RET = (1 - φ) * h.R + φ * (1 + d_nodes[m])
        excess_return = 1 + d_nodes[m] - h.R
        savings = W[i] + z[j] - c_interp[j](W[i])
        savings = clamp(savings, 1e-8, W[i] + z[j])
        W_ret = savings * RET

        for jp in eachindex(z)
            c_p = c_interp[jp](W_ret)
            c_p = clamp(c_p, 1e-8, W_ret + z[jp])
            V_p = V_interp[jp](W_ret)
            wedge = V_p^(h.ρ - h.θ)
            u_p = (1 - h.γ) * c_p^(-h.ρ)
            W_p = h.γ * h.A * (W_ret .+ h.W_bar)^(-h.η)
            resid += π[j, jp] * w_d[m] * (u_p + W_p) * wedge * excess_return
        end
    end
    return resid
end

function HH_Ayiagari_two_assets_endo(M, N, MD, ε, ε_κ, h)
    π = h.π
    z = h.z * h.w
    d_nodes = h.d_nodes
    w_d = h.w_d
    W = h.W

    c = h.c_init
    φ = h.φ_init
    V = h.V_init
    c_hat = similar(c)
    W_hat = similar(c)
    c_ret = similar(c)
    V_ret = similar(c)
    c_invest = similar(c)
    c_zero = similar(c)
    V_zero = similar(c)
    V_invest = similar(c)
    diff_array = similar(c)
    mask = zeros(N, M)
    V = value_function_two_assets(c, φ, h, V)

    @showprogress 1 "Solving EGM" for s in 1:1000
        c_prev = copy(c)

        itp_c = [extrapolate(interpolate((W,), c[:, j], Gridded(Linear())), Line()) for j in 1:M]
        itp_V = [extrapolate(interpolate((W,), V[:, j], Gridded(Linear())), Line()) for j in 1:M]

        # Step 1: choose optimal phi
        for i in 1:N
            for j in 1:M
                f0 = phi_resid(0.0, i, j, itp_c, itp_V, W, z, π, d_nodes, w_d, h)
                f1 = phi_resid(1.0, i, j, itp_c, itp_V, W, z, π, d_nodes, w_d, h)

                if f0 > 0 && f1 < 0
                    φ_sol = find_zero(φ -> phi_resid(φ, i, j, itp_c, itp_V, W, z, π, d_nodes, w_d, h),
                                        (0.0, 1.0);
                                        method = Roots.Brent(),
                                        x0 = φ[i, j],
                                        xtol = 1e-4)
                elseif f0 >= 0 && f1 >= 0
                    φ_sol = 1.0
                elseif f0 <= 0 && f1 <= 0
                    φ_sol = 0.0
                else
                    error("Unexpected FOC signs: f0 = $f0, f1 = $f1")
                end

                φ[i, j] = φ_sol
            end
        end

        # Step 2: portfolio returns in each scenario
        RET_nodes = [(1 .- φ) .* h.R .+ φ .* (1 .+ d) for d in d_nodes]

        itp_c = [extrapolate(interpolate((W,), c[:, j], Gridded(Linear())), Line()) for j in 1:M]
        itp_V = [extrapolate(interpolate((W,), V[:, j], Gridded(Linear())), Line()) for j in 1:M]
        EU_prim = zero(c)
        EV_Vpow_1mθ = zero(c)

        # Krok 3: consumption interpolation at W_ret = RET ⋅ W + z
        EU_prim = zero(c)
        for m in 1:MD
            RET = RET_nodes[m]
            W_ret = (W * ones(1, M)) .* RET

            for j in 1:M
                @inbounds for i in 1:N
                    c_ij = itp_c[j](W_ret[i, j])
                    c_ret[i, j] = clamp(c_ij, 1e-8, W_ret[i, j] + z[j] - h.κ)
                    V_ret[i, j] = itp_V[j](W_ret[i, j])
                end
            end

            EV_Vpow_1mθ .+= w_d[m] .* V_ret.^(1 - h.θ) * π'
            U_prim = h.A .* ((W_ret .+ h.W_bar).^(-h.η))

            u_prim_m = RET .* ((1 - h.γ) .* c_ret.^(-h.ρ) .+ h.γ * U_prim) .* V_ret.^(h.ρ - h.θ)
            EU_prim .+= w_d[m] .* u_prim_m * π'
        end

        # Step 4: Inverted Euler Equation
        c_hat .= (h.β .* EU_prim ./ EV_Vpow_1mθ.^((h.ρ - h.θ)/(1 - h.θ))) .^(-1 / h.ρ)

        # Step 5: Budget constraint: W = W' - z + c
        W_hat .= W .- h.κ .+ c_hat .- z'

        # Step 6: new consumption function interpolation
        itp_c = [extrapolate(interpolate((W_hat[:, j],), c_hat[:, j], Gridded(Linear())), Line()) for j in 1:M]
        for j in 1:M
            @inbounds for i in 1:N
                c_ij = itp_c[j](W[i])
                c_invest[i, j] = clamp(c_ij, 1e-12, W[i] + z[j])
            end
        end

        V_invest = value_function_one_update(c, φ, h, V, 1)

        # Step 7: computing policy and value function without risky asset investment

        W_ret = (W * ones(1, M)) .* h.R

        for j in 1:M
            @inbounds for i in 1:N
                c_ij = itp_c[j](W_ret[i, j])
                c_ret[i, j] = clamp(c_ij, 1e-8, W_ret[i, j] + z[j])
                V_ret[i, j] = itp_V[j](W_ret[i, j])
            end
        end

        EV_Vpow_1mθ = V_ret.^(1 - h.θ) * π'
        U_prim = h.A .* ((W_ret .+ h.W_bar).^(-h.η))

        u_prim = h.R .* ((1 - h.γ) .* c_ret.^(-h.ρ) .+ h.γ * U_prim) .* V_ret.^(h.ρ - h.θ)
        EU_prim = u_prim * π'
        c_hat .= (h.β .* EU_prim ./ EV_Vpow_1mθ.^((h.ρ - h.θ)/(1 - h.θ))) .^(-1 / h.ρ)

        W_hat .= W .+ c_hat .- z'

        itp_c = [extrapolate(interpolate((W_hat[:, j],), c_hat[:, j], Gridded(Linear())), Line()) for j in 1:M]
        for j in 1:M
            @inbounds for i in 1:N
                c_ij = itp_c[j](W[i])
                c_zero[i, j] = clamp(c_ij, 1e-8, W[i] + z[j])
            end
        end

        V_zero = value_function_one_update(c_zero, φ, h, V, 0)

        # Step 9: finding optimal policy and value function update

        mask_higher = V_zero .> V_invest
        mask_diff = abs.(V_zero .- V_invest) ./ V_invest .> ε_κ
        mask = mask_diff .* mask_higher + (1 .- mask_diff) .* mask
        c = mask .* c_zero .+ (1 .- mask) .* c_invest
        φ = (1 .- mask) .* φ

        V = value_function_two_assets(c, φ, h, V)

        # Step 10: checking the differences and convergence

        diff_array = abs.(c - c_prev)
        diff = maximum(diff_array)

        if diff < ε 
            println("Converged in $s iterations")
            break
        elseif s == 1000
            println("EGM not converged in $s iterations")
            idx = findall(x -> x > ε, diff_array)
            @show diff
            @show idx
        end
    end

    mask = φ .> 0.0
    k_prim = (1 .- φ) .* ((W .+ z') .- c .- mask .* h.κ)
    S = φ .* ((W .+ z') .- c .- mask .* h.κ)

    return k_prim, S, c, φ, V
end


