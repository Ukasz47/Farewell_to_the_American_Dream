include("toolbox.jl")

#|============================================================================================
#|--------------Functions for finding household optimum and general equilibrium---------------
#|============================================================================================

function HH_Ayiagari_vec(c_init, M, N, ε, h)
    π = h.π
    z = h.z * h.w

    W = h.W              
    c = c_init             
    c_hat = similar(c)
    W_hat = similar(c)
    k_prim = copy(c)

    @showprogress 1 "Solving EGM" for s in 1:2000
        c_prev = copy(c)

        W_eta = (W.^(-h.η)) * ones(1, M)
        E_u_prim = h.β .* h.R .* (π * ((1 .- h.γ) .* c.^(-h.θ) .+ h.γ .* h.A .* W_eta)')
        c_hat .= E_u_prim'.^(-1 / h.θ)

        W_hat = W ./ h.R .+ c_hat .- z'

        for j in 1:M
            itp_c = extrapolate(interpolate((W_hat[:, j],), c_hat[:, j], Gridded(Linear())), Line())
            for i in 1:N
                c_ij = itp_c(W[i])
                c[i, j] = clamp(c_ij, 1e-8, W[i] + z[j])
            end
        end

        k_prim .= (W .+ z') .- c

        if maximum(abs.(c - c_prev)) < ε
            println("Converged in $s iterations")
            break
        elseif s == 2000
            println("EGM not converged in $s iterations")
        end
    end

    return deepcopy(k_prim), deepcopy(c)
end

function phi_resid(
    φ::Float64,
    i::Int,
    j::Int,
    c_interp::AbstractVector,
    W::Vector{Float64},
    z::Vector{Float64},
    π::Matrix,
    d_nodes::Vector{Float64},
    w_d::Vector{Float64},
    h
)
    resid = 0.0

    for m in eachindex(d_nodes)
        RET = (1 - φ) * h.R + φ * (1 + d_nodes[m] )
        excess_return = (1 + d_nodes[m] - h.R)
        savings = W[i] + z[j] - c_interp[j](W[i])
        savings = clamp(savings, 1e-8, W[i] + z[j])
        W_ret = savings * RET

        for jp in eachindex(z)
            c_p = c_interp[jp](W_ret)
            c_p = clamp(c_p, 1e-8, W_ret + z[jp])
            u_p = (1 - h.γ) * c_p^(-h.θ)
            V_p = h.γ * h.A * W_ret^(-h.η)
            resid += π[j, jp] * w_d[m] * (u_p + V_p) * excess_return
        end
    end
    return resid
end

function HH_Ayiagari_two_assets_endo(M, N, MD, ε, h)
    π = h.π
    z = h.z * h.w
    d_nodes = h.d_nodes
    w_d = h.w_d
    W = h.W

    c = h.c_init
    φ = h.φ_init
    c_hat = similar(c)
    W_hat = similar(c)
    c_ret = similar(c)

    @showprogress 1 "Solving EGM" for s in 1:2000
        c_prev = copy(c)

        itp_c = [extrapolate(interpolate((W,), c[:, j], Gridded(Linear())), Line()) for j in 1:M]

        # Step 1: choose optimal phi
        for i in 1:N
            for j in 1:M
                f0 = phi_resid(0.0, i, j, itp_c, W, z, π, d_nodes, w_d, h)
                f1 = phi_resid(1.0, i, j, itp_c, W, z, π, d_nodes, w_d, h)

                if f0 > 0 && f1 < 0
                    φ_sol = find_zero(φ -> phi_resid(φ, i, j, itp_c, W, z, π, d_nodes, w_d, h),
                                      φ[i, j], Order1(), xtol=1e-4)
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

        # Krok 3: consumption interpolation at W_ret = RET ⋅ W + z
        EU_prim = zero(c)
        for m in 1:MD
            RET = RET_nodes[m]
            W_ret = (W * ones(1, M)) .* RET
            W_eta = W_ret.^(-h.η)

            for j in 1:M
                @inbounds for i in 1:N
                    c_ij = itp_c[j](W_ret[i, j])
                    c_ret[i, j] = clamp(c_ij, 1e-8, W_ret[i, j] + z[j])
                end
            end

            u_prim_m = RET .* ((1 .- h.γ) .* c_ret.^(-h.θ) .+ h.γ .* h.A .* W_eta)
            EU_prim .+= w_d[m] .* (π * u_prim_m')'
        end

        # Step 4: Inverted Euler Equation
        E_u_prim = h.β .* EU_prim
        c_hat .= E_u_prim.^(-1 / h.θ)

        # Step 5: Budget constraint: W = W' - z + c
        W_hat .= W .+ c_hat .- z'

        # Step 6: new consumption function interpolation
        itp_c = [extrapolate(interpolate((W_hat[:, j],), c_hat[:, j], Gridded(Linear())), Line()) for j in 1:M]
        for j in 1:M
            @inbounds for i in 1:N
                c_ij = itp_c[j](W[i])
                c[i, j] = clamp(c_ij, 1e-12, W[i] + z[j])
            end
        end

        if maximum(abs.(c - c_prev)) < ε
            println("Converged in $s iterations")
            break
        end
    end

    k_prim = (1 .- φ) .* ((W .+ z') .- c)
    S = φ .* ((W .+ z') .- c)

    return deepcopy(k_prim), deepcopy(S), deepcopy(c), deepcopy(φ)
end


