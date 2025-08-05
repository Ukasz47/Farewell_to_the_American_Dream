using Interpolations, Plots, FastGaussQuadrature, Parameters, LinearAlgebra, Distributions, Roots, SparseArrays, Measures

include("functions.jl")
include("analysis.jl")
include("plotting.jl")

#|==================================================================================
#|--------------------------Defining parameters-------------------------------------
#|==================================================================================

@with_kw mutable struct Household
    θ::Float64 # Utility function risk aversion parameter
    ρ::Float64 # Intertemporal elasticity of substitution (inverted)
    η::Float64 # Assets in utility parameter
    A::Float64 # Scaling factor for utility of bequests
    W_bar::Float64 # "Subsistence" level for bequests
    β::Float64 # Discount factor
    γ::Float64 # Mortality rate
    ξ::Float64 # "Fortune destruction" rate

    σ_z::Float64 # Productivity standard deviation
    ρ_z::Float64 # Productivity persistence
    ρ_top::Float64 # Persistence of top income state
    w::Float64 # Wage
    z_top::Float64 # Productivity of top 5% (comparing to average)
    R::Float64 # Gross interest rate
    E_d::Float64 # expected dividend
    σ_d::Float64 # dividend std. dev
    p_B::Float64 # probability of bankrupcy
    κ::Float64 # cost of participating in asset markets

    W::Vector{Float64} # Grid for wealth
    z::Vector{Float64} #Grid for productivity
    π::Array{Float64} #Transition matrix for productivity process
    d_nodes::Vector{Float64} #Grid for return
    w_d::Vector{Float64} #Probabilities for return

    c_init::Array{Float64} #Initial guess for consumption
    φ_init::Array{Float64} #Initial guess for risky share
    V_init::Array{Float64} #Initial guess for value function
end


ε_hh = 1e-6
N = 500
M = 4
MD = 6
mean_wealth = 6.9

h = Household(
    θ = 4.0,
    ρ = 4.0,
    η = 1.6,
    A = 4.0,
    W_bar = 10.0,
    β = 0.955,
    γ = 0.02,
    ξ = 0.01,

    σ_z = 0.1,
    ρ_z = 0.95,
    ρ_top = 0.95,
    w = 1,
    z_top = 6.0,
    R = 1.015,
    E_d = 0.065,
    σ_d = 0.2,
    p_B = 0.01,
    κ = 0.3,

    W = power_grid(0.01, 4000.0, N; ξ = 2.0),
    z = zeros(3),
    π = zeros(2,4),
    d_nodes = zeros(MD),
    w_d = zeros(MD),

    c_init = 0.2 * power_grid(0.01, 4000.0, N; ξ = 2.0) * ones(1, M),
    φ_init = 0.0332872 * ones(N, M),
    V_init = ones(N, M)
)

#|==================================================================================
#|----------------------------Computing results-------------------------------------
#|==================================================================================

h.π, h.z = income_process(M - 1, h.ρ_z, h.σ_z, h.ρ_top, h.z_top)
z_stat = stat_dist(h.π)

h.d_nodes, h.w_d = risk_process(MD, h)

@time k_prim, S, c, φ, V = HH_Ayiagari_two_assets_exo(M, N, MD, ε_hh, h)

W_stat, Wz_stat = get_distribution(mean_wealth, k_prim, S, h)

results_wealth = analyze_wealth_distribution(h, W_stat)
results_return = portfolio_stats_by_wealth_decile(h, Wz_stat, φ)
