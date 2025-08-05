using Interpolations, Plots, FastGaussQuadrature, Parameters, LinearAlgebra, Distributions, Roots, SparseArrays, Measures

include("functions_stylized.jl")
include("plotting.jl")

#|==================================================================================
#|--------------------------Defining parameters-------------------------------------
#|==================================================================================

@with_kw mutable struct Household
    θ::Float64 # Utility function risk aversion parameter
    η::Float64 # Assets in utility parameter
    A::Float64 # Scaling factor for utility of bequests
    β::Float64 # Discount factor
    γ::Float64 # Mortality rate

    σ_w::Float64 # Productivity standard deviation
    ρ_w::Float64 # Productivity persistence
    w::Float64 # Wage
    R::Float64 # Gross interest rate
    E_d::Float64 # expected dividend
    σ_d::Float64 # dividend std. dev

    W::Vector{Float64} # Grid for wealth
    z::Vector{Float64} #Grid for productivity
    π::Array{Float64} #Transition matrix for productivity process
    d_nodes::Vector{Float64} #Grid for return
    w_d::Vector{Float64} #Probabilities for return

    c_init::Array{Float64} #Initial guess for consumption
    φ_init::Array{Float64} #Initial guess for risky share
end

N = 1000
M = 3
MD = 11
ε_hh = 1e-8

h = Household(
    θ = 4.0,
    η = 1.1,
    A = 2.0,
    β = 0.97,
    γ = 0.0,

    σ_w = 0.2,
    ρ_w = 0.9,
    w = 1,
    R = 1.02,
    E_d = 0.02001,
    σ_d = 0.2,

    W = collect(LinRange(0.01, 1500.0, N)),
    z = zeros(3),
    π = zeros(2,4),
    d_nodes = zeros(MD),
    w_d = zeros(MD),

    c_init = 0.2 * collect(LinRange(0.01, 1500.0, N)* ones(1, M)),
    φ_init = zeros(N, M)
)

h.π, z = rouwenhorst(M, h.ρ_w, h.σ_w)
z = exp.(z)
z_stat = stat_dist(h.π)
h.z = z / dot(z, z_stat)

x_d, w_d = gausshermite(MD)
h.d_nodes = h.E_d * exp.(sqrt(2) * h.σ_d .* x_d .- 0.5 * h.σ_d^2)
h.w_d = w_d ./ sqrt(Base.MathConstants.pi)

#|==================================================================================
#|----------------------------Computing results-------------------------------------
#|==================================================================================

@time _, _, _, φ1 = HH_Ayiagari_two_assets_endo(M, N, MD, ε_hh, h)

h.γ = 0.02
@time _, _, _, φ2 = HH_Ayiagari_two_assets_endo(M, N, MD, ε_hh, h)

plt = fig_2(h, [φ1, φ2], custom_colors)
savefig(plt, "figures/figure_2.pdf")