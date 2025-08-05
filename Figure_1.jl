using Interpolations, Plots, FastGaussQuadrature, Parameters, LinearAlgebra, Distributions, Roots, SparseArrays

include("functions_stylized.jl")
include("plotting.jl")

pyplot()

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

    W::Vector{Float64} # Grid for wealth
    z::Vector{Float64} #Grid for productivity
    π::Array{Float64} #Transition matrix for productivity process
end

N = 800
M = 3
ε_hh = 1e-8

h = Household(
    θ = 4.0,
    η = 4.0,
    A = 0.0,
    β = 0.97,
    γ = 0.0,

    σ_w = 0.1,
    ρ_w = 0.9,
    w = 1,
    R = 1.02,

    W = collect(LinRange(0.01, 100.0, N)),
    z = zeros(3),
    π = zeros(2,4)
)

h.π, z = rouwenhorst(M, h.ρ_w, h.σ_w)
z = exp.(z)
z_stat = stat_dist(h.π)
h.z = z / dot(z, z_stat)

#|==================================================================================
#|----------------------------Computing results-------------------------------------
#|==================================================================================

c_init = 0.2 * collect(LinRange(0.01, 100.0, N)* ones(1, M))

@time _, c1 = HH_Ayiagari_vec(c_init, M, N, ε_hh, h)

h.γ = 0.05
@time _, c2 = HH_Ayiagari_vec(c_init, M, N, ε_hh, h)

h.A = 4.0
@time _, c3 = HH_Ayiagari_vec(c_init, M, N, ε_hh, h)

h.η = 1.1
@time _, c4 = HH_Ayiagari_vec(c_init, M, N, ε_hh, h)

plt = fig_1(h, [c1, c2, c3, c4], custom_colors)
savefig(plt, "figures/figure_1.pdf")
