using Bolt
using ForwardDiff
using PyPlot
using BenchmarkTools
using OrdinaryDiffEq
using LinearAlgebra
using ComponentArrays
using JLD2
##
T = Float64
𝕡 = ΛCDMParams(;
    opts = (
        η = (rtol=1e-5,),
        boltsolve = (alg=Rodas5(), reltol=0.005),
    )
)
bg = Background(𝕡)
ih = IonizationHistory(Peebles(), 𝕡, bg)
Nk = 100
k, = k_grid = quadratic_k(0.1bg.H₀, 1000bg.H₀, Nk)
ℓmax_γ = 8
ℓmax_ν = 10
hierarchy = Hierarchy(; integrator=BasicNewtonian(), 𝕡, bg, ih, k, ℓmax_γ, ℓmax_ν)

x₀ = first(bg.x_grid)
u₀ = Bolt.initial_conditions(x₀, hierarchy)

@show⌛ boltsolve(hierarchy)
sf = @show⌛ source_grid(𝕡, bg, ih, k_grid, BasicNewtonian())

##
ℓs = 100:50:1200
Cℓs = cltt(ℓs, 𝕡, bg, ih, sf)
##
clf()
plot(ℓs, Cℓs .* ℓs.^2, "-")
plot(ℓs, Cℓs₀ .* ℓs.^2, "-")
twinx()
plot(ℓs, @. Cℓs/Cℓs₀ - 1)
ylabel(raw"$\ell^2 C_{\ell}^{TT}$")
xlabel(raw"$\ell$")
gcf()
##
@load "cls.jld2"