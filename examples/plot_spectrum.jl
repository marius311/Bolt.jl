using Bolt
using ForwardDiff
using PyPlot
using BenchmarkTools
using OrdinaryDiffEq
using NaturallyUnitful
using StaticArrays
using LinearAlgebra
using Roots
using JLD2
using ComponentArrays
using UnPack
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
k = quadratic_k(0.1bg.H₀, 1000bg.H₀, Nk)
hierarchy = Hierarchy(;integrator=BasicNewtonian(), 𝕡, bg, ih, k)

Nk = length(k)
x₀ = first(bg.x_grid)
u₀ = Bolt.initial_conditions!(
    ComponentArray(
        Φ   = similar(k, Nk),
        δ   = similar(k, Nk),
        δ_b = similar(k, Nk),
        v   = similar(k, Nk),
        v_b = similar(k, Nk),
        Θ   = similar(k, Nk, hierarchy.ℓmax_γ+1),
        Θᵖ  = similar(k, Nk, hierarchy.ℓmax_γ+1),
        𝒩   = similar(k, Nk, hierarchy.ℓmax_ν+1),
    ),
    hierarchy,
    x₀
)
Bolt.initial_conditions!(u₀, hierarchy, first(bg.x_grid))

prob = ODEProblem(Bolt.hierarchy!, u₀, (x₀, zero(T)), hierarchy)
prob = ODEProblem(Bolt.hierarchy!, u₀, (-T(20), -T(19)), hierarchy)
sol = @profview solve(prob, 𝕡.opts.boltsolve.alg; 𝕡.opts.boltsolve.reltol, saveat=hierarchy.bg.x_grid, dense=false)

Bolt.hierarchy!(similar(u₀), u₀, hierarchy, first(bg.x_grid))

1

using BlockBandedMatrices

ForwardDiff.jacobian!




ODEFunction

sol

Bolt.boltsolve(h)



sf = @show⌛ source_grid(𝕡, bg, ih, k_grid, BasicNewtonian())
sf = @profview source_grid(𝕡, bg, ih, k_grid, BasicNewtonian())
ℓs = 100:50:1200
Cℓs = cltt(ℓs, 𝕡, bg, ih, sf)
# Cℓs₀ = Cℓs;
@save "cls.jld2" Cℓs₀
##
clf()
plt.plot(ℓs, Cℓs₀ .* ℓs.^2, "-")
plt.plot(ℓs, Cℓs .* ℓs.^2, "-")
twinx()
plt.plot(ℓs, @. Cℓs / Cℓs₀ - 1)
ylabel(L"\ell^2 C_{\ell}^{TT}")
xlabel(L"\ell")
gcf()
##

f(u,p,t) = 1.01*u
u0 = 1/2
tspan = (0.0,1.0)
prob = ODEProblem(f,u0,tspan)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)


##


reltol = find_zero(
    function (reltol)
        @show reltol
        𝕡 = ΛCDMParams(;
            opts = (
                η = (rtol=1e-5,),
                boltsolve = (;alg=TRBDF2(), reltol),
            )
        )
        bg = Background(𝕡)
        ih = IonizationHistory(Peebles(), 𝕡, bg)
        k_grid = quadratic_k(0.1bg.H₀, 1000bg.H₀, 100)
        global t = @elapsed(sf = source_grid(𝕡, bg, ih, k_grid, BasicNewtonian()))
        global Cℓs = cltt(ℓs, 𝕡, bg, ih, sf)
        norm(@. Cℓs / Cℓs₀ - 1) - 0.01
    end,
    (1e-10, 1e-1),
    atol = 0.001
)
@show reltol, t


# Rodas5 0.0028 4.7s
# KenCarp4 1.2e-6 4.8s
#





