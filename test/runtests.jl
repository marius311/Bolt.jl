using Bolt
using Test
using DelimitedFiles
using LinearAlgebra

@testset "FFTLog" begin
    N = 64
    μ = 0
    q = 0.0
    r₀ = 1.0
    L = 8.0
    Nhalf = N ÷ 2
    n = range(-Nhalf,Nhalf,length=N)
    r = r₀ .* 10 .^ (n .* L ./ N )
    pl = Bolt.plan_fftlog(r, μ, q, 1.0; kropt=true)
    aₙ = r .^ (μ + 1) .* exp.(-r.^2 / 2)
    y = similar(r, ComplexF64)
    fftdata = readdlm("data/fftlog_example.txt", ' ', Float64, '\n')

    # test forward
    mul!(y, pl, aₙ)
    f_ref = fftdata[:,2]
    @test all(abs.(y .- f_ref) .< 1e-15)
    @test isapprox(y, f_ref)

    # test backward
    y2 = similar(r, ComplexF64)
    ldiv!(y2, pl, y)
    @test all(abs.(y2 .- aₙ) .< 1e-15)

end

# test U with more naive result
# U_μ(μ, x) = 2^x * gamma((μ + 1 + x)/2) / gamma((μ + 1 - x)/2)

## hand-written ℋ derivative for testing
# function ℋ′(x, par::AbstractCosmoParams)
#     a = x2a(x)
#     return -H₀(par) * (2par.Ωr + (par.Ωb + par.Ωm) * a - 2ΩΛ(par) * a^4) /
#         (2 * a * √(par.Ωr + (par.Ωb + par.Ωm) * a + ΩΛ(par) * a^4))
# end
