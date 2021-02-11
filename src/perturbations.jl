# these types and functions integrate the Boltzmann hierarchy through time

abstract type PerturbationIntegrator end
struct BasicNewtonian <: PerturbationIntegrator end

# a container for everything needed to integrate a hierarchy at wavenumber k
@kwdef struct Hierarchy{
    ℓmax_γ, 
    ℓmax_ν,
    PI <: PerturbationIntegrator, 
    CP <: AbstractParams,
    BG <: AbstractBackground, 
    IH <: AbstractIonizationHistory,
    T,
}
    integrator :: PI
    𝕡          :: CP
    bg         :: BG
    ih         :: IH
    k          :: T
end

function Hierarchy(integrator::PI, 𝕡::CP, bg::BG, ih::IH, k::T; ℓmax_γ=8, ℓmax_ν=10) where {PI,CP,BG,IH,T}
    Hierarchy{ℓmax_γ,ℓmax_ν,PI,CP,BG,IH,T}(integrator, 𝕡, bg, ih, k)
end


@⌛ function boltsolve(hierarchy::Hierarchy)
    @unpack 𝕡, bg = hierarchy
    xᵢ = first(bg.x_grid)
    u₀ = initial_conditions(xᵢ, hierarchy)
    prob = ODEProblem(hierarchy_ode, u₀, (xᵢ, zero(xᵢ)), hierarchy)
    sol = solve(prob, 𝕡.opts.boltsolve.alg; 𝕡.opts.boltsolve.reltol, saveat=hierarchy.bg.x_grid, dense=false)
    return sol
end


function unpack(u)
    @unpack Φ, δ, δ_b, v, v_b, Θ, Θᵖ, 𝒩 = u
    # make these zero-based indexing, since ComponentArrays components
    # currently lose this information
    Θ  = OffsetArray(Θ,  -1)
    Θᵖ = OffsetArray(Θᵖ, -1)
    𝒩  = OffsetArray(𝒩,  -1)
    Φ, δ, δ_b, v, v_b, Θ, Θᵖ, 𝒩
end

# BasicNewtonian comes from Callin+06 and the Dodelson textbook (dispatches on hierarchy.integrator)
function hierarchy_ode(u, hierarchy::Hierarchy{ℓmax_γ,ℓmax_ν,BasicNewtonian,<:AbstractParams{T}}, x) where {ℓmax_γ,ℓmax_ν,T}
    
    Θ′  = OffsetVector(@SVector(zeros(T,ℓmax_γ+1)), -1)
    Θᵖ′ = OffsetVector(@SVector(zeros(T,ℓmax_γ+1)), -1)
    𝒩′  = OffsetVector(@SVector(zeros(T,ℓmax_ν+1)), -1)

    @unpack k, 𝕡, bg, ih = hierarchy
    Ωr, Ωb, Ωm, Nν, H₀² = 𝕡.Ωr, 𝕡.Ωb, 𝕡.Ωm, 𝕡.Nν, bg.H₀^2 #add Nν≡N_eff
    ℋₓ, ℋₓ′, ηₓ, τₓ′, τₓ″ = bg.ℋ(x), bg.ℋ′(x), bg.η(x), ih.τ′(x), ih.τ″(x)

    a = x2a(x)
    R = 4Ωr / (3Ωb * a)
    Ω_ν =  7Nν/8 *(4/11)^(4/3) *Ωr

    Φ, δ, δ_b, v, v_b, Θ, Θᵖ, 𝒩 = unpack(u)

    # metric perturbations
    #Ψ = -Φ - 12H₀² / k^2 / a^2 * (Ωr * Θ[2])
    Ψ = -Φ - 12H₀² / k^2 / a^2 * (Ωr * Θ[2] + Ω_ν * 𝒩[2]) #add rel quadrupole
    Φ′ = Ψ - k^2 / (3ℋₓ^2) * Φ + H₀² / (2ℋₓ^2) * (
        Ωm * a^(-1) * δ + Ωb * a^(-1) * δ_b + 4Ωr * a^(-2) * Θ[0]
        + 4Ω_ν * a^(-2) * 𝒩[0]) #add rel monopole on this line

    # matter
    δ′ = k / ℋₓ * v - 3Φ′
    v′ = -v - k / ℋₓ * Ψ
    δ_b′ = k / ℋₓ * v_b - 3Φ′
    v_b′ = -v_b - k / ℋₓ * Ψ + τₓ′ * R * (3Θ[1] + v_b)

    # relativistic neutrinos (massless)
    @set! 𝒩′[0] = -k / ℋₓ * 𝒩[1] - Φ′
    @set! 𝒩′[1] = k/(3ℋₓ)*𝒩[0] - 2*k/(3ℋₓ)*𝒩[2] + k/(3ℋₓ)*Ψ
    for ℓ in 2:(ℓmax_ν-1) #ℓmax_ν same as ℓᵧ for massless nu for now
        @set! 𝒩′[ℓ] =  k / ((2ℓ+1) * ℋₓ) *( ℓ*𝒩[ℓ-1] - (ℓ+1)*𝒩[ℓ+1])
    end
    #truncation
    @set! 𝒩′[ℓmax_ν] =  k / ℋₓ  * 𝒩[ℓmax_ν-1] - (ℓmax_ν+1)/(ℋₓ *ηₓ) * 𝒩[ℓmax_ν]#Callin 06

    # WIP: massive nu
    # # neutrinos (massive, MB 57) - change convention
    # #units not yet right
    # Ψ_ν′[0] = -q*k / ϵ * Ψ_ν[1] - Φ′ *dnlnf0dlnq #FIXME dnln, def Psi,IC, Einstein int
    # Ψ_ν′[1] = q*k/(3 ϵ)*(𝒩[0] - 2*Ψ_ν[2]) - This k*ϵ / (3*q)*Ψ *dnlnf0dlnq
    # I think can't mutate the u variable...
    # Ψ_ν[ℓmax_νm] = (2*ℓmax_νm+1)*ϵ/(q*k*ηₓ)*Ψ_ν[ℓmax_νm] - Ψ_ν[ℓmax_νm-1] #truncation of MB (51)
    # for ℓ in 2:(ℓmax_νm-1) #ℓmax_νm should be smaller than massless case
    #     Ψ_ν′[ℓ] =  q*k / ((2ℓ+1) * ϵ) *( ℓ*Ψ_ν[ℓ-1] - (ℓ+1)*Ψ_ν[ℓ+1])
    # end

    # photons
    Π = Θ[2] + Θᵖ[2] + Θᵖ[0]
    @set! Θ′[0] = -k / ℋₓ * Θ[1] - Φ′
    @set! Θ′[1] = k / (3ℋₓ) * Θ[0] - 2k / (3ℋₓ) * Θ[2] + k / (3ℋₓ) * Ψ + τₓ′ * (Θ[1] + v_b/3)
    for ℓ in 2:(ℓmax_γ-1)
        @set! Θ′[ℓ] = ℓ * k / ((2ℓ+1) * ℋₓ) * Θ[ℓ-1] -
            (ℓ+1) * k / ((2ℓ+1) * ℋₓ) * Θ[ℓ+1] + τₓ′ * (Θ[ℓ] - Π * δ_kron(ℓ, 2) / 10)
    end

    # polarized photons
    @set! Θᵖ′[0] = -k / ℋₓ * Θᵖ[1] + τₓ′ * (Θᵖ[0] - Π / 2)
    for ℓ in 1:(ℓmax_γ-1)
        @set! Θᵖ′[ℓ] = ℓ * k / ((2ℓ+1) * ℋₓ) * Θᵖ[ℓ-1] -
            (ℓ+1) * k / ((2ℓ+1) * ℋₓ) * Θᵖ[ℓ+1] + τₓ′ * (Θᵖ[ℓ] - Π * δ_kron(ℓ, 2) / 10)
    end

    # photon boundary conditions: diffusion damping
    @set! Θ′[ℓmax_γ]  = k / ℋₓ * Θ[ℓmax_γ-1]  - (ℓmax_γ + 1) / (ℋₓ * ηₓ) + τₓ′ * Θ[ℓmax_γ]
    @set! Θᵖ′[ℓmax_γ] = k / ℋₓ * Θᵖ[ℓmax_γ-1] - (ℓmax_γ + 1) / (ℋₓ * ηₓ) + τₓ′ * Θᵖ[ℓmax_γ]

    ComponentArray{SVector}(;Φ=Φ′,δ=δ′,δ_b=δ_b′,v=v′,v_b=v_b′,Θ=Θ′,Θᵖ=Θᵖ′,𝒩=𝒩′)
end

# BasicNewtonian Integrator (dispatches on hierarchy.integrator)
@⌛ function initial_conditions(xᵢ, hierarchy::Hierarchy{ℓmax_γ,ℓmax_ν,BasicNewtonian,<:AbstractParams{T}}) where {ℓmax_γ,ℓmax_ν,T}

    Θ  = OffsetVector(@SVector(zeros(T,ℓmax_γ+1)), -1)
    Θᵖ = OffsetVector(@SVector(zeros(T,ℓmax_γ+1)), -1)
    𝒩  = OffsetVector(@SVector(zeros(T,ℓmax_ν+1)), -1)
    
    @unpack 𝕡, k, bg, ih = hierarchy
    ℋₓ, ℋₓ′, ηₓ, τₓ′, τₓ″ = bg.ℋ(xᵢ), bg.ℋ′(xᵢ), bg.η(xᵢ), ih.τ′(xᵢ), ih.τ″(xᵢ)
    H₀², aᵢ² = bg.H₀^2, exp(xᵢ)^2

    # metric and matter perturbations
    Φ   = 1.0
    δ   = 3Φ / 2
    δ_b = δ
    v   = k / (2ℋₓ) * Φ
    v_b = v

    # photon hierarchy
    @set! Θ[0] = Φ / 2
    @set! Θ[1] = -k * Φ / (6ℋₓ)
    @set! Θ[2] = -8k / (15ℋₓ * τₓ′) * Θ[1]
    for ℓ in 3:ℓmax_γ
        @set! Θ[ℓ] = -ℓ/(2ℓ+1) * k/(ℋₓ * τₓ′) * Θ[ℓ-1]
    end

    @set! Θᵖ[0] = (5/4) * Θ[2]
    @set! Θᵖ[1] = -k / (4ℋₓ * τₓ′) * Θ[2]
    @set! Θᵖ[2] = (1/4) * Θ[2]
    for ℓ in 3:ℓmax_γ
        @set! Θᵖ[ℓ] = -ℓ/(2ℓ+1) * k/(ℋₓ * τₓ′) * Θᵖ[ℓ-1]
    end

    # neutrino hierarchy
    # for now we assume xᵢ is before neutrinos decouple
    f_ν = 1/(1 + 1/(7𝕡.Nν/8 *(4/11)^(4/3)))
    @set! 𝒩[0] = Θ[0]
    @set! 𝒩[1] = Θ[1]
    @set! 𝒩[2] = - (k^2 *aᵢ²*Φ) / (12H₀²) * 1 / (1 + 5f_ν/2) #Callin06 (71)
    for ℓ in 3:ℓmax_ν
        @set! 𝒩[ℓ] = k/((2ℓ+1)ℋₓ) * 𝒩[ℓ-1] #approximation of Callin06 (72)
    end

    #WIP: massive nu
    #FIXME: nonrelativistic transition for massive species, needs to go in bg
    #^this will have to wait for m_ν to be added to 𝕡s
    #below notation is not right yet, starting from MB
    #x_nr = m_ν/5.3e-4 -1 #m_ν in eV (PDG26-pg3)
    # same as photons for 0,1
    # # massive #FIXME ingegrate the q moments, get dlnf0dlnq,define ϵ
    # σ_ν= (k*ηₓ)^2 *Ψ/ 15
    # Ψ_ν[0] = -δ_ν *dlnf0dlnq
    # Ψ_ν[1] = -ϵ/(3*q*k) θ_ν *dlnf0dlnq  #change θ to -k/ℋₓ v
    # Ψ_ν[2] = -σ_ν/2  *dlnf0dlnq
    #ignore ℓ>2, small

    ComponentArray{SVector}(;Φ,δ,δ_b,v,v_b,Θ,Θᵖ,𝒩)

end

# TODO: this could be extended to any Newtonian gauge integrator if we specify the
# Bardeen potential Ψ and its derivative ψ′ for an integrator, or we saved them
@⌛ function source_function(du, u, hierarchy::Hierarchy{T, BasicNewtonian}, x) where T
    # compute some quantities
    k, ℓᵧ, 𝕡, bg, ih = hierarchy.k, hierarchy.ℓᵧ, hierarchy.𝕡, hierarchy.bg, hierarchy.ih
    H₀² = bg.H₀^2
    ℋₓ, ℋₓ′, ℋₓ″ = bg.ℋ(x), bg.ℋ′(x), bg.ℋ″(x)
    τₓ, τₓ′, τₓ″ = ih.τ(x), ih.τ′(x), ih.τ″(x)
    g̃ₓ, g̃ₓ′, g̃ₓ″ = ih.g̃(x), ih.g̃′(x), ih.g̃″(x)
    a = x2a(x)

    Θ, Θᵖ, 𝒩, Φ, δ, v, δ_b, v_b = unpack(u, hierarchy)  # the Θ, Θᵖ are mutable views (see unpack)
    Θ′, Θᵖ′, 𝒩′, Φ′, δ′, v′, δ_b′, v_b′ = unpack(du, hierarchy)

    # recalulate these since we didn't save them
    Ψ = -Φ - 12H₀² / k^2 / a^2 * 𝕡.Ωr * Θ[2]
    Ψ′ = -Φ′ - 12H₀² / k^2 / a^2 * 𝕡.Ωr * (Θ′[2] - 2 * Θ[2])
    Π = Θ[2] + Θᵖ[2] + Θᵖ[0]
    Π′ = Θ′[2] + Θᵖ′[2] + Θᵖ′[0]

    term1 =  g̃ₓ * (Θ[0] + Ψ + Π/4) + exp(-τₓ) * (Ψ′ - Φ′)
    term2 = (-1/k) * (ℋₓ′ * g̃ₓ * v_b + ℋₓ * g̃ₓ′ * v_b + ℋₓ * g̃ₓ * v_b′)
    Π″ = 2k / (5ℋₓ) * (-ℋₓ′ / ℋₓ * Θ[1] + Θ′[1]) + (3/10) * (τₓ″ * Π + τₓ′ * Π′) -
        3k / (5ℋₓ) * (-ℋₓ′ / ℋₓ * (Θ[3] + Θᵖ[1] + Θᵖ[3]) + (Θ′[3] + Θᵖ′[1] + Θᵖ′[3]))
    term3 = (3/(4k^2)) * (
        (ℋₓ′^2 + ℋₓ * ℋₓ″) * g̃ₓ * Π + 3 * ℋₓ * ℋₓ′ * (g̃ₓ′ * Π + g̃ₓ * Π′) +
        ℋₓ^2 * (g̃ₓ″ * Π + 2g̃ₓ′ * Π′ + g̃ₓ * Π″))
    return term1 + term2 + term3
end
