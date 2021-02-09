
TCA_condition(k, ℋₓ, τₓ′) = (abs(k / (ℋₓ * τₓ′)) < 0.1) & (abs(τₓ′) > 10.0)

# NOTE: NO NEUTRINOS 𝒩
function hierarchy!(du, u, p::AbstractCosmoParams, x)
    ℓᵧ, Ωr, Ωb, Ωm = p.ℓᵧ, p.Ωr, p.Ωb, p.Ωm
    H₀² = bg.H₀^2
    ℋₓ, ℋₓ′ = bg.ℋ(x), bg.ℋ′(x)
    τₓ′, τₓ′′ = ih.τ′(x), ih.τ′′(x)
    a = x2a(x)
    R = 4Ωr / (3Ωb * a)

    # get array views of photon perturbations
    Θ = OffsetVector(view(u, 1:(ℓᵧ+1)), 0:ℓᵧ)  # indicies 0 through ℓᵧ
    Θ′ = OffsetVector(view(du, 1:(ℓᵧ+1)), 0:ℓᵧ)
    Θᵖ = OffsetVector(view(u, (ℓᵧ+2):(2ℓᵧ+2)), 0:ℓᵧ)  # indicies 0 through ℓᵧ
    Θᵖ′ = OffsetVector(view(du, (ℓᵧ+2):(2ℓᵧ+2)), 0:ℓᵧ)
    Φ, δ, v, δ_b, v_b = u[(2ℓᵧ+3):(2ℓᵧ+7)]

    # metric perturbations
    Ψ = -Φ - 12H₀² / k^2 / a^2 * (Ωr * Θ[2])
    Φ′ = Ψ - k^2 / (3ℋₓ^2) * Φ + H₀² / (2ℋₓ^2) * (
        Ωm * a^(-1) * δ + Ωb * a^(-1) * δ_b + 4Ωr * a^(-2) * Θ[0])

    # matter
    δ′ = k / ℋₓ * v - 3Φ′
    v′ = -v - k / ℋₓ * Ψ
    δ_b′ = k / ℋₓ * v_b - 3Φ′

    Θ′[0] = -k / ℋₓ * Θ[1] - Φ′
    if TCA_condition(k, ℋₓ, τₓ′)
        Θ′[2] = 0.0  # could be solved for
        term_3Θ′_vb′ = (
            -((1-R)*τₓ′ + (1+R)*τₓ′′) * (3Θ[1] + v_b) - k * Ψ / ℋₓ
            + (1 - ℋₓ′ / ℋₓ) * (k / ℋₓ) * (-Θ[0] + 2Θ[2]) + k / ℋₓ * (-Θ′[0] + 2Θ′[2])
        ) / ((1+R)*τₓ′ + ℋₓ′ / ℋₓ - 1)
        v_b′ = (-v_b - k * Ψ / ℋₓ + R * (
            term_3Θ′_vb′ + k / ℋₓ * (-Θ[0] + 2Θ[2]) - k / ℋₓ * Ψ
        )) / (1+R)
        Θ′[1] = (term_3Θ′_vb′ - v_b′) / 3
    else
        v_b′ = -v_b - k / ℋₓ * Ψ + τₓ′ * R * (3Θ[1] + v_b)
        Θ′[1] = k / (3ℋₓ) * Θ[0] - 2k / (3ℋₓ) * Θ[2] + k / (3ℋₓ) * Ψ + τₓ′ * (Θ[1] + v_b/3)
    end

    # photons
    Π = Θ[2] + Θᵖ[2] + Θᵖ[0]
    for ℓ in 2:(ℓᵧ-1)
        Θ′[ℓ] = ℓ * k / ((2ℓ+1) * ℋₓ) * Θ[ℓ-1] -
            (ℓ+1) * k / ((2ℓ+1) * ℋₓ) * Θ[ℓ+1] + τₓ′ * (Θ[ℓ] - Π * δ_kron(ℓ, 2) / 10)
    end
    # polarized photons
    Θᵖ′[0] = -k / ℋₓ * Θᵖ[1] + τₓ′ * (Θᵖ[0] - Π / 2)
    for ℓ in 1:(ℓᵧ-1)
        Θᵖ′[ℓ] = ℓ * k / ((2ℓ+1) * ℋₓ) * Θᵖ[ℓ-1] -
            (ℓ+1) * k / ((2ℓ+1) * ℋₓ) * Θᵖ[ℓ+1] + τₓ′ * (Θᵖ[ℓ] - Π * δ_kron(ℓ, 2) / 10)
    end
    # photon hierarchy boundary conditions
    Θ′[ℓᵧ] = k / ℋₓ * Θ[ℓᵧ-1] - (ℓᵧ + 1) / (ℋₓ * bg.η(x)) + τₓ′ * Θ[ℓᵧ]
    Θᵖ′[ℓᵧ] = k / ℋₓ * Θᵖ[ℓᵧ-1] - (ℓᵧ + 1) / (ℋₓ * bg.η(x)) + τₓ′ * Θᵖ[ℓᵧ]

    du[(2ℓᵧ+3):(2ℓᵧ+7)] .= Φ′, δ′, v′, δ_b′, v_b′  # put non-photon perturbations back in
    return nothing
end
