
abstract type AbstractParams{T} end

@kwdef struct ΛCDMParams{T} <: AbstractParams{T}

    h  :: T = 0.7       # hubble factor
    Ωr :: T = 5.042e-5  # radiation density
    Ωb :: T = 0.046     # baryon density
    Ωm :: T = 0.224     # matter density
    n  :: T = 1.0       # spectral index
    Yp :: T = 0.0       # primordial helium fraction, currently unused
    Nν :: T = 3.046     # effective number of relativisic species (PDG25 value)

    # solver options, indexed by the function theyre used in
    opts = (
        η = (rtol=1e-5,),
        boltsolve = (alg=Rodas5(), reltol=1e-11),
    )
    
end
