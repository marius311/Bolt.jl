
# utilities for x ↔ scale factor ↔ redshift
a2z(a::T) where T = one(T)/a - one(T)
z2a(z::T) where T = one(T)/(z + one(T))
a2x(a) = log(a)
x2a(x) = exp(x)
z2x(z) = a2x(z2a(z))
x2z(x) = a2z(x2a(x))

# utility function for constructing an interpolator
spline(x, y) = scale(interpolate(y, BSpline(Cubic(Line(OnGrid())))), x)
spline_∂ₓ(f, x_grid) = spline(x_grid, [Interpolations.gradient(f, x)[1] for x in x_grid])
spline_∂ₓ²(f, x_grid) = spline(x_grid, [Interpolations.hessian(f, x)[1] for x in x_grid])

δ_kron(i, j) = (i == j) ? 1 : 0



"""

    @⌛ code ...
    @⌛ function_definition() = .... 

Label a section of code to be timed. The first form uses the code
itselfs as a label, the second uses the function name, and its the
body of the function which is timed. 

To run the timer and print output, returning the result of the
calculation, use

    @show⌛ run_code()

Timing uses `TimerOutputs.get_defaulttimer()`. 
"""
macro ⌛(ex)
    source_str = last(splitpath(string(__source__.file)))*":"*string(__source__.line)
    try
        # function definition
        sdef = splitdef(ex)
        sdef[:body] = quote
            $TimerOutputs.@timeit $("$(string(sdef[:name]))(…)  ($source_str)") $(sdef[:body])
        end
        esc(combinedef(sdef))
    catch
        # anything else
        :(@timeit $("$(Base._truncate_at_width_or_chars(string(prewalk(rmlines,ex)),26))  ($source_str)") $(esc(ex)))
    end
end


"""
See [`@⌛`](@ref)
"""
macro show⌛(ex)
    quote
        reset_timer!($TimerOutputs.get_defaulttimer())
        result = $(esc(ex))
        show($TimerOutputs.get_defaulttimer())
        result
    end
end


# needed until fix to https://github.com/jw3126/Setfield.jl/issues/153
@inline function Setfield.setindex(A::OffsetVector{T,S}, val, i::Int) where {T,S<:StaticArray}
    @boundscheck checkbounds(A, i)
    OffsetVector{T,S}(setindex(parent(A), val, OffsetArrays.parentindex(Base.axes1(A), i)), A.offsets)
end


# performance optimization for ComponentArray constructors when the
# names/indices of components are known at compile-time (ie when the
# components are Numbers, StaticArrays, or OffsetArrays of StaticArrays).
# this uses a generated function to move the (costly) construction of
# the Axis object to compile-time
@generated function ComponentArrays.ComponentArray{SVector}(nt::NamedTuple{Names,Tup}) where {Names, Tup<:NTuple{<:Any,Union{Number,StaticVector,OffsetVector{<:Any,<:StaticVector}}}}
    # construct a dummy ComponentArray at compile-time as an easy way
    # to get the Axis object (eltype doesnt matter so use 0)
    _zero(::Type{<:Number}) where {T<:Number} = 0
    _zero(::Type{<:Union{SVector{N},OffsetVector{<:Any,<:SVector{N}}}}) where {N} = fill(0,N)
    Axis = getaxes(ComponentArray(;(Names .=> map(_zero, Tup.parameters))...))
    # create constructor expression
    _splat(name, ::Type{<:Number}) = :(nt.$name)
    _splat(name, ::Type{<:Union{SVector,OffsetVector{<:Any,<:SVector}}}) = :((nt.$name)...)
    quote
        ComponentArray(SVector(tuple($(map(_splat, Names, Tup.parameters)...))), $Axis)
    end
end
