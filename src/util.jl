
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
