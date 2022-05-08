function CR3BPdynamics(rv) #Three body dynamics in Sun-Earth System

    μ = 3.040423398444176e-6
    
    r₁³= ((rv[1] + μ)^2.0     + rv[2]^2.0 + rv[3]^2.0)^1.5; # distance to m1, LARGER MASS
    r₂³= ((rv[1] - 1 + μ)^2.0 + rv[2]^2.0 + rv[3]^2.0)^1.5; # distance to m2, smaller mass
    # r₁³= ((x + μ)^2     + y^2 + z^2)^1.5; # distance to m1, LARGER MASS
    # r₂³= ((x - 1 + μ)^2 + y^2 + z^2)^1.5; # distance to m2, smaller mass

#     rvdot = zeros(6)
    rvdot = zeros(eltype(rv),6)
    rvdot[1:3] = [rv[4];rv[5];rv[6]]
    rvdot[4] = -((1.0 - μ)*(rv[1] + μ)/r₁³) - (μ*(rv[1] - 1.0 + μ)/r₂³) + 2.0*rv[5] + rv[1];
    rvdot[5] = -((1.0 - μ)*rv[2]      /r₁³) - (μ*rv[2]          /r₂³) - 2.0*rv[4] + rv[2];
    rvdot[6] = -((1.0 - μ)*rv[3]      /r₁³) - (μ*rv[3]          /r₂³);
    return rvdot
end

#rk4 function
function rk4(f, x, h)

    f1 = f(x)
    f2 = f(x + 0.5*h*f1)
    f3 = f(x + 0.5*h*f2)
    f4 = f(x + h*f3)
    # TODO: implement
    # xnext = zero(x)
    xnext = x + (h/6.0)*(f1 + 2*f2 + 2*f3 + f4)
    
    return xnext
end