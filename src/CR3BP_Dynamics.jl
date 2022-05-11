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

#dircol dynamics Hermite Simpson Spline
function dircol_dynamics(x1,x2,h)
    #Hermite-Simpson integration with first-order hold on u
    f1 = CR3BPdynamics(x1) #Timestep k
    f2 = CR3BPdynamics(x2) #Timestep k+1
    xm = 0.5.*(x1 + x2) .+ (h/8.0).*(f1 - f2) #
    ẋm = (-3.0/(2.0*h)).*(x1 - x2) .- 0.25.*(f1 + f2)
    fm = CR3BPdynamics(xm)
    return fm - ẋm
end

function CR3BPdynamicsWithUforA(X, U) #Three body dynamics in Sun-Earth System
    
    μ = 3.040423398444176e-6
    
    r₁³= ((X[1] + μ)^2.0     + X[2]^2.0 + X[3]^2.0)^1.5; # distance to m1, LARGER MASS
    r₂³= ((X[1] - 1.0 + μ)^2.0 + X[2]^2.0 + X[3]^2.0)^1.5; # distance to m2, smaller mass

    Xdot = zeros(eltype(X),6)
    Xdot[1:3] = [X[4];X[5];X[6]]
    Xdot[4] = -((1.0 - μ)*(X[1] + μ)/r₁³) - (μ*(X[1] - 1.0 + μ)/r₂³) + 2.0*X[5] + X[1]
    Xdot[4] = Xdot[4] + U[1]
    Xdot[5] = -((1.0 - μ)*X[2]      /r₁³) - (μ*X[2]          /r₂³) - 2.0*X[4] + X[2]
    Xdot[5] = Xdot[5] +  U[2]
    Xdot[6] = -((1.0 - μ)*X[3]      /r₁³) - (μ*X[3]          /r₂³)
    Xdot[6] = Xdot[6] + U[3]
    return Xdot
end


function CR3BPdynamicsWithUforB(X, U) #Three body dynamics in Sun-Earth System
    μ = 3.040423398444176e-6
    
    r₁³= ((X[1] + μ)^2.0     + X[2]^2.0 + X[3]^2.0)^1.5; # distance to m1, LARGER MASS
    r₂³= ((X[1] - 1.0 + μ)^2.0 + X[2]^2.0 + X[3]^2.0)^1.5; # distance to m2, smaller mass

    Xdot = zeros(eltype(U),6)
#     Xdot = zeros(6)
    Xdot[1:3] = [X[4];X[5];X[6]]
    Xdot[4] = -((1.0 - μ)*(X[1] + μ)/r₁³) - (μ*(X[1] - 1.0 + μ)/r₂³) + 2.0*X[5] + X[1]
    Xdot[4] = Xdot[4] + U[1]
    Xdot[5] = -((1.0 - μ)*X[2]      /r₁³) - (μ*X[2]          /r₂³) - 2.0*X[4] + X[2]
    Xdot[5] = Xdot[5] +  U[2]
    Xdot[6] = -((1.0 - μ)*X[3]      /r₁³) - (μ*X[3]          /r₂³)
    Xdot[6] = Xdot[6] + U[3]
    return Xdot
end

function rk4WithU(f, x, u, h) #the second input stays as u bc u is constant for the timestep, right?

    f1 = f(x,u)
    f2 = f(x + 0.5.*h.*f1,u)
    f3 = f(x + 0.5.*h.*f2,u)
    f4 = f(x + h.*f3,u)
    xnext = x + (h/6.0).*(f1 + 2.0 .*f2 + 2.0 .*f3 + f4)
    
    return xnext
end

function linearize!(A,B,xtraj, U)

    X = copy(xtraj)
    
    # loop over all the time steps in the reference trajectory
    for k = 1:Nt-1
        # evaluate the discrete jacobian at the current time step
        A[k] = ForwardDiff.jacobian(_X_->rk4WithU(CR3BPdynamicsWithUforA,_X_,U[k],utraj[k]),X[k])
        B[k] = ForwardDiff.jacobian(_U_->rk4WithU(CR3BPdynamicsWithUforB,X[k],_U_,utraj[k]),U[k])
    end
end

function calc_gains!(xtraj, A, B, P, R, K, Qf, Q) #input here is the full state output of ipopt
    
    #Implement Riccati recursion for TVLQR

    P[end] .= Qf
    for k = reverse(1:Nt-1) #does the recursion only work with uniform timesteps???
        K[k] .= (R + B[k]'P[k+1]*B[k])\(B[k]'P[k+1]*A[k])
        P[k] .= Q + A[k]'P[k+1]*A[k] - A[k]'P[k+1]*B[k]*K[k]
    end

    return nothing
end
