include("Non_lin_exp.jl")

Base.@kwdef struct CR3BPProblem{T}
    model::CR3BP{T} = CR3BP()
    n::Int = 6                               # num states
    m::Int = 1                               # num controls
    N::Int = 101                             # horizon
    tf::T = 2.0                              # final time (sec)
    x0::SVector{6,T} = @SVector zeros(6)     # initial state
    xf::SVector{6,T} = SA[0,0,0,0,0,0]         # goal state- same as initial state? 
    Q::Diagonal{T, SVector{6,T}} = Diagonal(@SVector fill(1e-2, 6))
    R::Diagonal{T, SVector{1,T}} = Diagonal(@SVector fill(5e-0, 1))
    Qf::Diagonal{T, SVector{6,T}} = Diagonal(@SVector fill(1e1, 6))
end

    
@autodiff struct CR3BP{T} <: ContinuousDynamics 
    μ::T
    function CR3BP(μ)
        T = eltype(promote(μ))
        new{T}(μ)
    end
end

CR3BP(;μ = 9.537e-4) = CR3BP(μ)

function dynamics(model::CR3BP, x, u)
    μ = model.μ #Ratio of two primaries 
    
    r₁³= ((x[1] + μ)^2.0     + x[2]^2.0 + x[3]^2.0)^1.5; # distance to m1, LARGER MASS
    r₂³= ((x[1] - 1 + μ)^2.0 + x[2]^2.0 + x[3]^2.0)^1.5; # distance to m2, smaller mass
    # r₁³= ((x + μ)^2     + y^2 + z^2)^1.5; # distance to m1, LARGER MASS
    # r₂³= ((x - 1 + μ)^2 + y^2 + z^2)^1.5; # distance to m2, smaller mass

#     xdot = zeros(6)
    xdot = zeros(eltype(x),6)
    xdot[1:3] = [x[4];x[5];x[6]]
    xdot[4] = -((1.0 - μ)*(x[1] + μ)/r₁³) - (μ*(x[1] - 1.0 + μ)/r₂³) + 2.0*x[5] + x[1];
    xdot[5] = -((1.0 - μ)*x[2]      /r₁³) - (μ*x[2]          /r₂³) - 2.0*x[4] + x[2];
    xdot[6] = -((1.0 - μ)*x[3]      /r₁³) - (μ*x[3]          /r₂³);
    
    return xdot
end

function dynamics!(model::CR3BP, xdot, x, u)
    xdot .= dynamics(model, x, u)
end


# function get_objective(prob::CartpoleProblem)
#     costfun = LQRCost(prob.Q, prob.R, prob.xf)
#     costterm = LQRCost(prob.Qf, prob.R, prob.xf)
#     obj = push!(fill(costfun, prob.N-1), costterm)
#     return obj
# end

"""
    get_initial_trajectory(prob)

Return the initial state and control trajectories for the cartpole. 
The state trajectory rotates the pendulum from 0 to pi with constant 
velocity, and the cart remains still.
"""
function get_initial_trajectory(prob::CR3BPProblem)
    tf = prob.tf
    x0 = prob.x0
    xf = prob.xf
    times = range(0,tf,length=prob.N)
    X =[@SVector zeros(prob.n) for i=1:prob.N];
    for k = 1:prob.N
        Pos = Non_lin_exp(times[k]);
        X[k][1] = Pos[1];
        X[k][2] = Pos[2];
        X[k][3] = Pos[3];
        Vels = ForwardDiff.derivative(Non_lin_exp,times[k]);
        X[k][4] = Vels[1];
        X[k][5] = Vels[2];
        X[k][6] = Vels[3];
    end
    Xref_[end] = xf
    # X = [x0 + (xf - x0)*t for t in range(0,1, length=T)]
    U = [@SVector zeros(prob.m) for k = 1:prob.N];
    return X, U
end