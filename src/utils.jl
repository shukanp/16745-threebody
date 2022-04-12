using LinearAlgebra

##############         Helper Functions         ##############
function skew(v)
    [0 -v[3] v[2]; v[3] 0 -v[1]; -v[2] v[1] 0]
end

function matrix_isapprox(A,B;atol = 1e-8)
    if (size(A) == size(B)) && isapprox(vec(A),vec(B),atol = atol)
        return true
    else
        return false
    end
end
##############################################################




##############   Controller Structs/Functions   ##############  (HW2, Q1.ipynb, Q3.ipynb)
struct LQRController
    K::Matrix{Float64}             # feedback gains ((m,n),N-1)
    Xref::Vector{Vector{Float64}}  # reference states
    Uref::Vector{Vector{Float64}}  # reference inputs
    xeq::Vector{Float64}           # equilibrium states
    ueq::Vector{Float64}           # equilibrium controls
    times::Vector{Float64}         # times          (N,)
    infinite_horizon::Bool         # use infinite horizon control
end

struct MPCController{S}
    P::SparseMatrixCSC{Float64,Int}
    q::Vector{Float64}
    A::SparseMatrixCSC{Float64,Int}
    lb::Vector{Float64}
    ub::Vector{Float64}
    Nmpc::Int
    solver::S
    Xref::Vector{Vector{Float64}}
    Uref::Vector{Vector{Float64}}
    times::Vector{Float64}
end
##############################################################




##############          Model Functions         ##############  (HW1, quadruped.jl; HW2, quadruped.jl)
RobotDynamics??
RigidBodyDynamics??

"""
    UnitreeA1 <: AbstractModel
A model of the UnitreeA1  ???
"""
struct UnitreeA1{C} <: AbstractModel
end

function UnitreeA1()
    UnitreeA1(??)
end
##############################################################




##############      Equilibrium Functions       ##############  (HW1, Q2.ipynb; HW2, quadruped.jl)
"""
    kkt_conditions(model, x, u, λ, A, B)
"""
function kkt_conditions()
end

"""
    kkt_jacobian(model, x, u, λ, A, B, [ρ])
"""
function kkt_jacobian()
end

"""
    newton_solve(model, x_guess, [mvis; verbose])
"""
function newton_solve()
end
##############################################################




##############        Simulate Functions        ##############  (HW2, rocket.jl)
"""
    simulate(model, x0, ctrl; [kwargs...])
"""
function simulate()
end
##############################################################




##############       Visualize Functions        ##############
"""
    set_mesh!(vis, model::BicycleModel)    
Add a simple model of the ??? model to the MeshCat visualizer `vis`
"""
function set_mesh!(vis, model::??)
end

"""
    visualize!(vis, model, x)
Visualize a single frame of `model` given the state vector `x`.
"""
function visualize!(vis, model::??, x::StaticVector)
end

"""
    initialize_visualizer(model)
Launch a MeshCat visualizer and insert the geometry for `model`.
"""
function initialize_visualizer(model::??)
end
##############################################################



