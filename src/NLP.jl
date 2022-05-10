#Our code 
function primal_bounds(n)
    x_l = -Inf*ones(n)
    x_u = Inf*ones(n) #+_ 20% dt_ref
    for k = 1:Nt
        x_l[u_idx[k]] .= h*0.01 #Something reasonable based on dt_ref 0.8*dt
        x_u[u_idx[k]] .= h*1.5
    end 
    return x_l, x_u
end
#Our code 

function constraint_bounds(m; idx_ineq=(1:0))
    c_l = zeros(m)
    c_u = zeros(m)
    return c_l, c_u
end

function row_col!(row,col,r,c)
    for cc in c
        for rr in r
            push!(row,convert(Int,rr))
            push!(col,convert(Int,cc))
        end
    end
    return row, col
end

####################-SPARSITY-####################
function initialize_sparsity!(prob::MOI.AbstractNLPEvaluator)
    # n,m = Nx, Nu
    n = 0 
    blocks = prob.blocks
    idx = x_idx[1]
    for k = 1:Nt-1
        idx = idx .+ n
        setblock!(blocks, idx, x_idx[k])
        setblock!(blocks, idx, u_idx[k])
        setblock!(blocks, idx, x_idx[k+1])
        n = Nx
    end
#Periodicity Constraint
    idx = idx .+ n 
    for i = 1:n
        setblock!(blocks, idx[i], x_idx[1][i])
        setblock!(blocks, idx[i], x_idx[end][i])
    end
end
####################-SPARSITY-####################

function sparsity_jacobian(n,m)
#######-Dense Jac-###############
    row = []
    col = []
    r = 1:m
    c = 1:n
    row_col!(row,col,r,c)
    return collect(zip(row,col))
end 

function jac_c!(prob::MOI.AbstractNLPEvaluator, jacvec::AbstractVector, Z) where {n,m}

    
    jac = NonzerosVector(jacvec, prob.blocks)

    X_ = [zeros(Nx) for k = 1:Nt]
    h_ = [zeros(Nu) for k = 1:Nt]

    for k=1:Nt
        X_[k] = Z[x_idx[k]]
        h_[k] .= Z[u_idx[k]]
    end
    
    n_ = 0
    idx = x_idx[1] 

    for k = 1:Nt-1

        idx = idx .+ n_
        x1,x2 = Z[x_idx[k]], Z[x_idx[k+1]]
        u1,u2 = Z[u_idx[k]], Z[u_idx[k+1]]

        jac_x1 = view(jac, idx, x_idx[k])
        jac_u1 = view(jac, idx, u_idx[k])
        jac_x2 = view(jac, idx, x_idx[k+1])
        
        f1,f2 = CR3BPdynamics(x1), CR3BPdynamics(x2)

        A1 = ForwardDiff.jacobian(CR3BPdynamics,x1)
        A2 = ForwardDiff.jacobian(CR3BPdynamics,x2)
        

        xm = 0.5.*(x1 + x2) .+ (h/8.0).*(f1 - f2) 

        Am = ForwardDiff.jacobian(CR3BPdynamics,xm)

        dx1 = Am*(0.5*I + (h_[k][1]/8)*A1) +(3/(2*h_[k][1]))*I + 0.25*A1 
        du1 = Am*((1/8)*(f1-f2)) - (3/(2*h_[k][1]^2))*(x1-x2)
        dx2 =Am*(0.5*I + (h_[k][1]/8)*(-A2)) - (3/(2*h_[k][1]))*I + 0.25*A2

        jac_x1 .= dx1
        jac_u1 .= du1
        jac_x2 .= dx2
        n_ = Nx
    end
    idx = idx .+ n_
    for i = 1:Nx
        jac[idx[i], x_idx[1][i]] = -1
        jac[idx[i], x_idx[end][i]] = 1
    end
end

function sparsity_hessian(n,m)

    row = []
    col = []

    r = 1:m
    c = 1:n

    row_col!(row,col,r,c)

    return collect(zip(row,col))
end

function dynamics_constraint!(c,ztraj)
    d = reshape(c,Nx,Nt-1) 
    z = reshape(ztraj,Nx+Nu,Nt) #ztraj[7,50]
    for k = 1:(Nt-1)
        x1 = z[1:Nx,k]
        u1 = z[(Nx+1):(Nx+Nu),k]
        x2 = z[1:Nx,k+1]
        d[:,k] = dircol_dynamics(x1,x2,u1)
    end
    return nothing
end

function con!(c,ztraj)
    z = reshape(ztraj,Nx+Nu,Nt)#7 x 50
    @views dynamics_constraint!(c[1:(end-Nx)],ztraj) #6 x 49 = 294
    #Periodicity Constraint 
    c[(end-Nx+1):end] .= z[1:Nx,end] - z[1:Nx,1]
end