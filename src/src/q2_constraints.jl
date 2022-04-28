
# TASK: Complete the following methods
#       eval_c!
#       jac_c!
"""
    eval_c!(nlp, c, Z)

Evaluate the equality constraints at `Z`, storing the result in `c`.
The constraints should be ordered as follows: 
1. Initial condition ``x_1 = x_\\text{init}``
2. Hermite Simpson Dynamics: ``\\frac{h}{6} (f(x_k, u_k) + 4 f(x_m, u_m) + f(x_{k+1}, u_{k+1})) + x_k - x_{k+1} = 0``
3. Terminal constraint ``x_N = x_\\text{goal}``

Consider leveraging the caches in `nlp` to evaluate the dynamics and the midpoints 
before the main loop, so that you can making redundant calls to the dynamics.

Remember, you will loose points if you make more dynamics calls than necessary. 
Start with something that works, then think about how to eliminate any redundant 
dynamics calls.
"""
function eval_c!(nlp::NLP{n,m}, c, Z) where {n,m}
    evaluate_dynamics!(nlp, Z)
    evaluate_midpoints!(nlp, Z)

    N = nlp.N
    xi,ui = nlp.xinds, nlp.uinds
    idx = xi[1]

    # TODO: initial condition
    c .= 0
    c[idx] = Z[idx] - nlp.x0
    
    # Dynamics
    for k = 1:N-1
        idx = idx .+ n
        x1,x2 = Z[xi[k]], Z[xi[k+1]]
        u1,u2 = Z[ui[k]], Z[ui[k+1]]
        t = nlp.times[k]
        h = nlp.times[k+1] - nlp.times[k]
        
        # TASK: Dynamics constraint
        c[idx] += (h/6)*(nlp.f[k] + 4*nlp.fm[k] + nlp.f[k + 1]) + x1 - x2
    end

    # TODO: terminal constraint
    idx = idx .+ n
    c[idx] .= Z[xi[end]] - nlp.xf
    
    return c
end

"""
    jac_c!(nlp, jac, Z)

Evaluate the constraint Jacobian, storing the result in the matrix `jac`.
You will need to apply the chain rule to calculate the derivative of the dynamics
constraints with respect to the states and controls at the current and next time 
steps.

### Use of automated differentiation tools
You are not allowed to use automatic differentiation methods for this function. 
You are only allowed to call `dynamics_jacobians` (which under the hood does use
ForwardDiff). You are allowed to check your answer with these tools, but your final 
solution should not use them.
"""
function jac_c!(nlp::NLP{n,m}, jac, Z) where {n,m}
    evaluate_dynamics_jacobians!(nlp, Z)
    evaluate_midpoint_jacobians!(nlp, Z)

    # TODO: Initial condition
    xi, ui = nlp.xinds, nlp.uinds
    idx = xi[1]
    jac .= 0
    jac[idx, idx] .= I(n)

    for k = 1:nlp.N - 1
        idx = idx .+ n
        x1, x2 = Z[xi[k]], Z[xi[k+1]]
        u1, u2 = Z[ui[k]], Z[ui[k+1]]
        t = nlp.times[k]
        h = nlp.times[k+1] - nlp.times[k]
        
        jac_x1 = view(jac, idx, xi[k])
        jac_x1 .= (h/6)*(4*nlp.Am[k, 1]*nlp.Am[k, 2] + nlp.A[k]) + I
        jac_x2 = view(jac, idx, xi[k + 1])
        jac_x2 .= (h/6)*(4*nlp.Am[k, 1]*nlp.Am[k, 3] + nlp.A[k + 1]) - I

        jac_u1 = view(jac, idx, ui[k])
        jac_u1 .= (h/6)*(4*(nlp.Bm[k, 1]*(1/2) + nlp.Am[k, 1]*nlp.Bm[k, 2]) + nlp.B[k])
        jac_u2 = view(jac, idx, ui[k + 1])
        jac_u2 .= (h/6)*(4*(nlp.Bm[k, 1]*(1/2) + nlp.Am[k, 1]*nlp.Bm[k, 3]) + nlp.B[k + 1])
    end        
    # TODO: Terminal constraint
    idx = idx .+ n 
    jac[idx, nlp.xinds[nlp.N]] .= I(n)
end






# EXTRA CREDIT: Specify the sparsity directly in the nonzeros vector
#               Read the MathOptInterface documentation!
function jac_c!(nlp::NLP{n,m}, jacvec::AbstractVector, Z) where {n,m}
end
