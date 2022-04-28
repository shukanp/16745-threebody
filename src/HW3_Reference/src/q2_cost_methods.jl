# TASK: Finish the following methods
#       eval_f
#       grad_f!
#       hess_f!

"""
    eval_f(nlp, Z)

Evaluate the objective, returning a scalar. The continuous time objective is of the form:

```math
\\int_{t0}^{tf} \\ell(x(t), u(t)) dt 
```

You need to approximate this with an integral of the form:
```math
\\sum_{k=1}^{N-1} \\frac{h}{6}(\\ell(x_k,u_k) + 4\\ell(x_m, u_m) + \\ell(x_{k+1}, u_{k+1}))
```

where
```math
x_m = \\frac{1}{2} (x_1 + x_2) + \\frac{h}{8}(f(x_1, u_1, t) - f(x_2, u_2, t + h))
```
and
```math
u_m = \\frac{1}{2} (u_1 + u_2)
```
"""
function eval_f(nlp::NLP, Z) 
    evaluate_dynamics!(nlp, Z)
    evaluate_midpoints!(nlp, Z)
    # TASK: compute the objective value (cost)
    J = NaN
    J = 0.0
    
    ix,iu = nlp.xinds, nlp.uinds

    Jstage = nlp.stagecost  # Stage cost object
    Jterm = nlp.termcost    # Term cost object

    for k = 1:nlp.N - 1
        x1, x2 = Z[ix[k]], Z[ix[k + 1]]
        u1, u2 = Z[iu[k]], Z[iu[k + 1]]
        xm = nlp.xm[k]
        um = nlp.um[k]
        t = nlp.times[k]
        h = nlp.times[k + 1] - nlp.times[k]
        
        J += (h/6)*(stagecost(Jstage, x1, u1) + 4*stagecost(Jstage, xm, um) + stagecost(Jstage, x2, u2))
    end
    J += termcost(Jterm, Z[ix[end]])

    return J
end

"""
    grad_f!(nlp, grad, Z)

Evaluate the gradient of the objective at `Z`, storing the result in `grad`.
"""
function grad_f!(nlp::NLP{n,m}, grad, Z) where {n,m}
    evaluate_dynamics!(nlp, Z)
    evaluate_midpoints!(nlp, Z)
    evaluate_dynamics_jacobians!(nlp, Z)
    ix,iu = nlp.xinds, nlp.uinds
    Jstage = nlp.stagecost
    Jterm = nlp.termcost
    grad .= 0
    Jx1JL = 0.0    # Current cost (k)
    Jx2JL = 0.0    # Midpoint cost (k + 1)/2
    Ju1JL = 0.0    # Current cost (k)
    Ju2JL = 0.0    # Midpoint cost (k + 1)/2
    Q = Jstage.Q
    q = Jstage.q
    R = Jstage.R
    r = Jstage.r
    for k = 1:nlp.N - 1
        x1,x2 = Z[ix[k]], Z[ix[k + 1]]
        u1,u2 = Z[iu[k]], Z[iu[k + 1]]
        xm = nlp.xm[k]
        um = nlp.um[k]
        t = nlp.times[k]
        h = nlp.times[k + 1] - nlp.times[k]
        
        # TASK: Compute the cost gradient
        # Cost gradient for states
        ∂lx1∂x1 = x1'*Q + q'
        ∂lx2∂x2 = x2'*Q + q'
        ∂lxm∂xm = xm'*Q + q'
        ∂xm∂x1 =  (h/8)*nlp.A[k] + I/2
        ∂xm∂x2 =  -(h/8)*nlp.A[k + 1] + I/2
        Jx1JL = (h/6)*(∂lx1∂x1 + 4*∂lxm∂xm*∂xm∂x1)'
        Jx2JL = (h/6)*(∂lx2∂x2 + 4*∂lxm∂xm*∂xm∂x2)'

        # Cost gradient for inputs
        ∂lu1∂u1 = u1'*R + r'
        ∂lu2∂u2 = u2'*R + r'
        ∂lum∂um = um'*R + r'
        ∂lum∂u1 = (xm'*Q + q')*((h/8)*nlp.B[k])
        ∂lum∂u2 = -(xm'*Q + q')*((h/8)*nlp.B[k + 1])
        Ju1JL = (h/6)*(∂lu1∂u1 + 4*∂lum∂u1 + ∂lum∂um*(I/2))'
        Ju2JL = (h/6)*(∂lu2∂u2 + 4*∂lum∂u2 + ∂lum∂um*(I/2))'
        grad[ix[k]]     += Jx1JL
        grad[ix[k + 1]] += Jx2JL
        grad[iu[k]]     += Ju1JL
        grad[iu[k + 1]] += Ju2JL
    end
    grad[ix[nlp.N]] += Jterm.Q*Z[ix[nlp.N]] + Jterm.q
    return nothing
end

