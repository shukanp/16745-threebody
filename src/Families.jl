
function Generating_Orbit_Families(z0, numIters, scaleValue, Up)
    # Initial seed X0 and P0 (Xref and Uref since Uref is the timesteps h)
    j = numIters
    z_fam = deepcopy(z0)
    scale = 1
    x_fam = zeros(Nx, Nt, j)
    for i = 1:j
        # Run IPOPT
        prob = ProblemMOI(n_nlp, m_nlp)
        z_sol = solve(z_fam, prob)

        # Extract states and inputs from z_sol
        ztraj = reshape(z_sol, Nx + Nu, Nt);
        xtraj = ztraj[1:Nx, :];               # 6 x 91

        x_fam[:,:,i] = xtraj

        xmean = mean(xtraj,dims = 2)

        # Set input to IPOPT to the output from the previous solution
        z_fam = z_sol
        if Up
            scale += scaleValue
        else
            scale -= scaleValue
        end

        for k = 1:Nt
            Xref[k] = z0[x_idx[k]]
            Xref[k][1:3] = xmean[1:3] + scale*(Xref[k][1:3] - xmean[1:3])
        end
    end
    return x_fam
end
