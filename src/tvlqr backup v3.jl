#make sure you use dimensionless
Nu = 3

#Establish Variables
Q = Diagonal([10,10,10,10,10,10])
R = Diagonal([1,1,1]) #we want to conserve fuel so penalizing the controls takes priority
Qf = Diagonal([100,100,100,100,100,100]) #this can be the same as Q since its a periodic orbit so the end isnt special
K = [zeros(Nu,Nx)*NaN for k = 1:Nt-1]
P = [zeros(Nx,Nx)*NaN for k = 1:Nt]

#Initialize jacobians
A = [zeros(Nx,Nx) for k = 1:Nt]
B = [zeros(Nu,Nu) for k = 1:Nt-1]

#Initialize control
U = [zeros(Nu) for k = 1:Nt-1] #under ideal circumstances, no controls need to be expended. (does this mean B will stay zeros?)


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


function calc_gains!(xtraj, A, B, P, K, Qf, Q) #input here is the full state output of ipopt
    
    #Implement Riccati recursion for TVLQR

    P[end] .= Qf
    for k = reverse(1:Nt-1) #does the recursion only work with uniform timesteps???
        K[k] .= (R + B[k]'P[k+1]*B[k])\(B[k]'P[k+1]*A[k])
        P[k] .= Q + A[k]'P[k+1]*A[k] - A[k]'P[k+1]*B[k]*K[k]
    end

    return nothing
end

function get_control(X,k,U,K,xtraj)
    #       should return a vector of size (m,), where m is the number of controls
    u = zeros(3)
    
    u = U[k] - ctrl.K[k]*(x - ctrl.X[k])
    return u 
end


###############################
#ACTUALLY RUN THE TVLQR

#xtraj is our halo orbit trajectory we are trying to stick to
#Xsat is the satellite

linearize!(A,B,xtraj, U)
calc_gains!(xtraj, A, B, P, K, Qf, Q)
Xsat = [zeros(Nx) for k = 1:Nt] 

#i still need to add some minor purturbation to this or the controller should jst do nothing
#I'll either do it at the beginning or at some random timestep(s)
Xsat[1] = copy(xtraj[1])*1.0001
# Xsat[1] = Xsat[1] + [.0001, .0001, .0001, 0, 0, 0]

for k = 1:Nt-1
    U[k] = U[k] - K[k]*(Xsat[k]-xtraj[k]) #in hw this was often the get_control function
    Xsat[k+1] = rk4WithU(CR3BPdynamicsWithUforA, Xsat[k], U[k], utraj[k])  #discrete_dynamics(RK4, model, X[k], U[k], times[k], dt) 
end

###############################



#plot error
sse = zeros(length(Xsat))
for k = 1:Nt
    sse[k] = norm(Xsat[k]-xtraj[k])
end
t = Array(range(1, Nt, step = 1));

xtrajX = zeros(length(xtraj))
xtrajY = zeros(length(xtraj))
xtrajZ = zeros(length(xtraj))

XsatX = zeros(length(xtraj))
XsatY = zeros(length(xtraj))
XsatZ = zeros(length(xtraj))

for k in 1:Nt
    xtrajX[k] = xtraj[k][1]
    xtrajY[k] = xtraj[k][2]
    xtrajZ[k] = xtraj[k][3]
    
    XsatX[k] = Xsat[k][1]
    XsatY[k] = Xsat[k][2]
    XsatZ[k] = Xsat[k][3]
end


Ux = zeros(length(U))
Uy = zeros(length(U))
Uz = zeros(length(U))
for k in 1:Nt-1
    Ux[k] = U[k][1]
    Uy[k] = U[k][2]
    Uz[k] = U[k][3]
end

print(sse[49])
plot(t, sse)
# plot(t[1:end-1], Uz)
# plot(xtraj, Xsat)
# print(U)

##############################################

p = Plots.plot(xtrajX, xtrajY, label = "IPOPT")
Plots.plot!(p, XsatX, XsatY, label = "TVLQR")

###############################################
p = Plots.plot(t[1:end-1], Ux, label = "Ux")
Plots.plot!(p, t[1:end-1], Uy, label = "Uy")
Plots.plot!(p, t[1:end-1], Uz, label = "Uz")

###############################################

p = Plots.plot3d(XsatX,XsatY,XsatZ, label = "Xsat")
Plots.plot3d!(p, xtrajX,xtrajY,xtrajZ, label = "xtraj")
# plot3d(XsatX,XsatY,XsatZ)

#############################################