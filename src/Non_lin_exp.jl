# Non linear Approximation to third order expansion
function Non_lin_exp(tau, LagrangePoint, Dimensionless, System)
    
    # Set parameters based on Lagrange point (need to implement parameters for other systems not Sun-Earth)
    x_L1 = 0.9899859900102089;
    x_L2 = 1.0100751933379162;
    if LagrangePoint == 1     # L1
        x_L = x_L1;
        y_L = 0;
        z_L = 0;
        k = 3.2292680962;
        c2 = 4.0610735668;
        c3 = 3.0200105081;
        c4 = 3.0305378797;
        w_p = 2.086453455;
    elseif LagrangePoint == 2 # L2
        x_L = x_L2;
        y_L = 0;
        z_L = 0;
        k = 3.18723;
        c2 = 3.94052;
        c3 = -2.97984;
        c4 = 2.97026;
        w_p = 2.05701;
    end

    # Calculate all necessary constants based on k, c2, c3, c4
    λ = sqrt((c2 + sqrt(9*c2^2 - 8*c2))/2);
    d1 = ((3*λ^2)/k)*(k*(6*λ^2 - 1) - 2*λ); 
    d2 = ((8*λ^2)/k)*(k*(11*λ^2 - 1) - 2*λ);
    a21 = 3*c3*(k^2 - 2)/(4*(1 + 2*c2)); 
    a22 = 3*c3/(4*(1 + 2*c2));
    a23 = -((3*c3*λ)/(4*k*d1))*(3*(k^3)*λ - 6*k*(k-λ) + 4); 
    a24 = -((3*c3*λ)/(4*k*d1))*(3*k*λ + 2);
    b21 = -3*c3*λ*(3*k*λ - 4)/(2*d1);
    b22 = -3*c3*λ/d1;
    a31 = (-9*λ/(4*d2))*(4*c3*(k*a23 - b21) + k*c4*(4 + k^2))  + ((9*λ^2 + 1 - c2)/(2*d2))*(3*c3*(2*a23 - k*b21) + c4*(2 + 3*k^2));
    d21 = -c3/(2*λ^2);
    d31 = (3/(64*λ^2))*(4*c3*a24 + c4);
    d32 = (3/(64*λ^2))*(4*c3*(a23 - d21) + c4*(4 + k^2));
    a32 = (-9*λ/(4*d2))*(4*c3*(3*k*a24 - b22) + k*c4) - 3*((9*λ^2 + 1 - c2)/(2*d2))*(c3*(k*b22 - d21 - 2*a24) - c4);
    b31 = (3/(8*d2))*8*λ*(3*c3*(k*b21 - 2*a23) - c4*(2 + 3*k^2)) + (3/(8*d2))*((9*λ^2 + 1 + 2*c2)*(4*c3*(k*a23 - b21)) + k*c4*(4 + k^2));
    b32 = ((9*λ)/d2)*(c3*(k*b22 + d21 + -2*a24) - c4)  + 3*((9*λ^2 + 1 + 2*c2)/(8*d2))*(4*c3*(k*a24 - b22) + k*c4);

    # Get scales based on system (Sun-Earth for now)
    [mu, LScale, VScale, TScale] = scales(System)

    # Set x, y, z of Lagrange Point based on if we want dimensionless system or not
    if(Dimensionless)
        xpoint = x_L;
        ypoint = y_L;
        zpoint = z_L;
    else
        xpoint = x_L*LScale;
        ypoint = y_L*LScale;
        zpoint = z_L*LScale;
    end

    phi = 0;
    tau_1 = w_p*tau + phi;

    Ax = 206000/LScale;
    Ay = 665000/LScale;
    Az = 110000/LScale;
    
    m = 1
    dm = 2 - m;
    x1 =   -Ax*cos(tau_1) + (a23*Ax^2 - a24*Az^2)*cos(2*tau_1) + (a31*Ax^3 - a32*Ax*Az^2)*cos(3*tau_1) + a21*Ax^2 + a22*Az^2;
    y1 =    Ay*sin(tau_1) + (b21*Ax^2 - b22*Az^2)*sin(2*tau_1) + (b31*Ax^3 - b32*Ax*Az^2)*sin(3*tau_1);
    z1 = dm*Az*cos(tau_1) + dm*d21*Ax*Az*cos(2*tau_1 - 3) + dm*(d32*Az*Ax^2 - d31*Az^3)*cos(3*tau_1);
    
    # Adjust x, y, z by actual location of Lagrange Point
    if(Dimensionless)
        x1 = x1 + xpoint;
        y1 = y1 + ypoint;
        z1 = z1 + zpoint;
    else
        x1 = x1*LScale + xpoint;
        y1 = y1*LScale + ypoint;
        z1 = z1*LScale + zpoint;
    end
    
    return [x1, y1, z1]
end 

# Linear Approximation of the CR3BP
# ------> NEED TO UPDATE LINEAR APPROX SAME WAY AS 3RD ORDER ABOVE
function lin_approx(tau)
    omega = 2.0152105515;
    scaling_factor = 1.495978714e8;
    w_p = 2.086453455;
    w_v = 2.0152105515;
    Ax = 206000/scaling_factor;
    Ay = 665000/scaling_factor;
    Az = 110000/scaling_factor;
    xpoint = 0.932385;
    ypoint = 0;
    zpoint = 0;
    mu = -9.537e-4;
    phi = 0;
    m = 1;
    psi = m*(pi/2) + phi;
    x1 = -Ax*cos(tau + phi) + xpoint;
    y1 = Ay*sin(tau + phi) + ypoint;
    z1 = Az*sin(tau + 1*(pi/2) + phi) + zpoint;
    return [x1, y1, z1]
end 
