#Non linear Approximation to third order expansion
function Non_lin_exp(tau)

    k = 3.2292680962;
    c2 = 4.0610735668;
    c3 = 3.0200105081;
    c4 = 3.0305378797;

    lam = sqrt((c2 + sqrt(9*c2^2 - 8*c2))/2);
    d1 = ((3*lam^2)/k)*(k*(6*lam^2 - 1) - 2*lam); 
    d2 = ((8*lam^2)/k)*(k*(11*lam^2 - 1) - 2*lam);
    a21 = 3*c3*(k^2 - 2)/(4*(1 + 2*c2)); 
    a22 = 3*c3/(4*(1 + 2*c2));
    a23 = -((3*c3*lam)/(4*k*d1))*(3*(k^3)*lam - 6*k*(k-lam) + 4); 
    a24 = -((3*c3*lam)/(4*k*d1))*(3*k*lam + 2);
    b21 = -3*c3*lam*(3*k*lam - 4)/(2*d1);
    b22 = -3*c3*lam/d1;
    a31 = (-9*lam/(4*d2))*(4*c3*(k*a23 - b21) + k*c4*(4 + k^2))  + ((9*lam^2 + 1 - c2)/(2*d2))*(3*c3*(2*a23 - k*b21) + c4*(2 + 3*k^2));
    d21 = -c3/(2*lam^2);
    d31 = (3/(64*lam^2))*(4*c3*a24 + c4);
    d32 = (3/(64*lam^2))*(4*c3*(a23 - d21) + c4*(4 + k^2));
    a32 = (-9*lam/(4*d2))*(4*c3*(3*k*a24 - b22) + k*c4) - 3*((9*lam^2 + 1 - c2)/(2*d2))*(c3*(k*b22 - d21 - 2*a24) - c4);
    b31 = (3/(8*d2))*8*lam*(3*c3*(k*b21 - 2*a23) - c4*(2 + 3*k^2)) + (3/(8*d2))*((9*lam^2 + 1 + 2*c2)*(4*c3*(k*a23 - b21)) + k*c4*(4 + k^2));
    b32 = ((9*lam)/d2)*(c3*(k*b22 + d21 + -2*a24) - c4)  + 3*((9*lam^2 + 1 + 2*c2)/(8*d2))*(4*c3*(k*a24 - b22) + k*c4);

    xpoint = 0.932385;
    ypoint = 0;
    zpoint = 0;

    phi = 0;
    w_p = 2.086453455;
    # tau = range(0,Tfinal,length = Nt)
    # tau_1 = w_p.*tau .+ phi;

    # tau = range(0,Tfinal,length = Nt)
    tau_1 = w_p*tau + phi;

    scaling_factor = 1.495978714e8;
    Ax = 206000/scaling_factor;
    Ay = 665000/scaling_factor;
    Az = 110000/scaling_factor;
    m=1 
    dm = 2 - m;
    x1 =    -Ax*cos(tau_1) + (a23*Ax^2 - a24*Az^2)*cos(2*tau_1) + (a31*Ax^3 - a32*Ax*Az^2)*cos(3*tau_1) + a21*Ax^2 + a22*Az^2 + xpoint;
    y1 =     Ay*sin(tau_1) + (b21*Ax^2 - b22*Az^2)*sin(2*tau_1) + (b31*Ax^3 - b32*Ax*Az^2)*sin(3*tau_1) + ypoint;
    z1 =  dm*Az*cos(tau_1) + dm*d21*Ax*Az*cos(2*tau_1 - 3) + dm*(d32*Az*Ax^2 - d31*Az^3)*cos(3*tau_1) + zpoint;
    # x1 =    -Ax.*cos.(tau_1) .+ (a23*Ax^2 - a24*Az^2).*cos.(2 .* tau_1) .+ (a31*Ax^3 - a32*Ax*Az^2).*cos.(3 .* tau_1) .+ a21*Ax^2 .+ a22*Az^2 .+ xpoint;
    # y1 =     Ay.*sin.(tau_1) .+ (b21*Ax^2 - b22*Az^2).*sin.(2 .* tau_1) .+ (b31*Ax^3 - b32*Ax*Az^2).*sin.(3 .* tau_1) .+ ypoint;
    # z1 =  dm*Az.*cos.(tau_1) .+ dm*d21*Ax*Az.*cos.(2 .* tau_1 .- 3) .+ dm*(d32*Az*Ax^2 - d31*Az^3).*cos.(3 .* tau_1) .+ zpoint;
    return [x1, y1, z1]
end 
