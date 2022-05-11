"The following file contains: 
-get_Period(LagrangePoint,System)-Calculate period for specified system
-scales(System)-Get scaling factors for the system of interest"

function get_Period(LagrangePoint, System)
    # Set parameters based on Lagrange point (need to implement parameters for other systems not Sun-Earth)
    if LagrangePoint == 1     # L1
        w_p = 2.086453455;
        k = 3.2292680962;
        c2 = 4.0610735668;
        c3 = 3.0200105081;
        c4 = 3.0305378797;
    elseif LagrangePoint == 2 # L2
        w_p = 2.05701;
        k = 3.18723;
        c2 = 3.94052;
        c3 = -2.97984;
        c4 = 2.97026;
    end

    # Get scales based on system (Sun-Earth for now)
    mu, LScale, VScale, TScale = scales(System)

    # Calculate constants s1 and s2
    λ = sqrt((c2 + sqrt(9*c2^2 - 8*c2))/2);
    d1 = ((3*λ^2)/k)*(k*(6*λ^2 - 1) - 2*λ); 
    a21 = 3*c3*(k^2 - 2)/(4*(1 + 2*c2)); 
    a22 = 3*c3/(4*(1 + 2*c2));
    a23 = -((3*c3*λ)/(4*k*d1))*(3*(k^3)*λ - 6*k*(k-λ) + 4); 
    a24 = -((3*c3*λ)/(4*k*d1))*(3*k*λ + 2);
    b21 = -3*c3*λ*(3*k*λ - 4)/(2*d1);
    b22 = -3*c3*λ/d1;
    d21 = -c3/(2*λ^2);
    s1 = (2*λ*(λ*(1 + k^2) - 2*k))^(-1)*((3/2)*c3*(2*a21*(k^2 - 2) - a23*(k^2 + 2) - 2*k*b21) - (3/8)*c4*(3*k^4 - 8*k^2 + 8));
    s2 = (2*λ*(λ*(1 + k^2) - 2*k))^(-1)*((3/2)*c3*(2*a22*(k^2 - 2) + a24*(k^2 + 2) + 2*k*b22 + 5*d21) + (3/8)*c4*(12 - k^2));
    
    # Calculate period based on w_p, s1, s2, Ax, Az, v
    Ax = 206000/LScale;
    Az = 110000/LScale;
    v1 = 0;
    v2 = s1*Ax^2 + s2*Az^2;
    v = 1 + v1 + v2;
    Tfinal = (2*pi)/(w_p*v);
    Tfinal_days = (Tfinal/(2*pi))*TScale/(24*60*60)
    return [Tfinal, Tfinal_days]
end


function scales(System)
    mu = 0
    LScale = 0
    VScale = 0
    TScale = 0
    if System == "Sun-Earth"
        msun = 1.98847e30;
        mmoon = 7.34767309e22;
        mearth = 5.9736e24 + mmoon;
        TOrbit = 1.53536256e7;        # L1 orbital period in seconds for halo orbits

        # Scaling factors
        mu = mearth/(msun + mearth);  # mu = 3.040423398444176e-6;
        LScale = 1.495978714e8;
        VScale = 29.784;
        TScale = 3.1536e7;
    elseif System == "Earth-Moon"
        mu = 1
        LScale = 1
        VScale = 1
        TScale = 1
    elseif System == "Sun-Jupiter"
        mu = 1
        LScale = 1
        VScale = 1
        TScale = 1
    end
    return [mu, LScale, VScale, TScale]
end
    