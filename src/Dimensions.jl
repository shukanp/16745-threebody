# Calculate period for specified system
function get_Period(LagrangePoint, Dimensionless, System)
    # Set parameters based on Lagrange point (need to implement parameters for other systems not Sun-Earth)
    if LagrangePoint == 1     # L1
        w_p = 2.086453455;
    elseif LagrangePoint == 2 # L2
        w_p = 2.05701;
    end

    # Get scales based on system (Sun-Earth for now)
    [mu, LScale, VScale, TScale] = scales(System)

    # Calculate period based on w_p, Ax, Az, v
    Ax = 206000/LScale;
    Az = 110000/LScale;
    v1 = 0;
    v2 = s1*Ax^2 + s2*Az^2;
    v = 1 + v1 + v2;
    Tfinal = (2*pi)/(w_p*v);
    if(Dimensionless)
        return Tfinal
    else
        return Tfinal/(2*pi))*TScale/(24*60*60)
    end
end

# Get scaling factors for the system of interest
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

    
    
    
    
    