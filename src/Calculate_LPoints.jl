# Calculate Lagrange Points for specified system
function get_LPoints(System)
    mu, LScale, VScale, TScale = scales(System)

    L1 = 0
    L2 = 0
    L3 = 0
    L4 = 0
    L5 = 0
    
    if System == "Sun-Earth"
        msun = 1.98847e30;
        mmoon = 7.34767309e22;
        mearth = 5.9736e24 + mmoon;

        L1 = (1 - (mu/3)^(1/3), 0)
        L2 = (1 + (mu/3)^(1/3), 0)
        L3 = (1 + (5/12)*mu, 0)
        L4 = ((msun - mu)/2, sqrt(3)/2)
        L5 = ((msun - mu)/2, -sqrt(3)/2)
    end
    return [L1, L2, L3, L4, L5]
end

    