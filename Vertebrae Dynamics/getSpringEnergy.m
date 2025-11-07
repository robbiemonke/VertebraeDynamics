function Es = getSpringEnergy(T, params)
    % Parameters
    P = params.P; 
    C1 = params.C1; 
    C2 = params.C2;
    ls_0 = params.ls_0; 
    lc_0 = params.lc_0;
    ls_max = params.ls_max; 
    lc_max = params.lc_max;
    k = params.k;

    % Solves equation for connection members vectors
    S = P * C1 - T * P * C2;

    % Calculation of spring energy
    Es = 0;
    for j = 1:12
        sj = S(1:3,j); 
        lj = norm(sj);

        if j == 1 || j == 4 || j == 9 || j == 12 % String is diagonal
            l0 = lc_0; 
        else % String is straight
            l0 = ls_0; 
        end

        if lj < l0 % String is slack
            e = 0;
        else % String is in tension
            e = 0.5 * k * (lj - l0) ^ 2;
        end
        Es = Es + e;
    end
end
