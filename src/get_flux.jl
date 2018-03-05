function getflux(scheme,UL,UR,gamma,R,n,Area)

    gmi = gamma-1.0;
    
    # process left state
    rL = UL[1]
    uL = UL[2]/rL
    vL = UL[3]/rL
    wL = UL[4]/rL
    unL = uL*n[1] + vL*n[2] + wL*n[3]
    qL = sqrt(UL[2]^2 + UL[3]^2 + UL[4]^2)/rL
    pL = (gamma-1)*(UL[5] - 0.5*rL*qL^2)
    if ((pL<=0) || (rL<=0))
        println("Negative density or pressure")
    end
    rHL = UL[5] + pL
    HL = rHL/rL
    cL = sqrt(gamma*pL/rL)
    
    # left flux
    FL = UL # for allocation
    FL[1] = rL*unL
    FL[2] = UL[2]*unL + pL*n[1]
    FL[3] = UL[3]*unL + pL*n[2]
    FL[4] = UL[4]*unL + pL*n[3]
    FL[5] = rHL*unL
    
    # process right state
    rR = UR[1]
    uR = UR[2]/rR
    vR = UR[3]/rR
    wR = UR[4]/rR
    unR = uR*n[1] + vR*n[2] + wR*n[3]
    qR = sqrt(UR[2]^2 + UR[3]^2 + UR[4]^2)/rR
    pR = (gamma-1)*(UR[5] - 0.5*rR*qR^2)
    if ((pR<=0) || (rR<=0))
        println("Negative density or pressure")
    end
    rHR = UR[5] + pR
    HR = rHR/rR
    cR = sqrt(gamma*pR/rR)
    
    # right flux
    FR = UR # for allocation
    FR[1] = rR*unR
    FR[2] = UR[2]*unR + pR*n[1]
    FR[3] = UR[3]*unR + pR*n[2]
    FR[4] = UR[4]*unR + pR*n[3]
    FR[5] = rHR*unR;
    
    # difference in states
    du = UR - UL;
    
    # Roe average
    di     = sqrt(rR/rL)
    d1     = 1.0/(1.0+di)
    
    ui     = (di*uR + uL)*d1
    vi     = (di*vR + vL)*d1
    wi     = (di*wR + wL)*d1
    Hi     = (di*HR + HL)*d1
    
    af     = 0.5*(ui*ui+vi*vi+wi*wi);
    ucp    = ui*n[1] + vi*n[2] + wi*n[3]
    c2     = gmi*(Hi - af);
    if (c2<=0)
        println("Imaginary Wavespeed")
    end
    ci     = sqrt(c2);
    ci1    = 1.0/ci;
    
    # eigenvalues
    l = zeros(3,1);
    l[1] = ucp+ci;
    l[2] = ucp-ci;
    l[3] = ucp;
    
    # entropy fix
    epsilon = ci*0.1
    for i=1:3
        if ((l[i]<epsilon) && (l[i]>-epsilon))
            l[i] = 0.5*(epsilon + l[i]*l[i]/epsilon)
        end
    end
    
    l = abs.(l)
    l3 = l[3]
    
    # average and half-difference of 1st and 2nd eigs
    s1    = 0.5*(l[1] + l[2]);
    s2    = 0.5*(l[1] - l[2]);
    
    # left eigenvector product generators (see Theory guide)
    G1    = gmi*(af*du[1] - ui*du[2] - vi*du[3] - wi*du[4] + du[5]);
    G2    = -ucp*du[1]+du[2]*n[1]+du[3]*n[2]+du[4]*n[3]
    
    # required functions of G1 and G2 (again, see Theory guide)
    C1    = G1*(s1-l3)*ci1*ci1 + G2*s2*ci1;
    C2    = G1*s2*ci1          + G2*(s1-l3);
     
    # flux assembly
    F = FL  # for allocation
    F[1]    = 0.5*(FL[1]+FR[1])-0.5*(l3*du[1] + C1   )
    F[2]    = 0.5*(FL[2]+FR[2])-0.5*(l3*du[2] + C1*ui + C2*n[1])
    F[3]    = 0.5*(FL[3]+FR[3])-0.5*(l3*du[3] + C1*vi + C2*n[2])
    F[4]    = 0.5*(FL[4]+FR[4])-0.5*(l3*du[4] + C1*vi + C2*n[3])
    F[5]    = 0.5*(FL[5]+FR[5])-0.5*(l3*du[5] + C1*Hi + C2*ucp )
    F[6]    = 0.0

    return F

end
