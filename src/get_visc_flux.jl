function getViscFlux(scheme,UL,UR,gUL,gUR,gamma,Pr,mu,n)
    
    flux    =   zeros(size(UL))

    if scheme=="Dirty"
        τ11 =   (4.0/3.0)*(gUL[1,1]+gUR[1,1])*0.5
        τ22 =   (4.0/3.0)*(gUL[2,2]+gUR[2,2])*0.5
        τ33 =   (4.0/3.0)*(gUL[3,3]+gUR[3,3])*0.5
        
        τ12 =   (gUL[1,2]+gUR[1,2]+gUL[2,1]+gUR[2,1])*0.5
        τ23 =   (gUL[2,3]+gUR[2,3]+gUL[3,2]+gUR[3,2])*0.5
        τ31 =   (gUL[1,3]+gUR[1,3]+gUL[3,1]+gUR[3,1])*0.5

        H_x =   (gUL[4,1]+gUR[4,1])*0.5
        H_y =   (gUL[4,2]+gUR[4,2])*0.5
        H_z =   (gUL[4,3]+gUR[4,3])*0.5

        u   =   (UL[2]+UR[2])*0.5
        v   =   (UL[3]+UR[3])*0.5
        w   =   (UL[4]+UR[4])*0.5

        Ef1 =   u*τ11+v*τ12+w*τ31 + (mu/Pr)*H_x
        Ef2 =   u*τ12+v*τ22+w*τ23 + (mu/Pr)*H_y
        Ef3 =   u*τ31+v*τ23+w*τ33 + (mu/Pr)*H_z

        F   =   [0.0            0.0             0.0
                 τ11            τ12             τ31
                 τ12            τ22             τ23
                 τ31            τ23             τ33
                 Ef1            Ef2             Ef3
                 0.0            0.0             0.0]

        flux    =   F*n
    end

    return flux

end
