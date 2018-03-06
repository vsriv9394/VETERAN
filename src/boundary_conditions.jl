function setbcflux(scheme,stateL,stateR,gamma,R,surf,ext_side)

    flux    =   [0.0,0.0,0.0,0.0,0.0,0.0]
    n       =   [0.0,0.0,0.0]
    Area    =   sqrt(surf[1]^2+surf[2]^2+surf[3]^2)
    n[1]    =   surf[1]/Area
    n[2]    =   surf[2]/Area
    n[3]    =   surf[3]/Area
    bc_type =   surf[4]

    state   =   [0.0,0.0,0.0,0.0,0.0,0.0]
    stateb  =   [0.0,0.0,0.0,0.0,0.0,0.0]

    # Full State
    if bc_type==0
        
        flux   =   getflux(scheme,stateL,stateR,gamma,R,n,Area)
    
    # Inviscid Wall
    elseif bc_type==1

        if ext_side==0
            state = stateR
        else
            state = stateL
        end

        stateb[1] = state[1]

        ρv      =   [0.0,0.0,0.0]

        ρv[1]   =   state[2]
        ρv[2]   =   state[3]
        ρv[3]   =   state[4]

        ρv_n    =   (ρv'*n)*n

        stateb[2] = ρv[1]-2*ρv_n[1]
        stateb[3] = ρv[2]-2*ρv_n[2]
        stateb[4] = ρv[3]-2*ρv_n[3]

        stateb[5] = state[5]

        if ext_side==0
            flux   =   getflux(scheme,stateb,state,gamma,R,n,Area)
        else
            flux   =   getflux(scheme,state,stateb,gamma,R,n,Area)
        end

    # Inflow
    elseif bc_type==2
    
        # stateb    c_t, n_in_x, n_in_y, n_in_z, p_t
        if ext_side==0
            state   =   stateR
            stateb  =   stateL
        else
            state   =   stateL
            stateb  =   stateR
        end
        
        c_t     =   stateb[1]

        v       =   [0.0,0.0,0.0]
        v[1]    =   state[2]/state[1]
        v[2]    =   state[3]/state[1]
        v[3]    =   state[4]/state[1]

        v_nm    =   v'*n

        ke      =   0.5*state[1]*(v[1]^2+v[2]^2+v[3]^2)

        c       =   sqrt(gamma*(gamma-1)*(state[5]-ke)/state[1])

        J       =   v_nm+2*c/(gamma-1)
        d_n     =   stateb[2]*n[1]+stateb[3]*n[2]+stateb[4]*n[3]

        aaa     =   c_t^2*d_n^2-0.5*(gamma-1)*J^2
        bbb     =   4*c_t^2*d_n/(gamma-1)
        ccc     =   4*c_t^2/(gamma-1)^2-J^2

        M_b     =   (sqrt(bbb*bbb-4*aaa*ccc)-bbb)/(2*aaa)
        c_b     =   c_t/sqrt(1+0.5*(gamma-1)M_b^2)

        p_t     =   stateb[5]
        p_b     =   p_t*(c_b/c_t)^(2*gamma/(gamma-1))
        ρ_b     =   gamma*p_b/c_b^2
        v_b     =   M_b*c_b*[stateb[2],stateb[3],stateb[4]]
        ρeb     =   p_b/(gamma-1)+0.5*ρ_b*(v_b[1]^2+v_b[2]^2+v_b[3]^2)

        F       =   [ρ_b*v_b[1]             ρ_b*v_b[2]              ρ_b*v_b[3]
                     ρ_b*v_b[1]*v_b[1]+p_b  ρ_b*v_b[1]*v_b[2]       ρ_b*v_b[1]*v_b[3]
                     ρ_b*v_b[2]*v_b[1]      ρ_b*v_b[2]*v_b[2]+p_b   ρ_b*v_b[2]*v_b[3]
                     ρ_b*v_b[3]*v_b[1]      ρ_b*v_b[3]*v_b[2]       ρ_b*v_b[3]*v_b[3]+p_b
                     ρeb*v_b[1]+v_b[1]*p_b  ρeb*v_b[2]+v_b[2]*p_b   ρ_b*v_b[3]+v_b[3]*p_b
                     0.0                    0.0                     0.0                  ]

        flux    =   F*n

    # Subsonic Outflow
    elseif bc_type==3
        
        if ext_side==0
            state   =   stateR
            stateb  =   stateL
        else
            state   =   stateL
            stateb  =   stateR
        end

        p       =   (gamma-1)*(state[5]-0.5*(state[2]^2+state[3]^2+state[4]^2)/state[1])
        p_b     =   (gamma-1)*(stateb[5]-0.5*(stateb[2]^2+stateb[3]^2+stateb[4]^2)/stateb[1])
        S       =   p/state[1]^gamma
        ρ_b     =   (p_b/S)^(1/gamma)
        c_b     =   sqrt(gamma*p_b/ρ_b)
        c       =   sqrt(gamma*p/state[1])
        
        v       =   [0.0,0.0,0.0]
        v[1]    =   state[2]/state[1]
        v[2]    =   state[3]/state[1]
        v[3]    =   state[4]/state[1]

        v_nm    =   v'*n
        v_t     =   v-v_nm*n

        v_bnm   =   v_nm+2*(c-c_b)/(gamma-1)
        v_b     =   v_t+v_bnm*n
        ρeb     =   p_b/(gamma-1)+0.5*ρ_b*(v_b[1]^2+v_b[2]^2+v_b[3]^2)

        F       =   [ρ_b*v_b[1]             ρ_b*v_b[2]              ρ_b*v_b[3]
                     ρ_b*v_b[1]*v_b[1]+p_b  ρ_b*v_b[1]*v_b[2]       ρ_b*v_b[1]*v_b[3]
                     ρ_b*v_b[2]*v_b[1]      ρ_b*v_b[2]*v_b[2]+p_b   ρ_b*v_b[2]*v_b[3]
                     ρ_b*v_b[3]*v_b[1]      ρ_b*v_b[3]*v_b[2]       ρ_b*v_b[3]*v_b[3]+p_b
                     ρeb*v_b[1]+v_b[1]*p_b  ρeb*v_b[2]+v_b[2]*p_b   ρ_b*v_b[3]+v_b[3]*p_b
                     0.0                    0.0                     0.0                  ]

        flux    =   F*n
    
    end

    return flux

end
