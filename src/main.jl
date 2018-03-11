function simulate(Files,Schemes,Time_specs,Gas_prop)

    #---------------------------------------------------------------------------
    # Initializing Arguments
    #---------------------------------------------------------------------------
    Meshfile    =   Files[1]
    Restartfile =   Files[2]
    wfile       =   Files[3]

    conv_scheme =   Schemes[1]
    visc_scheme =   Schemes[2]
    time_scheme =   Schemes[3]

    dt          =   Time_specs[1]
    stop        =   Time_specs[2]
    w_itr       =   Time_specs[3]

    gamma       =   Gas_prop[1]
    Pr          =   Gas_prop[2]
    mu          =   Gas_prop[3]

    ############################################################################
    # Setup the geometry
    ############################################################################
    Cells,SurfX,SurfY,SurfZ,X,Y,Z = create_Geometry(Meshfile,Restartfile)

    ############################################################################
    # Write the initial condition file
    ############################################################################
    println("Writing file for time-step ",0)
    write(string(wfile,lpad(0,6,0),".vtrn"),Cells)

    conv_Flow   =   zeros(size(Cells))
   
    ############################################################################
    # Initialize iterations
    ############################################################################
    it = 1

    while it*dt<=stop

        println("=================================================================================================")
        println("Iteration ",it," initialized")

        #------------------------------------------------------------------------
        # Time Integration schemes
        #------------------------------------------------------------------------
        if time_scheme=="Euler_explicit"

            setupbc(Meshfile,Cells,SurfX,SurfY,SurfZ)
            conv_Flow   =   flow_eval(conv_scheme,visc_scheme,Cells,SurfX,SurfY,SurfZ,gamma,Pr,mu)
            Cells       =   Cells   +   dt*conv_Flow
        
        elseif time_scheme=="RK2"
        
            setupbc(Meshfile,Cells,SurfX,SurfY,SurfZ)
            conv_Flow1  =   flow_eval(conv_scheme,visc_scheme,Cells,SurfX,SurfY,SurfZ,gamma,Pr,mu)
            Cells1      =   Cells   +   dt*conv_Flow1
            setupbc(Meshfile,Cells1,SurfX,SurfY,SurfZ)
            conv_Flow2  =   flow_eval(conv_scheme,visc_scheme,Cells1,SurfX,SurfY,SurfZ,gamma,Pr,mu)
            Cells       =   Cells   +   dt*(conv_Flow1+conv_Flow2)*0.5
        
        elseif time_scheme=="RK4"
            
            setupbc(Meshfile,Cells,SurfX,SurfY,SurfZ)
            conv_Flow1  =   flow_eval(conv_scheme,visc_scheme,Cells,SurfX,SurfY,SurfZ,gamma,Pr,mu)
            Cells1      =   Cells   +   0.5*dt*conv_Flow1
            setupbc(Meshfile,Cells1,SurfX,SurfY,SurfZ)
            conv_Flow2  =   flow_eval(conv_scheme,visc_scheme,Cells1,SurfX,SurfY,SurfZ,gamma,Pr,mu)
            Cells1      =   Cells   +   0.5*dt*conv_Flow2
            setupbc(Meshfile,Cells1,SurfX,SurfY,SurfZ)
            conv_Flow3  =   flow_eval(conv_scheme,visc_scheme,Cells1,SurfX,SurfY,SurfZ,gamma,Pr,mu)
            Cells1      =   Cells   +   dt*conv_Flow3
            setupbc(Meshfile,Cells1,SurfX,SurfY,SurfZ)
            conv_Flow4  =   flow_eval(conv_scheme,visc_scheme,Cells1,SurfX,SurfY,SurfZ,gamma,Pr,mu)
            Cells       =   Cells   +   dt/6*(conv_Flow1+2*conv_Flow2+2*conv_Flow3+conv_Flow2)
        
        elseif time_scheme=="LSRK4"

            aRK   =   [0.000000000000000, -0.417890474499852, -1.192151694642677, -1.697784692471528, -1.514183444257156]
            bRK   =   [0.149659021999229, 0.379210312999627, 0.822955029386982, 0.699450455949122, 0.153057247968152]

            for RKiter=1:5
    
                setupbc(Meshfile,Cells,SurfX,SurfY,SurfZ)
                conv_Flow   =   aRK[RKiter]*conv_Flow   +   dt*flow_eval(conv_scheme,visc_scheme,Cells,SurfX,SurfY,SurfZ,gamma,Pr,mu)
                Cells       =   Cells   +   bRK[RKiter]*conv_Flow
            
            end

        end
        #-------------------------------------------------------------------------



        #-------------------------------------------------------------------------
        # Write binary solution file
        #-------------------------------------------------------------------------
        if (it%w_itr)==0
            println("Writing file for time-step ",it)
            write(string(wfile,lpad(it,6,0),".vtrn"),Cells)
            #write_vtk(X,Y,Z,Cells,it,wfile)
        end
        #-------------------------------------------------------------------------

        println("Iteration ",it," finished")
        println("=================================================================================================\n")

        it  =   it+1

    end

end
