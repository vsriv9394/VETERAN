function simulate(Files,Schemes,Time_specs,Gas_prop)
    
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
    R           =   Gas_prop[2]

    Cells,SurfX,SurfY,SurfZ,X,Y,Z = create_Geometry(Meshfile,Restartfile)

    #write_vtk(X,Y,Z,Cells,0,wfile)

    setupbc(Meshfile,Cells,SurfX,SurfY,SurfZ)

    it = 1
    
    while it*dt<=stop

        println("=================================================================================================")
        println("Iteration ",it," initialized")

        if time_scheme=="Euler_explicit"
            setupbc(Meshfile,Cells,SurfX,SurfY,SurfZ)
            conv_Flow   =   flow_eval(conv_scheme,Cells,SurfX,SurfY,SurfZ,gamma)
            Cells       =   Cells   +   dt*conv_Flow
        elseif time_scheme=="RK2"
            setupbc(Meshfile,Cells,SurfX,SurfY,SurfZ)
            conv_Flow1  =   flow_eval(conv_scheme,Cells,SurfX,SurfY,SurfZ,gamma)
            Cells1      =   Cells   +   dt*conv_Flow1
            setupbc(Meshfile,Cells1,SurfX,SurfY,SurfZ)
            conv_Flow2  =   flow_eval(conv_scheme,Cells1,SurfX,SurfY,SurfZ,gamma)
            Cells       =   Cells   +   dt*(conv_Flow1+conv_Flow2)*0.5
        elseif time_scheme=="RK4"
            setupbc(Meshfile,Cells,SurfX,SurfY,SurfZ)
            conv_Flow1  =   flow_eval(conv_scheme,Cells,SurfX,SurfY,SurfZ,gamma)
            Cells1      =   Cells   +   0.5*dt*conv_Flow1
            setupbc(Meshfile,Cells1,SurfX,SurfY,SurfZ)
            conv_Flow2  =   flow_eval(conv_scheme,Cells1,SurfX,SurfY,SurfZ,gamma)
            Cells1      =   Cells   +   0.5*dt*conv_Flow2
            setupbc(Meshfile,Cells1,SurfX,SurfY,SurfZ)
            conv_Flow3  =   flow_eval(conv_scheme,Cells1,SurfX,SurfY,SurfZ,gamma)
            Cells1      =   Cells   +   dt*conv_Flow3
            setupbc(Meshfile,Cells1,SurfX,SurfY,SurfZ)
            conv_Flow4  =   flow_eval(conv_scheme,Cells1,SurfX,SurfY,SurfZ,gamma)
            Cells       =   Cells   +   dt/6*(conv_Flow1+2*conv_Flow2+2*conv_Flow3+conv_Flow2)
        end
        if (it%w_itr)==0
            println("Writing file for time-step ",it)
            write(string(wfile,lpad(it,6,0),".vtrn"),Cells)
            #write_vtk(X,Y,Z,Cells,it,wfile)
        end

        println("Iteration ",it," finished")
        println("=================================================================================================\n")

        it  =   it+1

    end

end
