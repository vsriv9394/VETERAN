function simulate(Meshfile,time_scheme,conv_scheme,gamma,R,dt,n_itr,w_itr,wfile)
    
    Cells,SurfX,SurfY,SurfZ,X,Y,Z = create_Geometry(Meshfile)

    println("Writing the initial file\n")
    write_vtk(X,Y,Z,Cells,0,wfile)

    setupbc(Meshfile,Cells,SurfX,SurfY,SurfZ)
    
    for it = 1:n_itr

        println("=================================================================================================")
        println("Iteration ",it," initialized")

        println("Boundary Conditions setup initialized")
        setupbc(Meshfile,Cells,SurfX,SurfY,SurfZ)
        println("Boundary Conditions setup finished")

        if time_scheme=="Euler_explicit"
            conv_Flow   =   flow_eval(conv_scheme,Cells,SurfX,SurfY,SurfZ,gamma,R)
            Cells       =   Cells   +   dt*conv_Flow
        end
        if (it%w_itr)==0
            println("Writing file for time-step ",it)
            write_vtk(X,Y,Z,Cells,it,wfile)
        end

        println("Iteration ",it," finished")
        println("=================================================================================================\n")

    end

end
