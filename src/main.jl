function simulate(Meshfile,time_scheme,conv_scheme,gamma,R,dt,n_itr)
    
    Cells   =   []
    SurfX   =   []
    SurfY   =   []
    SurfZ   =   []
    X       =   []
    Y       =   []
    Z       =   []
   
    create_Geometry(Meshfile,Cells,SurfX,SurfY,SurfZ)

    for it = 1:n_itr
        if time_scheme=="Euler_explicit"
            conv_Flow   =   flow_eval(conv_scheme,Cells,SurfX,SurfY,SurfZ,gamma,R)
            Cells       =   Cells   +   dt*conv_Flow
        end
    end

end
