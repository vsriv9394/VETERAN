function create_Geometry(Meshfile,Restartfile)
    
    Mesh    = readdlm(string(Meshfile,".mesh"))
    
    DimX    = convert(Int, Mesh[1,1])
    DimY    = convert(Int, Mesh[1,2])
    DimZ    = convert(Int, Mesh[1,3])
    
    X       = reshape(Mesh[2:end,1],(DimX,DimY,DimZ))
    Y       = reshape(Mesh[2:end,2],(DimX,DimY,DimZ))
    Z       = reshape(Mesh[2:end,3],(DimX,DimY,DimZ))
    
    println("Beginning Mesh Setup and Geometry Initialization")
    Cells   = reshape(reinterpret(Float64,read(Restartfile)),DimX+1,DimY+1,DimZ+1,6)
    println("Initial/Restart file loaded")

    SurfX   = zeros(DimX  , DimY-1, DimZ-1, 4) # Normal Vector & Area & Boundary Type
    SurfY   = zeros(DimX-1, DimY  , DimZ-1, 4) # Normal Vector & Area & Boundary Type
    SurfZ   = zeros(DimX-1, DimY-1, DimZ  , 4) # Normal Vector & Area & Boundary Type
    
    for i=2:DimX
        for j=2:DimY
            for k=2:DimZ
                x000 = X[i-1,j-1,k-1]
                y000 = Y[i-1,j-1,k-1]
                z000 = Z[i-1,j-1,k-1]
                
                x100 = X[i  ,j-1,k-1]
                y100 = Y[i  ,j-1,k-1]
                z100 = Z[i  ,j-1,k-1]
                
                x010 = X[i-1,j  ,k-1]
                y010 = Y[i-1,j  ,k-1]
                z010 = Z[i-1,j  ,k-1]
                
                x001 = X[i-1,j-1,k  ]
                y001 = Y[i-1,j-1,k  ]
                z001 = Z[i-1,j-1,k  ]
                
                x110 = X[i  ,j  ,k-1]
                y110 = Y[i  ,j  ,k-1]
                z110 = Z[i  ,j  ,k-1]
                
                x101 = X[i  ,j-1,k  ]
                y101 = Y[i  ,j-1,k  ]
                z101 = Z[i  ,j-1,k  ]
                
                x011 = X[i  ,j  ,k-1]
                y011 = Y[i  ,j  ,k-1]
                z011 = Z[i  ,j  ,k-1]
                
                x111 = X[i  ,j  ,k  ]
                y111 = Y[i  ,j  ,k  ]
                z111 = Z[i  ,j  ,k  ]
    
                Vol1 = abs(1/6*det([x100-x000 x010-x000 x001-x000; y100-y000 y010-y000 y001-y000; z100-z000 z010-z000 z001-z000]))
                Vol2 = abs(1/6*det([x100-x110 x010-x110 x111-x110; y100-y110 y010-y110 y111-y110; z100-z110 z010-z110 z111-z110]))
                Vol3 = abs(1/6*det([x001-x101 x111-x101 x100-x101; y001-y101 y111-y101 y100-y101; z001-z101 z111-z101 z100-z101]))
                Vol4 = abs(1/6*det([x001-x011 x111-x011 x010-x011; y001-y011 y111-y011 y010-y011; z001-z011 z111-z011 z010-z011]))
                Vol5 = abs(1/6*det([x100-x111 x010-x111 x001-x111; y100-y111 y010-y111 y001-y111; z100-z111 z010-z111 z001-z111]))
    
                Cells[i,j,k,6] = Vol1+Vol2+Vol3+Vol4+Vol5
            end
        end
    end
    
    for i=1:DimX
        for j=1:DimY-1
            for k=1:DimZ-1
                x00 = X[i  ,j  ,k  ]
                y00 = Y[i  ,j  ,k  ]
                z00 = Z[i  ,j  ,k  ]
    
                x01 = X[i  ,j  ,k+1]
                y01 = Y[i  ,j  ,k+1]
                z01 = Z[i  ,j  ,k+1]
    
                x10 = X[i  ,j+1,k  ]
                y10 = Y[i  ,j+1,k  ]
                z10 = Z[i  ,j+1,k  ]
    
                x11 = X[i  ,j+1,k+1]
                y11 = Y[i  ,j+1,k+1]
                z11 = Z[i  ,j+1,k+1]
    
                Area1   = 0.5*(det([y10-y00 z10-z00; y11-y10 z11-z10])+det([y01-y11 z01-z11; y00-y01 z00-z01]))
                Area2   = 0.5*(det([z10-z00 x10-x00; z11-z10 x11-x10])+det([z01-z11 x01-x11; z00-z01 x00-x01]))
                Area3   = 0.5*(det([x10-x00 y10-y00; x11-x10 y11-y10])+det([x01-x11 y01-y11; x00-x01 y00-y01]))
                SurfX[i,j,k,1]  = Area1
                SurfX[i,j,k,2]  = Area2
                SurfX[i,j,k,3]  = Area3
            end
        end
    end
    
    for i=1:DimX-1
        for j=1:DimY
            for k=1:DimZ-1
                x00 = X[i  ,j  ,k  ]
                y00 = Y[i  ,j  ,k  ]
                z00 = Z[i  ,j  ,k  ]
    
                x01 = X[i  ,j  ,k+1]
                y01 = Y[i  ,j  ,k+1]
                z01 = Z[i  ,j  ,k+1]
    
                x10 = X[i+1,j  ,k  ]
                y10 = Y[i+1,j  ,k  ]
                z10 = Z[i+1,j  ,k  ]
    
                x11 = X[i+1,j  ,k+1]
                y11 = Y[i+1,j  ,k+1]
                z11 = Z[i+1,j  ,k+1]
    
                Area1   = -0.5*(det([y10-y00 z10-z00; y11-y10 z11-z10])+det([y01-y11 z01-z11; y00-y01 z00-z01]))
                Area2   = -0.5*(det([z10-z00 x10-x00; z11-z10 x11-x10])+det([z01-z11 x01-x11; z00-z01 x00-x01]))
                Area3   = -0.5*(det([x10-x00 y10-y00; x11-x10 y11-y10])+det([x01-x11 y01-y11; x00-x01 y00-y01]))
                SurfY[i,j,k,1]  = Area1
                SurfY[i,j,k,2]  = Area2
                SurfY[i,j,k,3]  = Area3
            end
        end
    end
    
    for i=1:DimX-1
        for j=1:DimY-1
            for k=1:DimZ
                x00 = X[i  ,j  ,k  ]
                y00 = Y[i  ,j  ,k  ]
                z00 = Z[i  ,j  ,k  ]
    
                x01 = X[i  ,j+1,k  ]
                y01 = Y[i  ,j+1,k  ]
                z01 = Z[i  ,j+1,k  ]
    
                x10 = X[i+1,j  ,k  ]
                y10 = Y[i+1,j  ,k  ]
                z10 = Z[i+1,j  ,k  ]
    
                x11 = X[i+1,j+1,k  ]
                y11 = Y[i+1,j+1,k  ]
                z11 = Z[i+1,j+1,k  ]
    
                Area1   = 0.5*(det([y10-y00 z10-z00; y11-y10 z11-z10])+det([y01-y11 z01-z11; y00-y01 z00-z01]))
                Area2   = 0.5*(det([z10-z00 x10-x00; z11-z10 x11-x10])+det([z01-z11 x01-x11; z00-z01 x00-x01]))
                Area3   = 0.5*(det([x10-x00 y10-y00; x11-x10 y11-y10])+det([x01-x11 y01-y11; x00-x01 y00-y01]))
                SurfZ[i,j,k,1]  = Area1
                SurfZ[i,j,k,2]  = Area2
                SurfZ[i,j,k,3]  = Area3
            end
        end
    end

    println("Mesh Setup and Geometry Initialization Finished")
    println("Maximum X-Area\t\t\t",maximum(SurfX[:,:,:,1]),"\nMinimum X-Area\t\t\t",minimum(SurfX[:,:,:,1]),"\n")
    println("Maximum Y-Area\t\t\t",maximum(SurfY[:,:,:,2]),"\nMinimum Y-Area\t\t\t",minimum(SurfY[:,:,:,2]),"\n")
    println("Maximum Z-Area\t\t\t",maximum(SurfZ[:,:,:,3]),"\nMinimum Z-Area\t\t\t",minimum(SurfZ[:,:,:,3]),"\n")
    println("Maximum Cell Volume\t\t",maximum(Cells[2:end-1,2:end-1,2:end-1,6]),"\nMinimum Cell Volume\t\t",minimum(Cells[2:end-1,2:end-1,2:end-1,6]),"\n")

    return Cells,SurfX,SurfY,SurfZ,X,Y,Z
    
end
