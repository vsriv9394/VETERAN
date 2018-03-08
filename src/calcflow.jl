function flow_eval(scheme,Cells,SurfX,SurfY,SurfZ,gamma)

    Flow    =   zeros(size(Cells))

    stateL  =   [0.0,0.0,0.0,0.0,0.0,0.0]
    stateR  =   [0.0,0.0,0.0,0.0,0.0,0.0]

    surf    =   [0.0,0.0,0.0,0.0]

    for i=1:size(SurfX,1)
        for j=1:size(SurfX,2)
            for k=1:size(SurfX,3)
                ext_side        =   1
                if i==1
                    ext_side    =   0
                end
                Area    =   sqrt(SurfX[i,j,k,1]^2+SurfX[i,j,k,2]^2+SurfX[i,j,k,3]^2)
                for m=1:5
                    stateL[m]   =   Cells[i  ,j+1,k+1,m]
                    stateR[m]   =   Cells[i+1,j+1,k+1,m]
                end
                for m=1:4
                    surf[m]     =   SurfX[i,j,k,m]
                end
                Flow_local  =   Area*setbcflux(scheme,stateL,stateR,gamma,surf,ext_side)
                for m=1:5
                    if i!=1
                        Flow[i  ,j+1,k+1,m] =   Flow[i  ,j+1,k+1,m]-Flow_local[m]/Cells[i  ,j+1,k+1,6]
                    end
                    if i!=size(SurfX,1)
                        Flow[i+1,j+1,k+1,m] =   Flow[i+1,j+1,k+1,m]+Flow_local[m]/Cells[i+1,j+1,k+1,6]
                    end
                end
            end
        end
    end

    for i=1:size(SurfY,1)
        for j=1:size(SurfY,2)
            for k=1:size(SurfY,3)
                ext_side        =   1
                if j==1
                    ext_side    =   0
                end
                Area    =   sqrt(SurfY[i,j,k,1]^2+SurfY[i,j,k,2]^2+SurfY[i,j,k,3]^2)
                for m=1:5
                    stateL[m]   =   Cells[i+1,j  ,k+1,m]
                    stateR[m]   =   Cells[i+1,j+1,k+1,m]
                end
                for m=1:4
                    surf[m]     =   SurfY[i,j,k,m]
                end
                Flow_local  =   Area*setbcflux(scheme,stateL,stateR,gamma,surf,ext_side)
                for m=1:5
                    if j!=1
                        Flow[i+1,j  ,k+1,m] =   Flow[i+1,j  ,k+1,m]-Flow_local[m]/Cells[i+1,j  ,k+1,6]
                    end
                    if j!=size(SurfY,2)
                        Flow[i+1,j+1,k+1,m] =   Flow[i+1,j+1,k+1,m]+Flow_local[m]/Cells[i+1,j+1,k+1,6]
                    end
                end
            end
        end
    end

    for i=1:size(SurfZ,1)
        for j=1:size(SurfZ,2)
            for k=1:size(SurfZ,3)
                ext_side        =   1
                if k==1
                    ext_side    =   0
                end
                Area    =   sqrt(SurfZ[i,j,k,1]^2+SurfZ[i,j,k,2]^2+SurfZ[i,j,k,3]^2)
                for m=1:5
                    stateL[m]   =   Cells[i+1,j+1,k  ,m]
                    stateR[m]   =   Cells[i+1,j+1,k+1,m]
                end
                for m=1:4
                    surf[m]     =   SurfZ[i,j,k,m]
                end
                Flow_local  =   Area*setbcflux(scheme,stateL,stateR,gamma,surf,ext_side)
                for m=1:5
                    if k!=1
                        Flow[i+1,j+1,k  ,m] =   Flow[i+1,j+1,k  ,m]-Flow_local[m]/Cells[i+1,j+1,k  ,6]
                    end
                    if k!=size(SurfZ,3)
                        Flow[i+1,j+1,k+1,m] =   Flow[i+1,j+1,k+1,m]+Flow_local[m]/Cells[i+1,j+1,k+1,6]
                    end
                end
            end
        end
    end

    return Flow
end
