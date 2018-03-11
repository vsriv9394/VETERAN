function reconstruct(Cells,SurfX,SurfY,SurfZ,gamma)

    siz     =   collect(size(Cells))
    Grads   =   zeros(siz[1],siz[2],siz[3],4,3)
    stats   =   zeros(siz[1],siz[2],siz[3],4)

    for i=1:siz[1]
        for j=1:siz[2]
            for k=1:siz[3]
                stats[i,j,k,1]     =    Cells[i,j,k,2]/Cells[i,j,k,1]
                stats[i,j,k,2]     =    Cells[i,j,k,3]/Cells[i,j,k,1]
                stats[i,j,k,3]     =    Cells[i,j,k,4]/Cells[i,j,k,1]
                stats[i,j,k,4]     =    gamma*(Cells[i,j,k,5] - 0.5*(Cells[i,j,k,2]^2+Cells[i,j,k,3]^2+Cells[i,j,k,4]^2)/Cells[i,j,k,1])/Cells[i,j,k,1]
            end
        end
    end

    for i=2:siz[1]-1
        for j=2:siz[2]-1
            for k=2:siz[3]-1
                
                Grads[i,j,k,1,1]   =    Grads[i,j,k,1,1] - SurfX[i-1,j-1,k-1,1] * ( stats[i,j,k,1] + stats[i-1,j  ,k  ,1] )*0.5
                Grads[i,j,k,1,1]   =    Grads[i,j,k,1,1] - SurfY[i-1,j-1,k-1,1] * ( stats[i,j,k,1] + stats[i  ,j-1,k  ,1] )*0.5
                Grads[i,j,k,1,1]   =    Grads[i,j,k,1,1] - SurfZ[i-1,j-1,k-1,1] * ( stats[i,j,k,1] + stats[i  ,j  ,k-1,1] )*0.5
                Grads[i,j,k,1,1]   =    Grads[i,j,k,1,1] + SurfX[i  ,j-1,k-1,1] * ( stats[i,j,k,1] + stats[i+1,j  ,k  ,1] )*0.5
                Grads[i,j,k,1,1]   =    Grads[i,j,k,1,1] + SurfY[i-1,j  ,k-1,1] * ( stats[i,j,k,1] + stats[i  ,j+1,k  ,1] )*0.5
                Grads[i,j,k,1,1]   =    Grads[i,j,k,1,1] + SurfZ[i-1,j-1,k  ,1] * ( stats[i,j,k,1] + stats[i  ,j  ,k+1,1] )*0.5
                Grads[i,j,k,1,2]   =    Grads[i,j,k,1,2] - SurfX[i-1,j-1,k-1,2] * ( stats[i,j,k,1] + stats[i-1,j  ,k  ,1] )*0.5
                Grads[i,j,k,1,2]   =    Grads[i,j,k,1,2] - SurfY[i-1,j-1,k-1,2] * ( stats[i,j,k,1] + stats[i  ,j-1,k  ,1] )*0.5
                Grads[i,j,k,1,2]   =    Grads[i,j,k,1,2] - SurfZ[i-1,j-1,k-1,2] * ( stats[i,j,k,1] + stats[i  ,j  ,k-1,1] )*0.5
                Grads[i,j,k,1,2]   =    Grads[i,j,k,1,2] + SurfX[i  ,j-1,k-1,2] * ( stats[i,j,k,1] + stats[i+1,j  ,k  ,1] )*0.5
                Grads[i,j,k,1,2]   =    Grads[i,j,k,1,2] + SurfY[i-1,j  ,k-1,2] * ( stats[i,j,k,1] + stats[i  ,j+1,k  ,1] )*0.5
                Grads[i,j,k,1,2]   =    Grads[i,j,k,1,2] + SurfZ[i-1,j-1,k  ,2] * ( stats[i,j,k,1] + stats[i  ,j  ,k+1,1] )*0.5
                Grads[i,j,k,1,3]   =    Grads[i,j,k,1,3] - SurfX[i-1,j-1,k-1,3] * ( stats[i,j,k,1] + stats[i-1,j  ,k  ,1] )*0.5
                Grads[i,j,k,1,3]   =    Grads[i,j,k,1,3] - SurfY[i-1,j-1,k-1,3] * ( stats[i,j,k,1] + stats[i  ,j-1,k  ,1] )*0.5
                Grads[i,j,k,1,3]   =    Grads[i,j,k,1,3] - SurfZ[i-1,j-1,k-1,3] * ( stats[i,j,k,1] + stats[i  ,j  ,k-1,1] )*0.5
                Grads[i,j,k,1,3]   =    Grads[i,j,k,1,3] + SurfX[i  ,j-1,k-1,3] * ( stats[i,j,k,1] + stats[i+1,j  ,k  ,1] )*0.5
                Grads[i,j,k,1,3]   =    Grads[i,j,k,1,3] + SurfY[i-1,j  ,k-1,3] * ( stats[i,j,k,1] + stats[i  ,j+1,k  ,1] )*0.5
                Grads[i,j,k,1,3]   =    Grads[i,j,k,1,3] + SurfZ[i-1,j-1,k  ,3] * ( stats[i,j,k,1] + stats[i  ,j  ,k+1,1] )*0.5
                
                
                Grads[i,j,k,2,1]   =    Grads[i,j,k,2,1] - SurfX[i-1,j-1,k-1,1] * ( stats[i,j,k,2] + stats[i-1,j  ,k  ,2] )*0.5
                Grads[i,j,k,2,1]   =    Grads[i,j,k,2,1] - SurfY[i-1,j-1,k-1,1] * ( stats[i,j,k,2] + stats[i  ,j-1,k  ,2] )*0.5
                Grads[i,j,k,2,1]   =    Grads[i,j,k,2,1] - SurfZ[i-1,j-1,k-1,1] * ( stats[i,j,k,2] + stats[i  ,j  ,k-1,2] )*0.5
                Grads[i,j,k,2,1]   =    Grads[i,j,k,2,1] + SurfX[i  ,j-1,k-1,1] * ( stats[i,j,k,2] + stats[i+1,j  ,k  ,2] )*0.5
                Grads[i,j,k,2,1]   =    Grads[i,j,k,2,1] + SurfY[i-1,j  ,k-1,1] * ( stats[i,j,k,2] + stats[i  ,j+1,k  ,2] )*0.5
                Grads[i,j,k,2,1]   =    Grads[i,j,k,2,1] + SurfZ[i-1,j-1,k  ,1] * ( stats[i,j,k,2] + stats[i  ,j  ,k+1,2] )*0.5
                Grads[i,j,k,2,2]   =    Grads[i,j,k,2,2] - SurfX[i-1,j-1,k-1,2] * ( stats[i,j,k,2] + stats[i-1,j  ,k  ,2] )*0.5
                Grads[i,j,k,2,2]   =    Grads[i,j,k,2,2] - SurfY[i-1,j-1,k-1,2] * ( stats[i,j,k,2] + stats[i  ,j-1,k  ,2] )*0.5
                Grads[i,j,k,2,2]   =    Grads[i,j,k,2,2] - SurfZ[i-1,j-1,k-1,2] * ( stats[i,j,k,2] + stats[i  ,j  ,k-1,2] )*0.5
                Grads[i,j,k,2,2]   =    Grads[i,j,k,2,2] + SurfX[i  ,j-1,k-1,2] * ( stats[i,j,k,2] + stats[i+1,j  ,k  ,2] )*0.5
                Grads[i,j,k,2,2]   =    Grads[i,j,k,2,2] + SurfY[i-1,j  ,k-1,2] * ( stats[i,j,k,2] + stats[i  ,j+1,k  ,2] )*0.5
                Grads[i,j,k,2,2]   =    Grads[i,j,k,2,2] + SurfZ[i-1,j-1,k  ,2] * ( stats[i,j,k,2] + stats[i  ,j  ,k+1,2] )*0.5
                Grads[i,j,k,2,3]   =    Grads[i,j,k,2,3] - SurfX[i-1,j-1,k-1,3] * ( stats[i,j,k,2] + stats[i-1,j  ,k  ,2] )*0.5
                Grads[i,j,k,2,3]   =    Grads[i,j,k,2,3] - SurfY[i-1,j-1,k-1,3] * ( stats[i,j,k,2] + stats[i  ,j-1,k  ,2] )*0.5
                Grads[i,j,k,2,3]   =    Grads[i,j,k,2,3] - SurfZ[i-1,j-1,k-1,3] * ( stats[i,j,k,2] + stats[i  ,j  ,k-1,2] )*0.5
                Grads[i,j,k,2,3]   =    Grads[i,j,k,2,3] + SurfX[i  ,j-1,k-1,3] * ( stats[i,j,k,2] + stats[i+1,j  ,k  ,2] )*0.5
                Grads[i,j,k,2,3]   =    Grads[i,j,k,2,3] + SurfY[i-1,j  ,k-1,3] * ( stats[i,j,k,2] + stats[i  ,j+1,k  ,2] )*0.5
                Grads[i,j,k,2,3]   =    Grads[i,j,k,2,3] + SurfZ[i-1,j-1,k  ,3] * ( stats[i,j,k,2] + stats[i  ,j  ,k+1,2] )*0.5
                
                
                Grads[i,j,k,3,1]   =    Grads[i,j,k,3,1] - SurfX[i-1,j-1,k-1,1] * ( stats[i,j,k,3] + stats[i-1,j  ,k  ,3] )*0.5
                Grads[i,j,k,3,1]   =    Grads[i,j,k,3,1] - SurfY[i-1,j-1,k-1,1] * ( stats[i,j,k,3] + stats[i  ,j-1,k  ,3] )*0.5
                Grads[i,j,k,3,1]   =    Grads[i,j,k,3,1] - SurfZ[i-1,j-1,k-1,1] * ( stats[i,j,k,3] + stats[i  ,j  ,k-1,3] )*0.5
                Grads[i,j,k,3,1]   =    Grads[i,j,k,3,1] + SurfX[i  ,j-1,k-1,1] * ( stats[i,j,k,3] + stats[i+1,j  ,k  ,3] )*0.5
                Grads[i,j,k,3,1]   =    Grads[i,j,k,3,1] + SurfY[i-1,j  ,k-1,1] * ( stats[i,j,k,3] + stats[i  ,j+1,k  ,3] )*0.5
                Grads[i,j,k,3,1]   =    Grads[i,j,k,3,1] + SurfZ[i-1,j-1,k  ,1] * ( stats[i,j,k,3] + stats[i  ,j  ,k+1,3] )*0.5
                Grads[i,j,k,3,2]   =    Grads[i,j,k,3,2] - SurfX[i-1,j-1,k-1,2] * ( stats[i,j,k,3] + stats[i-1,j  ,k  ,3] )*0.5
                Grads[i,j,k,3,2]   =    Grads[i,j,k,3,2] - SurfY[i-1,j-1,k-1,2] * ( stats[i,j,k,3] + stats[i  ,j-1,k  ,3] )*0.5
                Grads[i,j,k,3,2]   =    Grads[i,j,k,3,2] - SurfZ[i-1,j-1,k-1,2] * ( stats[i,j,k,3] + stats[i  ,j  ,k-1,3] )*0.5
                Grads[i,j,k,3,2]   =    Grads[i,j,k,3,2] + SurfX[i  ,j-1,k-1,2] * ( stats[i,j,k,3] + stats[i+1,j  ,k  ,3] )*0.5
                Grads[i,j,k,3,2]   =    Grads[i,j,k,3,2] + SurfY[i-1,j  ,k-1,2] * ( stats[i,j,k,3] + stats[i  ,j+1,k  ,3] )*0.5
                Grads[i,j,k,3,2]   =    Grads[i,j,k,3,2] + SurfZ[i-1,j-1,k  ,2] * ( stats[i,j,k,3] + stats[i  ,j  ,k+1,3] )*0.5
                Grads[i,j,k,3,3]   =    Grads[i,j,k,3,3] - SurfX[i-1,j-1,k-1,3] * ( stats[i,j,k,3] + stats[i-1,j  ,k  ,3] )*0.5
                Grads[i,j,k,3,3]   =    Grads[i,j,k,3,3] - SurfY[i-1,j-1,k-1,3] * ( stats[i,j,k,3] + stats[i  ,j-1,k  ,3] )*0.5
                Grads[i,j,k,3,3]   =    Grads[i,j,k,3,3] - SurfZ[i-1,j-1,k-1,3] * ( stats[i,j,k,3] + stats[i  ,j  ,k-1,3] )*0.5
                Grads[i,j,k,3,3]   =    Grads[i,j,k,3,3] + SurfX[i  ,j-1,k-1,3] * ( stats[i,j,k,3] + stats[i+1,j  ,k  ,3] )*0.5
                Grads[i,j,k,3,3]   =    Grads[i,j,k,3,3] + SurfY[i-1,j  ,k-1,3] * ( stats[i,j,k,3] + stats[i  ,j+1,k  ,3] )*0.5
                Grads[i,j,k,3,3]   =    Grads[i,j,k,3,3] + SurfZ[i-1,j-1,k  ,3] * ( stats[i,j,k,3] + stats[i  ,j  ,k+1,3] )*0.5
                
                
                Grads[i,j,k,4,1]   =    Grads[i,j,k,4,1] - SurfX[i-1,j-1,k-1,1] * ( stats[i,j,k,4] + stats[i-1,j  ,k  ,4] )*0.5
                Grads[i,j,k,4,1]   =    Grads[i,j,k,4,1] - SurfY[i-1,j-1,k-1,1] * ( stats[i,j,k,4] + stats[i  ,j-1,k  ,4] )*0.5
                Grads[i,j,k,4,1]   =    Grads[i,j,k,4,1] - SurfZ[i-1,j-1,k-1,1] * ( stats[i,j,k,4] + stats[i  ,j  ,k-1,4] )*0.5
                Grads[i,j,k,4,1]   =    Grads[i,j,k,4,1] + SurfX[i  ,j-1,k-1,1] * ( stats[i,j,k,4] + stats[i+1,j  ,k  ,4] )*0.5
                Grads[i,j,k,4,1]   =    Grads[i,j,k,4,1] + SurfY[i-1,j  ,k-1,1] * ( stats[i,j,k,4] + stats[i  ,j+1,k  ,4] )*0.5
                Grads[i,j,k,4,1]   =    Grads[i,j,k,4,1] + SurfZ[i-1,j-1,k  ,1] * ( stats[i,j,k,4] + stats[i  ,j  ,k+1,4] )*0.5
                Grads[i,j,k,4,2]   =    Grads[i,j,k,4,2] - SurfX[i-1,j-1,k-1,2] * ( stats[i,j,k,4] + stats[i-1,j  ,k  ,4] )*0.5
                Grads[i,j,k,4,2]   =    Grads[i,j,k,4,2] - SurfY[i-1,j-1,k-1,2] * ( stats[i,j,k,4] + stats[i  ,j-1,k  ,4] )*0.5
                Grads[i,j,k,4,2]   =    Grads[i,j,k,4,2] - SurfZ[i-1,j-1,k-1,2] * ( stats[i,j,k,4] + stats[i  ,j  ,k-1,4] )*0.5
                Grads[i,j,k,4,2]   =    Grads[i,j,k,4,2] + SurfX[i  ,j-1,k-1,2] * ( stats[i,j,k,4] + stats[i+1,j  ,k  ,4] )*0.5
                Grads[i,j,k,4,2]   =    Grads[i,j,k,4,2] + SurfY[i-1,j  ,k-1,2] * ( stats[i,j,k,4] + stats[i  ,j+1,k  ,4] )*0.5
                Grads[i,j,k,4,2]   =    Grads[i,j,k,4,2] + SurfZ[i-1,j-1,k  ,2] * ( stats[i,j,k,4] + stats[i  ,j  ,k+1,4] )*0.5
                Grads[i,j,k,4,3]   =    Grads[i,j,k,4,3] - SurfX[i-1,j-1,k-1,3] * ( stats[i,j,k,4] + stats[i-1,j  ,k  ,4] )*0.5
                Grads[i,j,k,4,3]   =    Grads[i,j,k,4,3] - SurfY[i-1,j-1,k-1,3] * ( stats[i,j,k,4] + stats[i  ,j-1,k  ,4] )*0.5
                Grads[i,j,k,4,3]   =    Grads[i,j,k,4,3] - SurfZ[i-1,j-1,k-1,3] * ( stats[i,j,k,4] + stats[i  ,j  ,k-1,4] )*0.5
                Grads[i,j,k,4,3]   =    Grads[i,j,k,4,3] + SurfX[i  ,j-1,k-1,3] * ( stats[i,j,k,4] + stats[i+1,j  ,k  ,4] )*0.5
                Grads[i,j,k,4,3]   =    Grads[i,j,k,4,3] + SurfY[i-1,j  ,k-1,3] * ( stats[i,j,k,4] + stats[i  ,j+1,k  ,4] )*0.5
                Grads[i,j,k,4,3]   =    Grads[i,j,k,4,3] + SurfZ[i-1,j-1,k  ,3] * ( stats[i,j,k,4] + stats[i  ,j  ,k+1,4] )*0.5

                Grads[i,j,k,1,1]   =    Grads[i,j,k,1,1]/Cells[i,j,k,6]
                Grads[i,j,k,1,2]   =    Grads[i,j,k,1,1]/Cells[i,j,k,6]
                Grads[i,j,k,1,3]   =    Grads[i,j,k,1,1]/Cells[i,j,k,6]
                
                Grads[i,j,k,2,1]   =    Grads[i,j,k,2,1]/Cells[i,j,k,6]
                Grads[i,j,k,2,2]   =    Grads[i,j,k,2,2]/Cells[i,j,k,6]
                Grads[i,j,k,2,3]   =    Grads[i,j,k,2,3]/Cells[i,j,k,6]
                
                Grads[i,j,k,3,1]   =    Grads[i,j,k,3,1]/Cells[i,j,k,6]
                Grads[i,j,k,3,2]   =    Grads[i,j,k,3,2]/Cells[i,j,k,6]
                Grads[i,j,k,3,3]   =    Grads[i,j,k,3,3]/Cells[i,j,k,6]
                
                Grads[i,j,k,4,1]   =    Grads[i,j,k,4,1]/Cells[i,j,k,6]
                Grads[i,j,k,4,2]   =    Grads[i,j,k,4,2]/Cells[i,j,k,6]
                Grads[i,j,k,4,3]   =    Grads[i,j,k,4,3]/Cells[i,j,k,6]

            end
        end
    end

    return Grads

end
