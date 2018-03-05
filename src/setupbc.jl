function setupbc(Meshfile,Cells,SurfX,SurfY,SurfZ)

    # Compile Periodic Boundary Conditions
    BC  =   readdlm(string(Meshfile,"bc"))      # Surface Element Indices,BC_Marker
    IP  =   readdlm(string(Meshfile,"inp"))     # line_number: BC_Marker,states,BC_type
    
    DimX    =   size(Cells,1)-1
    DimY    =   size(Cells,2)-1
    DimZ    =   size(Cells,3)-1
    
    DimSX   =   (DimY-1)*(DimZ-1)
    DimSY   =   (DimZ-1)*(DimX-1)
    DimSZ   =   (DimX-1)*(DimY-1)
    
    index   =   0
    
    for i=1:DimSX
    
        index   =   index+1
        i       =   BC[index,1]
        j       =   BC[index,2]
        k       =   BC[index,3]
        m       =   BC[index,4]
        
        if m<0
            Cells[1,j+1,k+1,1]  =   IP[-m,1]
            Cells[1,j+1,k+1,2]  =   IP[-m,2]
            Cells[1,j+1,k+1,3]  =   IP[-m,3]
            Cells[1,j+1,k+1,4]  =   IP[-m,4]
            Cells[1,j+1,k+1,5]  =   IP[-m,5]
            SurfX[1,j  ,k  ,4]  =   IP[-m,6]
        else
            Cells[1,j+1,k+1,1]  =   Cells[BC[m,1],BC[m,2],BC[m,3],1]
            Cells[1,j+1,k+1,2]  =   Cells[BC[m,1],BC[m,2],BC[m,3],2]
            Cells[1,j+1,k+1,3]  =   Cells[BC[m,1],BC[m,2],BC[m,3],3]
            Cells[1,j+1,k+1,4]  =   Cells[BC[m,1],BC[m,2],BC[m,3],4]
            Cells[1,j+1,k+1,5]  =   Cells[BC[m,1],BC[m,2],BC[m,3],5]
            SurfX[1,j  ,k  ,4]  =   0
        end
    
    end
    
    for i=1:DimSX
    
        index   =   index+1
        i       =   BC[index,1]
        j       =   BC[index,2]
        k       =   BC[index,3]
        m       =   BC[index,4]
        
        if m<0
            Cells[end,j+1,k+1,1]  =   IP[-m,1]
            Cells[end,j+1,k+1,2]  =   IP[-m,2]
            Cells[end,j+1,k+1,3]  =   IP[-m,3]
            Cells[end,j+1,k+1,4]  =   IP[-m,4]
            Cells[end,j+1,k+1,5]  =   IP[-m,5]
            SurfX[end,j  ,k  ,4]  =   IP[-m,6]
        else
            Cells[end,j+1,k+1,1]  =   Cells[BC[m,1],BC[m,2],BC[m,3],1]
            Cells[end,j+1,k+1,2]  =   Cells[BC[m,1],BC[m,2],BC[m,3],2]
            Cells[end,j+1,k+1,3]  =   Cells[BC[m,1],BC[m,2],BC[m,3],3]
            Cells[end,j+1,k+1,4]  =   Cells[BC[m,1],BC[m,2],BC[m,3],4]
            Cells[end,j+1,k+1,5]  =   Cells[BC[m,1],BC[m,2],BC[m,3],5]
            SurfX[end,j  ,k  ,4]  =   0
        end
    
    end
    
    for i=1:DimSY
    
        index   =   index+1
        i       =   BC[index,1]
        j       =   BC[index,2]
        k       =   BC[index,3]
        m       =   BC[index,4]
        
        if m<0
            Cells[i+1,1,k+1,1]  =   IP[-m,1]
            Cells[i+1,1,k+1,2]  =   IP[-m,2]
            Cells[i+1,1,k+1,3]  =   IP[-m,3]
            Cells[i+1,1,k+1,4]  =   IP[-m,4]
            Cells[i+1,1,k+1,5]  =   IP[-m,5]
            SurfY[i  ,1,k  ,4]  =   IP[-m,6]
        else
            Cells[i+1,1,k+1,1]  =   Cells[BC[m,1],BC[m,2],BC[m,3],1]
            Cells[i+1,1,k+1,2]  =   Cells[BC[m,1],BC[m,2],BC[m,3],2]
            Cells[i+1,1,k+1,3]  =   Cells[BC[m,1],BC[m,2],BC[m,3],3]
            Cells[i+1,1,k+1,4]  =   Cells[BC[m,1],BC[m,2],BC[m,3],4]
            Cells[i+1,1,k+1,5]  =   Cells[BC[m,1],BC[m,2],BC[m,3],5]
            SurfY[i  ,1,k  ,4]  =   0
        end
    
    end
    
    for i=1:DimSY
    
        index   =   index+1
        i       =   BC[index,1]
        j       =   BC[index,2]
        k       =   BC[index,3]
        m       =   BC[index,4]
        
        if m<0
            Cells[i+1,end,k+1,1]  =   IP[-m,1]
            Cells[i+1,end,k+1,2]  =   IP[-m,2]
            Cells[i+1,end,k+1,3]  =   IP[-m,3]
            Cells[i+1,end,k+1,4]  =   IP[-m,4]
            Cells[i+1,end,k+1,5]  =   IP[-m,5]
            SurfY[i  ,end,k  ,4]  =   IP[-m,6]
        else
            Cells[i+1,end,k+1,1]  =   Cells[BC[m,1],BC[m,2],BC[m,3],1]
            Cells[i+1,end,k+1,2]  =   Cells[BC[m,1],BC[m,2],BC[m,3],2]
            Cells[i+1,end,k+1,3]  =   Cells[BC[m,1],BC[m,2],BC[m,3],3]
            Cells[i+1,end,k+1,4]  =   Cells[BC[m,1],BC[m,2],BC[m,3],4]
            Cells[i+1,end,k+1,5]  =   Cells[BC[m,1],BC[m,2],BC[m,3],5]
            SurfY[i  ,end,k  ,4]  =   0
        end
    
    end
    
    for i=1:DimSZ
    
        index   =   index+1
        i       =   BC[index,1]
        j       =   BC[index,2]
        k       =   BC[index,3]
        m       =   BC[index,4]
        
        if m<0
            Cells[i+1,j+1,1,1]  =   IP[-m,1]
            Cells[i+1,j+1,1,2]  =   IP[-m,2]
            Cells[i+1,j+1,1,3]  =   IP[-m,3]
            Cells[i+1,j+1,1,4]  =   IP[-m,4]
            Cells[i+1,j+1,1,5]  =   IP[-m,5]
            SurfZ[i  ,j  ,1,4]  =   IP[-m,6]
        else
            Cells[i+1,j+1,1,1]  =   Cells[BC[m,1],BC[m,2],BC[m,3],1]
            Cells[i+1,j+1,1,2]  =   Cells[BC[m,1],BC[m,2],BC[m,3],2]
            Cells[i+1,j+1,1,3]  =   Cells[BC[m,1],BC[m,2],BC[m,3],3]
            Cells[i+1,j+1,1,4]  =   Cells[BC[m,1],BC[m,2],BC[m,3],4]
            Cells[i+1,j+1,1,5]  =   Cells[BC[m,1],BC[m,2],BC[m,3],5]
            SurfZ[i  ,j  ,1,4]  =   0
        end
    
    end
    
    for i=1:DimSZ
    
        index   =   index+1
        i       =   BC[index,1]
        j       =   BC[index,2]
        k       =   BC[index,3]
        m       =   BC[index,4]
        
        if m<0
            Cells[i+1,j+1,end,1]  =   IP[-m,1]
            Cells[i+1,j+1,end,2]  =   IP[-m,2]
            Cells[i+1,j+1,end,3]  =   IP[-m,3]
            Cells[i+1,j+1,end,4]  =   IP[-m,4]
            Cells[i+1,j+1,end,5]  =   IP[-m,5]
            SurfZ[i  ,j  ,end,4]  =   IP[-m,6]
        else
            Cells[i+1,j+1,end,1]  =   Cells[BC[m,1],BC[m,2],BC[m,3],1]
            Cells[i+1,j+1,end,2]  =   Cells[BC[m,1],BC[m,2],BC[m,3],2]
            Cells[i+1,j+1,end,3]  =   Cells[BC[m,1],BC[m,2],BC[m,3],3]
            Cells[i+1,j+1,end,4]  =   Cells[BC[m,1],BC[m,2],BC[m,3],4]
            Cells[i+1,j+1,end,5]  =   Cells[BC[m,1],BC[m,2],BC[m,3],5]
            SurfZ[i  ,j  ,end,4]  =   0
        end
    
    end

end
