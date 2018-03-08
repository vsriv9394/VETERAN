using WriteVTK

function write_vtk(Files)

    Meshfile    =   Files[1]
    Restartfile =   Files[2]
    wfile       =   Files[3]

    nothing,nothing,nothing,nothing,X,Y,Z = create_Geometry(Meshfile,Restartfile)

    for index=parse(Int,ARGS[1]):parse(Int,ARGS[2]):parse(Int,ARGS[3])
        vtkfile =   vtk_grid(string(wfile,lpad(index,6,0)), X, Y, Z)
        Cells   =   reshape(reinterpret(Float64,read(string(wfile,lpad(index,6,0),".vtrn"))),size(X,1)+1,size(X,2)+1,size(X,3)+1,6)
        vtk_cell_data(vtkfile,Cells[2:end-1,2:end-1,2:end-1,1],"rho")
        vtk_cell_data(vtkfile,Cells[2:end-1,2:end-1,2:end-1,2]./Cells[2:end-1,2:end-1,2:end-1,1],"U")
        vtk_cell_data(vtkfile,Cells[2:end-1,2:end-1,2:end-1,3]./Cells[2:end-1,2:end-1,2:end-1,1],"V")
        vtk_cell_data(vtkfile,Cells[2:end-1,2:end-1,2:end-1,4]./Cells[2:end-1,2:end-1,2:end-1,1],"W")
        vtk_cell_data(vtkfile,Cells[2:end-1,2:end-1,2:end-1,5],"rhoE")
        vtk_save(vtkfile)
        println(index," Done")
    end

end
