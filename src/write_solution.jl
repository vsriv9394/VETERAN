using WriteVTK

function write_vtk(X,Y,Z,Cells,index,wfile)

    vtkfile =   vtk_grid(string(wfile,index), X, Y, Z)
    vtk_cell_data(vtkfile,Cells[2:end-1,2:end-1,2:end-1,1],"rho")
    vtk_cell_data(vtkfile,Cells[2:end-1,2:end-1,2:end-1,2]./Cells[2:end-1,2:end-1,2:end-1,1],"U")
    vtk_cell_data(vtkfile,Cells[2:end-1,2:end-1,2:end-1,3]./Cells[2:end-1,2:end-1,2:end-1,1],"V")
    vtk_cell_data(vtkfile,Cells[2:end-1,2:end-1,2:end-1,4]./Cells[2:end-1,2:end-1,2:end-1,1],"W")
    vtk_cell_data(vtkfile,Cells[2:end-1,2:end-1,2:end-1,5],"rhoE")
    vtk_save(vtkfile)

end
