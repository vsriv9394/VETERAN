using WriteVTK

function write_vtk(X,Y,Z,Cells,index,wfile)

    vtkfile =   vtk_grid(string(wfile,index), X, Y, Z)
    vtk_cell_data(vtkfile,Cells[:,:,:,1],"rho")
    vtk_cell_data(vtkfile,Cells[:,:,:,2]./Cells[:,:,:,1],"U")
    vtk_cell_data(vtkfile,Cells[:,:,:,3]./Cells[:,:,:,1],"V")
    vtk_cell_data(vtkfile,Cells[:,:,:,4]./Cells[:,:,:,1],"W")
    vtk_cell_data(vtkfile,Cells[:,:,:,5],"rhoE")

end
