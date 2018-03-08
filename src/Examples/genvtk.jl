using VETERAN

Mesh_file_name      =   "new"
Restart_file_name   =   "initcond.vtrn"
Solution_file_name  =   "Outputs/new"

#################################################################################################

Files               =   [Mesh_file_name,Restart_file_name,Solution_file_name]

VETERAN.write_vtk(Files)
