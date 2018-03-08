using VETERAN

Mesh_file_name      =   "new"
Restart_file_name   =   "initcond.vtrn"
Solution_file_name  =   "Outputs/new"

Conv_scheme         =   "Roe"
Visc_scheme         =   ""
#Time_scheme         =   "Euler_explicit"
Time_scheme         =   "RK4"

dt                  =   1e-4
stop_time           =   1
write_freq          =   1

Gas_prop_gamma      =   1.4
Gas_prop_R          =   287

#################################################################################################

Files               =   [Mesh_file_name,Restart_file_name,Solution_file_name]
Schemes             =   [Conv_scheme,Visc_scheme,Time_scheme]
Time_specs          =   [dt,stop_time,write_freq]
Gas_prop            =   [Gas_prop_gamma,Gas_prop_R]

VETERAN.simulate(Files,Schemes,Time_specs,Gas_prop)
