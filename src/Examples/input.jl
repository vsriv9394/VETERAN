using VETERAN

Mesh_file_name      =   "new"
Restart_file_name   =   "initcond.vtrn"
Solution_file_name  =   "Outputs/new"

Conv_scheme         =   "Roe"
Visc_scheme         =   "Dirty"
Time_scheme         =   "LSRK4"

dt                  =   1e-2
stop_time           =   30
write_freq          =   10

Gas_prop_gamma      =   1.4
Gas_prop_Pr         =   0.72
Gas_prop_mu         =   1.81e-5

#################################################################################################

Files               =   [Mesh_file_name,Restart_file_name,Solution_file_name]
Schemes             =   [Conv_scheme,Visc_scheme,Time_scheme]
Time_specs          =   [dt,stop_time,write_freq]
Gas_prop            =   [Gas_prop_gamma,Gas_prop_Pr,Gas_prop_mu]

VETERAN.simulate(Files,Schemes,Time_specs,Gas_prop)
