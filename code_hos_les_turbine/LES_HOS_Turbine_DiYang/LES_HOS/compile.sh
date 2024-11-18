# compile using ifort
#mpif77 -r8 -heap-arrays -o hos_les.exe wind_wave_les_wallmodel_turbine_v11_ln_2us.f 

# compile using gfortran
mpif77 -fdefault-real-8  -o hos_les.exe wind_wave_les_wallmodel_turbine_v11_ln_2us.f 
