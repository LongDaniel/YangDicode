rm -f *~ *.o hos.exe
# For Intel Fortran (ifort)
# mpif77 -r8 -heap-arrays -o hos.exe jonswap_spectrum.f hos.f pressure.f spectral_pack.f check.f smooth_v2.f spectrk2d.f

# For GNU Fortran (gfortran)
mpif77 -fdefault-real-8 -g -o hos.exe jonswap_spectrum.f hos.f pressure.f spectral_pack.f check.f smooth_v2.f spectrk2d.f
