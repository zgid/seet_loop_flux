
set(LIB_SOURCES  0_gf2_modu.f90          fetch.f90      integral1.f90  mu_dens.f90       print_g_se.f90      read_input.f90        sigma.f90          transform_ints.f90
aux_library.f90         fockgauss.f90  integrals.f90  pbc_nao.f90       print_se_real.f90   read_proc_k.f90       sigma_pbc.f90
compare_sigma_fock.f90  ft.f90         legendre.f90   pbc_sao.f90       print_sigma_om.f90  retired_routines.f90  sigma_pbc_rev.f90
energy.f90              inigrids.f90   mlpn.f90       pbc_specific.f90  read_gl.f90         seet.f90              timest.f90)
             
add_library(gf2-common-lib ${LIB_SOURCES})

#target_lib_library(gf2-common-lib )
