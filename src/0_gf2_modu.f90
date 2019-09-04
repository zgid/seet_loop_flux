module gf2
  implicit none 
  save
  integer :: nk, nao, ncell, ncellgf2, ncellden, nel_cell
  integer :: uniform_tau, power_tau, max_points_mesh, nuniform
  integer :: nleg
  integer :: iwmax
  integer :: itermax
  integer :: n_threads
  integer :: hf_temp
  integer :: debug_p
  integer :: rst
  integer :: rst_old_sig


  logical,allocatable, dimension(:) :: vertex

  double precision :: beta, mu

  integer,allocatable, dimension(:,:) :: indices_cell, diffcell, sumcell
  integer,allocatable, dimension(:,:,:) :: i_cell_p
  integer,allocatable, dimension(:) :: n_cell_p

  double precision, allocatable, dimension(:) :: k_weights

  double precision, allocatable, dimension(:) :: power_mesh, power_mesh_small
  integer,allocatable, dimension(:) :: power_mesh_small_int

  complex*16, allocatable, dimension(:,:,:) :: S_k, F_k
  complex*16, allocatable, dimension(:,:) :: four_coef, tnl
  complex*16, allocatable, dimension(:) :: omegas

  complex*16, allocatable, dimension(:,:,:) :: Sigma
 
  double precision, allocatable, dimension(:,:,:,:) :: Gr_tau 
  double precision, allocatable, dimension(:,:,:,:) :: Gl
  
  double precision, allocatable, dimension(:,:,:,:) :: vi,vi_full
  double precision, allocatable, dimension(:) :: vi_s
  double precision, allocatable, dimension(:,:) :: vi_snc, vi_sne
  integer, allocatable, dimension(:,:,:) :: i_vi_snc, i_vi_sne 

  double precision, allocatable, dimension(:) :: vi_sncl, vi_snel 
  integer, allocatable, dimension(:,:) :: i_vi_sncl, i_vi_snel

 
  integer, allocatable, dimension(:) :: i_len_c, i_len_e
  integer, allocatable, dimension(:,:) :: v_ind, v_ind_2
  integer, allocatable, dimension(:) :: v_ind_last, v_ind_1
  integer, allocatable, dimension(:) :: n_vi_c, n_vi_e
  integer :: np, ind
 
  double precision, allocatable, dimension(:,:,:,:) :: Sigma_tau
  complex*16, allocatable, dimension(:,:,:,:) :: Sigma_rw, Sigma_kw, Sigma_kt

  integer :: nc, logm ! maybe these are temporary
  integer :: ncs

  double precision :: E_thr ! energy threshold for convergence
  double precision :: damp, dampd  ! self-energy and density damping

  double precision, allocatable, dimension(:) :: sig11
  double precision, allocatable, dimension(:,:,:,:,:) :: Sigma_tau_tmp

  double precision :: nkpw
  double precision, allocatable, dimension(:) :: w_func
  integer seet_gf2_iter  



end module gf2

module seet
  implicit none
  logical ortho_basis,dm_basis,fock_reeval,HF_SEET
  integer cells_seet
  integer niter_seet 
  double precision thr_seet
  integer num_imp   !number of impurities
  integer, dimension(:), allocatable::num_orb_in_imps !number of orbitals in every impurity
  integer, dimension(:,:), allocatable::imp_orb_numbers !orbtial numbers for every impurity (imp, orbital numbers)
  integer, dimension(:,:), allocatable::imp_bath_numbers !bath numbers for every impurity orbital within an impurity (imp, orbital numbers)
  integer max_imp_orb !number of orb in the biggest impurity
  double precision :: damps, dampsinf
end module seet
