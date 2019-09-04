subroutine read_params
  use gf2, only : nk, nao, ncell, ncellgf2, ncellden, nel_cell, seet_gf2_iter
  use gf2, only : nleg
  use gf2, only : iwmax
  use gf2, only : beta
  use gf2, only : uniform_tau, power_tau, nuniform
  use gf2, only : itermax
  use gf2, only : n_threads
  use gf2, only : E_thr
  use gf2, only : damp
  use gf2, only : dampd
  use gf2, only : hf_temp
  use gf2, only : debug_p
  use gf2, only : rst,rst_old_sig

  implicit none
  logical dummy_L
  integer ip, dummy_i
  double precision dummy_r
    
  logical exist

  ip = 77708

  inquire(file="gf2.inp", exist=exist) 
  If(.not.exist) then
    call timest('GF2 input file is missing. Calculation aborted')
    stop
  Endif
  open(ip,file='gf2.inp',status='old')
  call fetch_i(nel_cell,    'nel_cell   ',ip)
  call fetch_i(nao,         'nao        ',ip)
  call fetch_i(nk,          'nk         ',ip)
  call fetch_r(beta,        'beta       ',ip)
  call fetch_i(ncell,       'ncell      ',ip)
  call fetch_i(ncellgf2,       'ncellgf2   ',ip)
  call fetch_i(uniform_tau, 'uniform    ',ip)
  call fetch_i(nuniform, 'nuniform    ',ip)
  call fetch_i(power_tau,   'power      ',ip)
  call fetch_i(nleg,        'nleg       ',ip)
  call fetch_i(iwmax,       'iwmax      ',ip)
  call fetch_i(itermax,       'itermax      ',ip)
  call fetch_i(n_threads,       'n_threads      ',ip)
  call fetch_r(E_thr,        'E_thr       ',ip)
  call fetch_r(damp,        'damp       ',ip)
  call fetch_r(dampd,        'dampd      ',ip)
  call fetch_i(hf_temp,       'hf_temp    ',ip)
  call fetch_i(debug_p,       'debug_p    ',ip)
  call fetch_i(rst,         'rst        ',ip)
  call fetch_i(seet_gf2_iter,         'seet_gf2_iter        ',ip)
  call fetch_i(rst_old_sig,         'rst_old_sig        ',ip)
  
  call flush(6)
  close(ip)
 
end subroutine read_params
