subroutine calc_gl_from_gtau(nts,uniform,Gt)
  use gf2, only:nao,max_points_mesh,power_mesh,power_mesh_small,power_mesh_small_int,uniform_tau, power_tau,Gl,nleg,beta,ncell
  implicit none
  integer::nts,uniform,ic
  double precision:: x,h,sqrtl,x1,x2,x3
  double precision,dimension(:),allocatable:: pn1,pn2,pn3,pd
  double precision,dimension(nao,nao,ncell,nts) :: Gt
  integer:: ntaul,ntaur,n_int
  integer:: ntau1,ntau2,ntau3
  integer:: l,p,t,i,j

  n_int = 2*power_tau + 2
  

  allocate(pn1(0:nleg))
  allocate(pn2(0:nleg))
  allocate(pn3(0:nleg))
  allocate(pd(0:nleg))

  Gl(:,:,:,:)=0.0d0
  do p=1,n_int
     ntaul = uniform*(p-1)
     ntaur = uniform*(p)
     h=power_mesh_small(ntaul+2)-power_mesh_small(ntaul+1)
     do t=1,uniform/2
        ntau1 = ntaul + 2*t-1
        x1=(power_mesh_small(ntau1)/beta-0.5d0)*2.0d0
        call lpn(nleg,x1,pn1(0),pd(0))
        ntau2 = ntaul + 2*t
        x2=(power_mesh_small(ntau2)/beta-0.5d0)*2.0d0
        call lpn(nleg,x2,pn2(0),pd(0))
        ntau3 = ntaul + 2*t+1
        x3=(power_mesh_small(ntau3)/beta-0.5d0)*2.0d0
        call lpn(nleg,x3,pn3(0),pd(0))
        do l=0,nleg
           sqrtl=sqrt(2.0d0*l+1.0d0)
           Gl(:,:,:,l)=Gl(:,:,:,l)+h/3.0d0*sqrtl*(Gt(:,:,:,ntau1)*pn1(l)+4.0d0*Gt(:,:,:,ntau2)*pn2(l)+Gt(:,:,:,ntau3)*pn3(l))
        end do
     end do
  end do
  

   deallocate(pn1)
   deallocate(pn2)
   deallocate(pn3)
   deallocate(pd)
end subroutine calc_gl_from_gtau
