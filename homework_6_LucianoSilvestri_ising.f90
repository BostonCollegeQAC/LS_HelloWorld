!=========================================================!
!			MODULE				  !
!=========================================================!
 
 module systemvariables

 implicit none
 
 integer		:: L, volume, Nsizes
 real(8) 		:: pflips(-2:2,-4:4)
 integer	 	:: choices(0:1,-1:1) ! see comment below for choices array
 integer,allocatable	:: spin(:), sizes(:)

 end module systemvariables
 
!**********************************************************!
!			MAIN PROGRAM			   !
!**********************************************************!
 
 program ising2d
 
 use systemvariables
 implicit none
  
 integer	:: nt, bins, binsteps, i, j, k
 real(8)	:: tmax, tmin, dt, temp, time0, time1
 real(8) 	:: energy,energy2,magnet,magnet2,nspins,nspins2
 common/obs/energy,energy2,magnet,magnet2,nspins,nspins2
 
!--------------------------------------------------------------!
! choices array is a matrix containing the possible values     !
! each spin can take when it has to flip. 		       !
! the column is indicated by the spin value {-1,0,1}           !
! the row by the {0,1} where 0 and 1 are the values the random !
! number can take                                              !
!--------------------------------------------------------------!

 choices(0,-1)=1; choices(0,0)=-1; choices(0,1)=0
 choices(1,-1)=0; choices(1,0)=1; choices(1,1)=-1

 write(*,*) 'How many sizes do you want to calculate ? We start from L=4'
 read(*,*) Nsizes
 
 allocate(sizes(Nsizes))
 do i=1,Nsizes
    sizes(i)=2**(i+1)
 enddo

 open(10,file='read.in', status='old')
 read(10,*) nt, tmax, dt
 read(10,*) bins, binsteps
 close(10)

 tmin=tmax-nt*dt
 write(*,*)tmin

 call cpu_time(time0)
 do k=1,Nsizes
    L=sizes(k);  volume=L*L
    temp=tmax
    do while (temp.ge.tmin)
       call initialize(temp)
       do i=1,binsteps
          call mcstep
       enddo
       do j=1,bins
          call resetdatasums
          do i=1,binsteps
             call mcstep
             call measurement
          end do
          call writebindata(binsteps,temp,j)
       end do
       temp=temp-dt
       deallocate(spin)
      
     end do	! while loop on temperatures
     call cpu_time(time1)
    23 format ('L = ', i3, ' Time elapsed = ', f10.2, ' sec')
    write(*,23) L, time1-time0
 end do ! loop on sizes
 
 deallocate(sizes)

 end program ising2d

!**********************************************************!

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!		SUBROUTINES & FUNCTIONS 		   !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

 subroutine initialize(t)
 
 use systemvariables
 implicit none

 integer	:: i,j
 real(8)	:: t,rnd
 
 pflips=0
 do j=-4,4
    do i=-2,2
       pflips(i,j)=exp(-dble(i)*dble(j)/t)
    enddo
 enddo

 !56 format(8f12.7)
! do i=-2,2
! write(*,56) pflips(i,:)
! enddo
 
 call initran(1)
 
 allocate(spin(0:volume-1))
 do j=0,volume-1
    spin(j)=int(3.*rnd())-1
 end do

 end subroutine initialize

!===========================================================! 

 subroutine mcstep
 
 use systemvariables
 
 implicit none
 
 integer	 :: i,s,x,y,s1,s2,s3,s4,c,cprime
 real(8)	 :: rnd

 do i=1,volume
   s=int(rnd()*volume)
   x=mod(s,L); y=s/L
   s1=spin(mod(x+1,L)+y*L)
   s2=spin(x+mod(y+1,L)*L)
   s3=spin(mod(x-1+L,L)+y*L)
   s4=spin(x+mod(y-1+L,L)*L)
   c=spin(s)
   cprime=choices(int(2.*rnd()),spin(s))
   if (rnd().lt.pflips(c-cprime,s1+s2+s3+s4)) spin(s)=cprime
 end do
 
 end subroutine mcstep
 
!============================================================!  
 
 subroutine measurement
 
 use systemvariables
 implicit none
  
 integer	:: e,m,ns,x,y,s
 real(8) 	:: energy,energy2,magnet,magnet2,nspins,nspins2
 common/obs/energy,energy2,magnet,magnet2,nspins,nspins2
 
 e=0; m=0; ns=0 
 do s=0,volume-1
    x=mod(s,L); y=s/L
       e=e-spin(s)*(spin(mod(x+1,L)+y*L)+spin(x+mod(y+1,L)*L))	! ENERGY
       ns=ns+abs(spin(s))					! NUMBER OF SPINS
 end do
 m=sum(spin) 
 m=abs(m)
 energy=energy+dble(e)
 energy2=energy2+dble(e)**2
 magnet=magnet+dble(m)
 magnet2=magnet2+dble(m)**2
 nspins=nspins+dble(ns)
 nspins2=nspins2+dble(ns)**2
 
 end subroutine measurement

!===========================================================! 

 subroutine resetdatasums

 implicit none

 real(8) 	:: energy,energy2,magnet,magnet2,nspins,nspins2
 common/obs/energy,energy2,magnet,magnet2,nspins,nspins2

 energy=0.d0; energy2=0.d0; magnet=0.d0; magnet2=0.d0; nspins=0.d0; nspins2=0.d0

 end subroutine resetdatasums

!===========================================================! 

 real(8) function rnd()
!----------------------------------------------!
! 64-bit congruental generator                 !
! iran64=oran64*2862933555777941757+1013904243 !
!----------------------------------------------!
 implicit none

 real(8)    :: dmu64
 integer(8) :: ran64,mul64,add64
 common/bran64/dmu64,ran64,mul64,add64

 ran64=ran64*mul64+add64
 rnd=0.5d0+dmu64*dble(ran64)

 end function rnd

!===========================================================! 

 subroutine initran(w)

 implicit none

 integer(8) :: irmax
 integer(4) :: w,nb,b

 real(8)    :: dmu64
 integer(8) :: ran64,mul64,add64
 common/bran64/dmu64,ran64,mul64,add64
      
 irmax=2_8**31
 irmax=2*(irmax**2-1)+1
 mul64=2862933555777941757_8
 add64=1013904243
 dmu64=0.5d0/dble(irmax)

 open(10,file='seed.in',status='old')
 read(10,*)ran64
 close(10)
 if (w.ne.0) then
    open(10,file='seed.in',status='unknown')
    write(10,*)abs((ran64*mul64)/5+5265361)
    close(10)
 endif

 end subroutine initran

!=============================================================! 
 
 subroutine writebindata(steps,temp,bin)
 
 use systemvariables
 implicit none
 
 integer	:: steps,bin
 character(len=10) :: ci
 real(8)	:: temp
 real(8) 	:: energy,energy2,magnet,magnet2,nspins,nspins2
 common/obs/energy,energy2,magnet,magnet2,nspins,nspins2
 
 write(ci,'(I5)') L
 ci=adjustl(ci)

 open(1,file='bindata_L'//trim(ci)//'.dat',status='unknown',position='append')
 energy=energy/(dble(steps)*dble(volume))
 energy2=energy2/(dble(steps)*dble(volume)**2)
 magnet=magnet/(dble(steps)*dble(volume))
 magnet2=magnet2/(dble(steps)*dble(volume)**2)
 nspins=nspins/(dble(steps)*dble(volume))
 nspins2=nspins2/(dble(steps)*dble(volume)**2)

 write(1,1) temp, bin, energy, energy2, magnet, magnet2, nspins, nspins2
 1 format(f8.4,'   ',i3,'   ',6f12.7)
 close(1)
 
 end subroutine writebindata
