 program calculations
 
 implicit none
 
 integer	:: Nbins, bin, i, ctrl, Ntemps, Nsizes, k
 real(8)	:: temp,e,e2,m,m2,n,n2,avg_e,avg_m,avg_n,binder_ratio,c_v,chi,kappa
 real(8)	:: sigma_e, sigma_c, sigma_m, sigma_chi, sigma_n, sigma_k, sigma_binder, volume
 integer,allocatable	:: sizes(:)
 character(len=10)	:: ci

 write(*,*) 'How many bins?'; read(*,*) Nbins
 write(*,*) 'How many temperatures?'; read(*,*) Ntemps
 write(*,*) 'How many sizes did you do?'; read(*,*) Nsizes
  
 allocate(sizes(Nsizes))
 do i=1,Nsizes
    sizes(i)=2**(i+1)
 enddo

 do k=1,Nsizes
    write(ci,'(I5)') sizes(k)
    ci=adjustl(ci)
    volume=sizes(k)*sizes(k)

    open(2,file='bindata_L'//trim(ci)//'.dat',status='old')
  
    open(10,file='e_L'//trim(ci)//'.dat',status='replace')
    open(20,file='m_L'//trim(ci)//'.dat',status='replace')
    open(30,file='n_L'//trim(ci)//'.dat',status='replace')
    open(40,file='q_L'//trim(ci)//'.dat',status='replace')

    do i=1,Ntemps
       ctrl=0;
       avg_e=0; c_v=0
       avg_m=0; chi=0
       avg_n=0; kappa=0
       binder_ratio=0

       sigma_e=0; sigma_c=0
       sigma_m=0; sigma_chi=0
       sigma_n=0; sigma_k=0
       sigma_binder=0
   
       do while (ctrl.eq.0)
          read(2,*) temp, bin, e, e2, m, m2, n, n2
          avg_e=avg_e+e; sigma_e=sigma_e+e2
          avg_m=avg_m+m; sigma_m=sigma_m+m2
          avg_n=avg_n+n; sigma_n=sigma_n+n2
       
          c_v=c_v+(e2-e**2)/(temp**2)	
          kappa=kappa+(n2-n**2)/temp
          chi=chi+(m2-m**2)/temp
          binder_ratio=binder_ratio+m2/m**2
       
          sigma_c=sigma_c+((e**2-e2)/temp**2)**2
          sigma_k=sigma_k+((n**2-n2)/temp)**2
          sigma_chi=sigma_chi+((m**2-m2)/temp)**2
          sigma_binder=sigma_binder+(m2/m**2)**2

          if (bin.eq.Nbins) then
	     avg_e=avg_e/dble(Nbins); avg_m=avg_m/dble(Nbins); avg_n=avg_n/dble(Nbins)
	     c_v=c_v*dble(volume)/dble(Nbins); chi=chi*dble(volume)/dble(Nbins); kappa=kappa*dble(volume)/dble(Nbins)
             binder_ratio=binder_ratio/dble(Nbins)
             
             sigma_e=(sigma_e-dble(Nbins)*avg_e**2)/dble(Nbins*(Nbins-1)); sigma_e=sqrt(sigma_e)
             sigma_m=(sigma_m-dble(Nbins)*avg_m**2)/dble(Nbins*(Nbins-1)); sigma_m=sqrt(sigma_m)
             sigma_n=(sigma_n-dble(Nbins)*avg_n**2)/dble(Nbins*(Nbins-1)); sigma_n=sqrt(sigma_n)
             
             sigma_c=(sigma_c*dble(volume)**2-dble(Nbins)*c_v**2)/dble(Nbins*(Nbins-1)); sigma_c=sqrt(sigma_c)
             sigma_k=(sigma_k*dble(volume)**2-dble(Nbins)*kappa**2)/dble(Nbins*(Nbins-1)); sigma_k=sqrt(sigma_k)
             sigma_chi=(sigma_chi*dble(volume)**2-dble(Nbins)*chi**2)/dble(Nbins*(Nbins-1)); sigma_chi=sqrt(sigma_chi)

             sigma_binder=(sigma_binder-dble(Nbins)*binder_ratio**2)/dble(Nbins*(Nbins-1)); sigma_binder=sqrt(sigma_binder)

    
             write(10,*) temp, avg_e, sigma_e, c_v, sigma_c
             write(20,*) temp, avg_m, sigma_m, chi, sigma_chi
	     write(30,*) temp, avg_n, sigma_n, kappa, sigma_k
             write(40,*) temp, binder_ratio, sigma_binder
   
             ctrl=1
             end if
         end do 	! while loop to collect/calculate avgs and sigmas
     end do 	! loop on Ntemps
     close(2); close(10); close(20); close(30); close(40)
 end do 	! loop on sizes

 deallocate(sizes)
 end program calculations
