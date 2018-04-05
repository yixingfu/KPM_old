!Created=Wed 13 Dec 2017 01:56:55 PM STD
!Last Modified=Thu 05 Apr 2018 02:38:06 PM EDT

        write(inputfile,'(a,i4.4)')trim(outputfile)//'_',k
      open(11,file=trim(inputfile)//".dat",status="old",access="stream")
      read(11) Nc
      allocate(mu_avg(0:Nc-1),mu2_avg(0:Nc-1))
      read(11) norm_a,norm_b
      read(11) mu_avg,mu2_avg
      close(11)
        if (SetNtilde .eq. 0) then
      Ntilde = 2*Nc
        else 
      Ntilde = SetNtilde
        endif
        write(*,*) "Nc=",Nc,"Forced to ",ForceNc
        

      if (k.eq.RLZmin) then
              ! prepare grid
              allocate(Egrid(1:Ntilde),rho_tot(1:Ntilde), &
                                rho2_tot(1:Ntilde))
              rho_tot = 0
              rho2_tot= 0
      ! This is the same for all. 
      do i=1,Ntilde
              Egrid(i) = i/real(Ntilde)*2d0*outputEmax-outputEmax
      End do
      ! and prepare kernel
        
        include "get_kernel.f90"



      endif
        

      !test
!      write(*,*)Nc
!      write(*,*)norm_a,norm_b
!      write(*,*)mu_avg,mu2_avg
