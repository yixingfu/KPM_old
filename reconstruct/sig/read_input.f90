!Created=Wed 13 Dec 2017 01:56:55 PM STD
!Last Modified=Mon 26 Mar 2018 06:15:19 PM EDT

        write(inputfile,'(a,i4.4)')trim(outputfile)//'_',k
      open(11,file=trim(inputfile)//".dat",&
                iostat=ierr,status="old",access="stream")
     read(11,iostat=ierr) Nc
        write(*,*)"Nc=",Nc
      allocate(mu_avg(0:Nc-1,0:Nc-1),mu2_avg(0:Nc-1,0:Nc-1))
      read(11,iostat=ierr) norm_a,norm_b
      read(11,iostat=ierr) mu_avg,mu2_avg
      close(11,iostat=ierr)
         if (ierr .ne. 0) then 
                badfiles = badfiles+1
                close(11)
                cycle
        endif
        write(*,*)"Done reading"
      Ntilde = 2*Nc
        if (Noutput.eq.0) then
                Noutput=Ntilde
        endif

      if (k.eq.RLZmin) then
              ! prepare grid
              allocate(xygrid(1:Noutput),Jxy_tot(1:Noutput,1:Noutput))
                Jxy_tot = 0
      ! This is the same for all. 
        !but different for FFT or non FFT
        if (useFFT) then 
              call cos_trans_grid(Noutput,xygrid)
                allocate(xygrid_ref(1:Noutput*ref_lev))
              call cos_trans_grid(Noutput*ref_lev,xygrid_ref)

         else
      do i=1,Noutput
              xygrid(i) = 0.99d0*(i/real(Noutput)*2d0-1d0)*(EmaxNONFFT)! xy grid is final. Discarding a boundary.
      End do
        endif
        write(*,*)"grid prepared"
      ! and prepare kernel
        
        include "get_kernel.f90"

        ! read auxiliary file for DC conductivity
!        allocate(intgammamn(0:Nc-1,0:Nc-1))
!      open(18,file="/scratch/yf160/intgammamn.dat",&
!                form="unformatted",access="stream")
!      read(18)intgammamn
!      close(18)
        DCtot=0


      endif
        
        
        

      !test
!      write(*,*)Nc
!      write(*,*)norm_a,norm_b
!      write(*,*)mu_avg,mu2_avg
