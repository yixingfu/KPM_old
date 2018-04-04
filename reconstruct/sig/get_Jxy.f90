! Created=Thu 14 Dec 2017 02:49:11 AM STD
! Last Modified=Wed 21 Mar 2018 03:50:57 PM EDT
      ! This file computes Jxy
      ! This is slow
      allocate(Jxy(1:Noutput,1:Noutput))
        allocate(mu_tilde(0:Nc-1,0:Nc-1))
        do m=0,Nc-1
        do n=0,Nc-1
            mu_tilde(m,n) = mu_avg(m,n)*hm(m)*hm(n)*gJ(m)*gJ(n)
        enddo
        enddo

        if (useFFT) then
      write(*,*) "This if the FFT version of Jxy"

        allocate(Jxy_ref(1:Noutput*ref_lev,1:Noutput*ref_lev))
        Jxy_ref=0d0
        Jxy = 0d0
        allocate(Jtemp_in(0:Nc-1))
        allocate(Jtemp_out(0:Noutput*ref_lev-1))
        do i=0,Nc-1
                Jtemp_in(0:Nc-1) = mu_tilde(i,0:Nc-1)
                call cos_trans(Jtemp_in,Jtemp_out,Nc,Noutput*ref_lev)
                Jxy_ref(i+1,1:Noutput*ref_lev) = &
                                Jtemp_out(0:Noutput*ref_lev-1)
        enddo

        do i=0,Noutput*ref_lev-1
                Jtemp_in(0:Nc-1) = Jxy_ref(1:Nc,i+1)
                call cos_trans(Jtemp_in,Jtemp_out,Nc,Noutput*ref_lev)
                Jxy_ref(1:Noutput*ref_lev,i+1) = &
                                Jtemp_out(0:Noutput*ref_lev-1)
        enddo
        do i=1,Noutput*ref_lev
                xi=xygrid_ref(i)
        do j=1,i
                yj=xygrid_ref(j)
                Jxy_ref(i,j) = Jxy_ref(i,j)/((norm_a**2)*&
                        (pi**2)*dsqrt((1d0-xi**2)*(1d0-yj**2)))
                Jxy_ref(j,i) = Jxy_ref(i,j)
        enddo
        enddo
        call interp2d(xygrid_ref,xygrid_ref,Jxy_ref,&
                        xygrid,xygrid,Jxy)
        
        deallocate(Jtemp_in,Jtemp_out,Jxy_ref)
        else 
      write(*,*) "slow version of Jxy"
      
      Jxy = 0
      do i=1,Noutput
              xi=(xygrid(i)-norm_b)/norm_a
      do j=1,i
              yj=(xygrid(j)-norm_b)/norm_a
              Jij = 0
        do m=0,Nc-1
        do n=0,Nc-1
          Jij = Jij+(mu_tilde(m,n)*ChebyT(n,xi)*ChebyT(m,yj))
        End do
        End do
        Jxy(i,j) = Jij/((norm_a**2)*&
                (pi**2)*dsqrt((1-xi**2)*(1-yj**2)))
        Jxy(j,i) = Jxy(i,j)

      End do
      End do

        endif
        
        include "dc_cond.f90"
        DCtot = DCtot + dc_conductivity

        deallocate(mu_avg,mu2_avg,mu_tilde)
        ! after getting Jxy, we need to rescale xygrid. This is done
        ! before saving it.
        

        Jxy_tot = Jxy_tot+dabs(Jxy)
        deallocate(Jxy)
