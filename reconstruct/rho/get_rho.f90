! Created=Wed 13 Dec 2017 03:05:25 PM STD
! Last Modified=Thu 05 Apr 2018 05:25:01 PM EDT
      ! 
      allocate(rho(1:Ntilde))
      ! it need to be scaled here, and scaled back later
        allocate(Egrid_t(1:Ntilde))
      Egrid_t = (Egrid-norm_b)/norm_a


      if (useFFT) then
              write(*,*) "Not implemented yet"
      else 
        if (Der .eq. 0) then
              do i=1,Ntilde
                !1/(pi sqrt(1-x^2)) (g0 mu0+2 sum (g_n mu_n Tn(x)))
                      Ei = Egrid_t(i)
                      rhoi = gJ(0)*mu_avg(0)
                      do j=1,ForceNc-1
                              rhoi = rhoi+&
                2d0*gJ(j)*mu_avg(j)*ChebyT(j,Ei)
                      End do
                      rho(i) = rhoi/(norm_a*pi*dsqrt(1d0-Ei*Ei))
              End do
        else if (Der .eq. 1) then
              do i=1,Ntilde
                        Ei = Egrid_t(i)
                        rhoi = gJ(0)*mu_avg(0)*(1d0/(1d0-Ei*Ei))
                ! T_0(x)=1,U_-1(x)=0
                        do j=1,ForceNc-1
                                rhoi = rhoi+&
                2d0*gJ(j)*mu_avg(j)*(ChebyT(j,Ei)/(1d0-Ei*Ei) &
                                        + ChebyU(j-1,Ei)*j)
                        enddo
                      rho(i) = rhoi/(norm_a*pi*dsqrt(1d0-Ei*Ei))
              enddo
        endif
    
      endif
      rho_tot = rho_tot+(rho)
      rho2_tot = rho2_tot+rho*rho
      deallocate(rho,Egrid_t)
      deallocate(mu_avg,mu2_avg)! no longer used

