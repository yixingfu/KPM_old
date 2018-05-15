! Created =Wed 13 Dec 2017 12:14:32 PM STD
! Last Modified=Thu 10 May 2018 03:02:02 PM DST
      ! This file computes mu

      allocate(mu_tot(0:Nc-1),mu2_tot(0:Nc-1))
      allocate(mu_avg(0:Nc-1),mu2_avg(0:Nc-1))
      allocate(mu(0:Nc-1))
      allocate(psi0R(N),psi0(N),psi1(N),psi_tmp(N))
      allocate(psi_j(N),psi_j_p(N),psi_j_pp(N))
      mu_tot = 0
      mu2_tot = 0
!        call ResetRandSeed()
      idum=-(my_id+1)*64
!      randomtest1 = ran2(idum)
!      randomtest2 = ran2(idum)
!      write(*,*) randomtest1,randomtest2

      do i=1,Rep
      ! for each repetition
      ! first, create random vector psi0(random phase)
      do ipsi=1,N
      psirand = ran2(idum)*pi*2d0
      psi0(ipsi) = dcos(psirand)+III*dsin(psirand)
      enddo
      ! psi0 gives mu0
      mu = 0
      mu(0) = dot_product(psi0,psi0)

      ! Then calculate psi1
      call CSRmultVc16(N,NNZ,A,rp,col,psi0,psi1)
      ! psi1 gives mu1
      mu(1) = dot_product(psi0,psi1)

      psi_j_pp = psi0
      psi_j_p  = psi1
      j0 = 2

      do j=j0,Nc-1
      ! \psi_j = 2H \psi_{j-1} - \psi_{j-2}
      call CSRmultVc16(N,NNZ,A,rp,col,psi_j_p,psi_tmp)
      psi_j = 2d0*psi_tmp-psi_j_pp

      ! mu(j)
      mu(j) = dot_product(psi0,psi_j)

      ! renew psi_j_pp and psi_j_p
      psi_j_pp = psi_j_p
      psi_j_p = psi_j
      enddo
      mu = mu/N
      if (task .eq. RHODER) then
          write(*,*) "THIS IS ABANDONED"
!                allocate(orig_mu(0:Nc-1))
!                orig_mu = mu
!                mu = 0d0
!                ! For even Uk
!                tempTsum = 0d0
!                do iT=0,Nc-3,2
!                        tempTsum = tempTsum+orig_mu(iT)
!                        mu(iT+1) = (iT+1)*(2d0*tempTsum-1d0)
!                enddo
!                if (modulo(Nc,2).eq.0) then
!                        tempTsum = tempTsum+orig_mu(Nc-2)
!                        mu(Nc-1) = (Nc-1)*(2d0*tempTsum-1d0)
!                endif
!
!                ! For odd Uk
!                tempTsum = 0d0
!                do iT=1,Nc-3,2
!                        tempTsum = tempTsum+orig_mu(iT)
!                        mu(iT+1) = (iT+1)*(2d0*tempTsum)
!                enddo
!                if (modulo(Nc,2).eq.1) then
!                        tempTsum = tempTsum+orig_mu(Nc-2)
!                        mu(Nc-1) = (Nc-1)*(2d0*tempTsum)
!                endif
!
!                deallocate(orig_mu)
      endif

!        write(*,*)mu(0:5)
      mu_tot = mu_tot+mu
      mu2_tot = mu2_tot+mu**2
      enddo
      mu_avg = mu_tot/(Rep)
      mu2_avg = mu2_tot/(Rep)
      ! divide by N as required by DoS
      ! This is because the delta function rep of rho has 1/N,
      ! dividing by total # of states.
      deallocate(mu_tot,mu2_tot,mu)
      deallocate(psi0R,psi0,psi1)
      deallocate(psi_j,psi_j_p,psi_j_pp,psi_tmp)
