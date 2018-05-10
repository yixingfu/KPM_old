!Created=Thu 14 Dec 2017 09:07:46 AM STD
!Last Modified=Thu 10 May 2018 03:05:35 PM DST
      ! This is made to be incorporated in other places

      subroutine time_evolve(N,NNZ,A,col,rp,psi,psi_out,t,Ncutoff)
          integer*8,intent(in)::N,NNZ
          complex*16,dimension(NNZ),intent(in)::A
          integer*8,dimension(NNZ),intent(in)::col
          integer*8,dimension(N+1),intent(in)::rp
          complex*16,dimension(N),intent(in)::psi
          real*8,intent(in)::t
          integer,intent(in)::Ncutoff
          complex*16,dimension(N),intent(out)::psi_out
          real*8::Emax,Emin,norm_a,norm_b
          complex*16,dimension(NNZ)::Atilde
          integer::i,j,k
          complex*16,dimension(N)::psi_tmp,psi_p,psi_pp

          call LanczosBound(N,NNZ,A,rp,col,50,Emax,Emin)
          norm_a = dabs(Emax-Emin)/(2d0-0.01d0)
          norm_b = 0.5d0*(Emax+Emin)
          Atilde = (A-norm_b)/norm_a


          psi_out = 0
          psi_pp = psi ! psi0
          call CSRmultVc16(N,NNZ,Atilde,rp,col,psi_pp,psi_p)!psi1

          do i=2,Ncutoff-1
          call CSRmultVc16(N,NNZ,Atilde,rp,col,psi_p,psi_tmp)
          psi_out = 2.0d0*psi_tmp - psi_pp

          psi_pp = psi_p
          psi_p = psi_out
          enddo



          return


      End subroutine time_evolve







