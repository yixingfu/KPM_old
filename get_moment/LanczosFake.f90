!File=Lanczos.f95
!Author=yxfu
!Created=Fri 03 Nov 2017 10:48:07 PM DST
!Last Modified=Thu 10 May 2018 03:05:34 PM DST
      !Fake Lanczos

      subroutine LanczosBound(N,NNZ,A,rp,col,Lmax,Emax,Emin)
          integer*8::N,NNZ
          complex*16,dimension(NNZ)::A
          integer*8,dimension(NNZ)::col
          integer*8,dimension(N+1)::rp
          integer::Lmax
          real*8::Emax,Emin
          real*8,parameter::rel_error=1.0d-4

          integer::i,j,k ! index is still 32 bit
          real*8::EmaxTEMP,EminTEMP
          ! constructing Lanczos subspace
          complex*16,dimension(N)::phiL0,phiL1,phiL2
          real*8,dimension(N)::phiL0_r
          real*8::normL
          complex*16,dimension(Lmax+1)::alpha, beta
          ! note: they are real. 
          ! Diagonalizing 
          complex*16,dimension(:,:),allocatable::aL
          real*8,dimension(:),allocatable::wL,rworkL
          complex*16,dimension(:),allocatable::workL
          integer::info

          ! misc
          integer::ierr


          write(*,*) "Fake Lanczos for testing"
          Emin = -3
          Emax = 3


          return
      End subroutine LanczosBound



