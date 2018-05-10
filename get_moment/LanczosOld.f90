
!File=Lanczos.f95
!Author=yxfu
!Created=Fri 03 Nov 2017 10:48:07 PM DST
!Last Modified=Thu 10 May 2018 03:05:34 PM DST
      ! Performs Lanczos
      ! Takes CSR matrix
      ! Set Lmax as maximum iteration
      ! output: Emax Emin

      subroutine LanczosBound(N,NNZ,A,rp,col,&
              Lmax,&
              Emax,Emin,&
              rel_err_in)
          ! Input: matrix
          integer*8,intent(in)::N,NNZ
          integer*8,dimension(NNZ),intent(in)::col
          complex*16,dimension(NNZ),intent(in)::A
          integer*8,dimension(N+1),intent(in)::rp
          real*8,intent(in),optional::rel_err_in

          ! Input: Lanczos subspace size
          integer*4,intent(in)::Lmax

          ! Output: Emax, Emin
          real*8,intent(out)::Emax,Emin

          real*8::rel_error

          ! local var
          ! iteration
          integer::i,j,k ! index is still 32 bit
          real*8::EmaxTEMP,EminTEMP,EmaxL,EminL
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

          if (present(rel_err_in)) then
              rel_error = rel_err_in
          else 
              rel_error = 1d-10
          endif



          ! First generate random vector and normalize
          call random_number(phiL0_r)
          phiL0 = phiL0_r
          normL = dot_product(phiL0,phiL0)
          normL = dsqrt(normL)
          phiL0 = phiL0/normL

          ! Then generate phiL1
          ! phiL1 = H|phiL0>-<phiL0|phiL1>|phiL0>
          call CSRmultVc16(N,NNZ,A,rp,col,phiL0,phiL1)
          alpha(1) = dot_product(phiL0,phiL1)
          phiL1 = phiL1 - alpha(1)*phiL0
          ! and normalize
          normL = dot_product(phiL1, phiL1)
          normL = dsqrt(normL)
          beta(1) = normL
          phiL1 = phiL1/normL


          ! Initiate Emin Emax
          Emin = 0.0d0
          EminTEMP = 1d70
          Emax = 0.0d0
          EmaxTEMP = -1d70

          ! up to Lmax steps
          do j = 1,Lmax
          ! get phiL2 from phiL1
          ! |phiLj+1>=H|phiLj>-a(j)|phiLj>-b(j-1)|phiLj-1>
          call CSRmultVc16(N,NNZ,A,rp,col,phiL1,phiL2)
          alpha(j+1) = dot_product(phiL1,phiL2)
          phiL2 = phiL2 - alpha(j+1)*phiL1 - beta(j)*phiL0
          ! and normalize
          normL = dot_product(phiL2,phiL2)
          normL = dsqrt(normL)
          beta(j+1) = normL
          phiL2 = phiL2/normL

          ! update for next iteration
          phiL0 = phiL1
          phiL1 = phiL2

          ! now we have constructed subspace up to j+1
          ! Diagonalize.


          ! allocate needed array and reset
          allocate(aL(j+1,j+1))
          allocate(wL(j+1))
          allocate(workL(2*j+1))
          allocate(rworkL(3*j+1))
          ierr = 0
          aL = 0.d0
          wL = 0d0
          workL=0d0
          rworkL=0d0

          do i=1,j+1
          ! for every dimension in Lanczos up to j+1,
          ! diagonal term is alpha
          aL(i,i)=(alpha(i))
          End do

          do i=1,j
          ! for every dimension up to j, next to diagonal 
          ! term is beta
          aL(i,i+1)=beta(i)
          aL(i+1,i)=beta(i) ! this line not needed because 
          ! only upper triangle used,
          End do


          ! now aL is H corresponding to Lanczos subspace.
          ! diagonalize it
          ! use zheev for comp*16,real*8
          ! 'N' for only computing eigenvalue
          ! 'U' for udsing upper 
          ! j+1 for size of aL
          ! j+1 for LDA
          ! wL for result
          ! workL for dim of wL
          !write(*,*) aL
          call zheev('N','U',j+1,aL,j+1,wL,workL,2*j+1,rworkL,info)
          !write(*,*) "zheev result."
          !write(*,*) "info: ",info
          !write(*,*) wL

          ! assign minimum to EminTEMP
          EminL = minval(wL)
          EmaxL = maxval(wL)

          ! now the arrays can be deallocated
          deallocate(aL,wL,workL,rworkL)

          ! check convergence -> if converge, goto 
          if ((dabs(EminTEMP-EminL) .lt. rel_error) .and. &
              (dabs(EmaxTEMP-EmaxL) .lt. rel_error)) then 
              Emin=EminL
              Emax=EmaxL
              write(*,*)"Good convergence" ,Emin,Emax
              goto 2223
          end if
          ! if not, still update
          EminTEMP = EminL
          EmaxTEMP = EmaxL


          ! (end) up to Lmax steps
          End do
          Emin = EminTEMP
          Emax = EmaxTEMP
          write(*,*)'not good convergence', Emin,Emax

2223    continue

          return
      End subroutine LanczosBound


