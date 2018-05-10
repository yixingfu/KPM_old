!File=Sparse.f95
!Author=yxfu
!Created=Sun 08 Oct 2017 06:25:23 PM DST
!Last Modified=Thu 10 May 2018 03:05:35 PM DST
! implementation of sparse multiplication and storage.
      subroutine CSRmultVc16(N,NNZ,A,rp,col,v_in, v_out)
          !!! INPUTS
          ! N: dimension of matrix
          ! NNZ: # of nunzero
          ! A: value
          ! rp: row pointer
          ! col: column index
          ! v_in: vector to multiply
          ! v_out: result vector(OUTPUT)
          ! Operation: v_out = M*v_in

          integer*8,intent(in)::N,NNZ
          integer*8,dimension(NNZ),intent(in)::col
          integer*8,dimension(N+1),intent(in)::rp
          complex*16,dimension(NNZ),intent(in)::A
          complex*16,dimension(N),intent(in)::v_in
          complex*16,dimension(N),intent(out)::v_out
          ! local variables
          integer*8::i,j
          complex*16::v_out_i

          v_out = 0d0
          do i=1,N
          v_out_i = 0d0
          do j=rp(i),rp(i+1)-1 
          v_out_i = v_out_i+A(j)*v_in(col(j))
          End do
          v_out(i)=v_out_i
          End do
          return


      End subroutine CSRmultVc16



      subroutine Sparse2Dense(N,NNZ,A,rp,col,Mat_out)
          ! input: same as first few in CSRmultVc16
          integer*8,intent(in)::N,NNZ
          integer*8,dimension(NNZ),intent(in)::col
          integer*8,dimension(N+1),intent(in)::rp
          complex*16,dimension(NNZ),intent(in)::A
          ! output: N*N matrix
          complex*16,dimension(N,N),intent(out)::Mat_out

          ! local variables
          integer*8::i,j

          Mat_out = 0
          do i=1,N
          do j=rp(i),rp(i+1)-1
          Mat_out(i,col(j)) = A(j)
          enddo
          enddo
          return
      end subroutine Sparse2Dense
