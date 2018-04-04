! Created=Mon 22 Jan 2018 12:49:18 PM EST
! Last Modified=Wed 24 Jan 2018 09:02:39 AM EST

! This file includes function for cos transformation.
! The transformation computes 
! \gamma_k = \sum_{n=0}^{N-1} cos(n \phi_k) \mu_n
! Where gamma_k is output, mu_k is input. The input has
! moments up to dimension N, and should be already prepared
! with h and/or g, the kernels, before feeding into this.

        subroutine cos_trans(mu_n,gamma_k,Nc,Ntilde)
                integer,intent(in)::Nc,Ntilde
                real*8,intent(in),dimension(0:Nc-1)::mu_n
                real*8,intent(out),dimension(0:Ntilde-1)::gamma_k
                complex*16,dimension(0:Ntilde-1)::lambda
                integer::n,j,k
                complex*16::i
                type(DFTI_DESCRIPTOR),POINTER::My_Desc1_Handle
                integer::status
                real*8,parameter::pi=3.1415926535897932384626433832795d0

                i = dcmplx(0d0,1d0)
                lambda=0d0
                if (Ntilde .le. Nc) then
                do n=0,Ntilde-1
                        lambda(n) = mu_n(n) * zexp(i*pi*n/2d0/Ntilde)
                enddo 
                else ! Ntilde .gt Nc
                do n=0,Nc-1
                        lambda(n) = mu_n(n) * zexp(i*pi*n/2d0/Ntilde)
                enddo
                endif

                

                ! Then do standard FFT
                status = DftiCreateDescriptor(My_Desc1_Handle,&
                        DFTI_DOUBLE,DFTI_COMPLEX,&
                        1,Ntilde)
                status = DftiCommitDescriptor(My_Desc1_Handle)
                status = DftiComputeForward(My_Desc1_Handle,lambda)
                status = DftiFreeDescriptor(My_Desc1_Handle)
                
                ! put into order
                do j=0,Ntilde/2-1
                        gamma_k(2*j) = real(lambda(j))
                        gamma_k(2*j+1) = real(lambda(Ntilde-1-j))
                enddo

                return
        end subroutine cos_trans


        
        subroutine cos_trans_grid(Ntilde, grid_k)
                integer,intent(in)::Ntilde
                real*8,intent(out),dimension(0:Ntilde-1)::grid_k
                integer::k
                real*8,parameter::pi=3.1415926535897932384626433832795d0
                
                do k=0,Ntilde-1
                        grid_k(k) = dcos(pi*(k+0.5d0)/Ntilde)
                enddo

                return
        end subroutine cos_trans_grid
