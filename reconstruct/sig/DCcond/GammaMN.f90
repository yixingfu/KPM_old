! Created=Tue 16 Jan 2018 04:24:11 AM STD
! Last Modified=Thu 25 Jan 2018 08:24:54 PM EST
      ! This file computes gamma_mn as prescribed in PRL Garcia,
      ! Covaci & Rappoport 

      ! We also prepare a int_gamma_mn that includes the integral part. 
      ! For this function we specify m,n and total number of points 
      ! for computing numerical integral

      complex*16 function gamma_mn(m,n,eps)
              real*8,intent(in)::eps
              integer,intent(in)::m,n
              complex*16::i
              
              i = dcmplx(0d0,1d0)

              gamma_mn = ((eps-i*n*dsqrt(1d0-eps**2d0))&
                      * zexp( i*n*dacos(eps))*ChebyT(m,eps))&
                      + ((eps+i*m*dsqrt(1d0-eps**2d0))&
                      * zexp(-i*m*dacos(eps))*ChebyT(n,eps))
      end function gamma_mn

      complex*16 function int_gamma_mn(m,n,mu,beta,pts)
              ! NOTE: mu should be rescaled accordingly
              integer,intent(in)::m,n,pts
              real*8,intent(in)::mu,beta
              real*8::eps
              complex*16::result_tmp,result_step
              integer::i
              result_tmp=0d0
              do i=1,pts
                eps = 2d0*((i-0.5d0)/real(pts)-0.5d0)
                result_step = gamma_mn(m,n,eps)&
                        * FermiFunction(eps,mu,beta)&
                        / ((1d0-eps**2d0)**2d0)
                result_tmp =result_tmp+result_step
              enddo
              int_gamma_mn = result_tmp/real(pts)
              return
      end function int_gamma_mn

