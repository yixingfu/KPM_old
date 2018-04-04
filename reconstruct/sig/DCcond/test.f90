        program main
                implicit none
                real*8::eps
                integer::i,m,n
                integer,parameter::Nc=256
                complex*16,dimension(:,:),allocatable::intgammamn
                allocate(intgammamn(0:Nc-1,0:Nc-1))
                intgammamn = 0
        read(*,*)m,n
          intgammamn(m,n)=int_gamma_mn(m,n,0d0,10000000d0,400)
        write(*,*)"Gamma(",m,",",n,")=",intgammamn(m,n)

                contains
                include "GammaMN.f90"
                include "Chebyshev.f90"
                include "StatMech.f90"


        end program main
