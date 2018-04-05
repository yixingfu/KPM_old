! Created=Wed 13 Dec 2017 01:47:49 PM STD
! Last Modified=Thu 05 Apr 2018 02:38:34 PM EDT

      program main
              implicit none
              include "header.h"


              call getarg(1,inputfile)
              call getarg(2,arg_tmp)
              read(arg_tmp,*) useFFT
              call getarg(3,arg_tmp)
              read(arg_tmp,*) outputEmax
              call getarg(4,arg_tmp)
                read(arg_tmp,*) RLZmin
              call getarg(5,arg_tmp)
              read(arg_tmp,*) RLZmax
                call getarg(6, arg_tmp)
                read(arg_tmp,*) Der
                call getarg(7,arg_tmp)
        read(arg_tmp,*) ForceNc
                call getarg(8,arg_tmp)
        read(arg_tmp,*) SetNtilde

              outputfile = inputfile
        write(*,*)RLZmin,RLZmax,Der

                do k = RLZmin,RLZmax
                write(*,*)k

              !first read, and prepare Egrid
              include "read_input.f90"

              !transform
              include "get_rho.f90"


                enddo
!              do i=1,Ntilde
!                      write(*,*)Egrid(i),rho_tot(i)
!              End do
              if (Der.eq.1) then
                outputfile=trim(outputfile)//'d'
                endif
             write(outputfile,'(a,i6.6)')trim(outputfile),ForceNc
              open(13,file=trim(outputfile)//".dat",status="replace",&
                      form="unformatted",access="stream")
              write(13)RLZmax+1-RLZmin,Ntilde,Egrid,rho_tot,rho2_tot
              close(13)
                deallocate(Egrid,rho_tot,rho2_tot)
              contains 
                      include "Chebyshev.f90"
      End program main



