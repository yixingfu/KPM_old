! Created=Wed 13 Dec 2017 01:47:49 PM STD
! Last Modified=Tue 27 Mar 2018 04:59:03 PM EDT

      program main

                use mkl_dfti
                use lapack95
                use f95_precision
              implicit none
                include "mkl.fi"
                include "fftw3.f"
                include "fftw3_mkl.f"
              include "header.h"



              call getarg(1,inputfile)
              call getarg(2,arg_tmp)
              read(arg_tmp,*) useFFT
              call getarg(3,arg_tmp)
              read(arg_tmp,*) RLZmin
              call getarg(4,arg_tmp)
              read(arg_tmp,*) RLZmax
        write(*,*)RLZmax
        badfiles=0
              outputfile = inputfile
              call getarg(5,arg_tmp)
              read(arg_tmp,*) EmaxNONFFT
              call getarg(6,arg_tmp) !optional arg
        if (arg_tmp.eq."") then
                Noutput = 0
        else
              read(arg_tmp,*) Noutput
        endif

                do k = RLZmin,RLZmax
                write(*,*)k

              !first read, and prepare Egrid
              include "read_input.f90"

                write(*,*)k, "read done"
              !transform
              include "get_Jxy.f90"
                write(*,*)k, "get Jxy done"
!              include "get_sig.f90"

                write(*,*)'norm_b:',norm_b
                enddo
!rescale xygrid now
        if (useFFT) then        
        xygrid = xygrid*norm_a+norm_b
        endif

              open(13,file=trim(outputfile)//".dat",status="replace",&
                      form="unformatted",access="stream")
              write(13)RLZmax+1-badfiles-RLZmin,&
                        Noutput,xygrid,Jxy_tot,DCtot
              close(13)
              open(22,file="output.log",status="replace")
                write(22,*)"RLZ #, Nmax,xygrid,Jxy"
                close(22)
                deallocate(xygrid,Jxy_tot)
              contains 
                      include "Chebyshev.f90"
                      include "CosTrans.f90"
                      include "GammaMN.f90"
                      include "StatMech.f90"
                      include "interp.f90"
      End program main



