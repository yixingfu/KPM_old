! Created=Wed 13 Dec 2017 01:47:49 PM STD
! Last Modified=Thu 05 Apr 2018 02:38:34 PM EDT

      program main
              implicit none
              include "mpif.h"
              include "header.h"

        call MPI_INIT(ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD,my_id,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,num_procs,ierr)

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
        if (my_id.eq.0) then
           write(*,*)RLZmin,RLZmax,Der
        endif
                N_rlz_actual = 0
                do k_base = 0,((RLZmax-RLZmin+1)/num_procs)
                k = k_base*num_procs+RLZmin+my_id
                if (k.gt. RLZmax) then
                        cycle
                endif
                write(*,*)k,"@",my_id,"max",RLZmax

              !first read, and prepare Egrid
              include "read_input.f90"

              !transform
              include "get_rho.f90"


                enddo
!!!no longer using
!!!              if (Der.eq.1) then
!!!                outputfile=trim(outputfile)//'d'
!!!                endif
        call MPI_REDUCE(rho_tot,rho_tot_collect,Ntilde,&
                MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(rho2_tot,rho2_tot_collect,Ntilde,&
                MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(N_rlz_actual,N_rlz_collect,1,&
                MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)

     if (my_id .eq. 0) then
        write(*,*) N_rlz_collect
             write(outputfile,'(a,i6.6)')trim(outputfile),ForceNc
              open(13,file=trim(outputfile)//".dat",status="replace",&
                      form="unformatted",access="stream")
                write(13)N_rlz_collect,Ntilde,Egrid,&
                rho_tot_collect,rho2_tot_collect
!              write(13)RLZmax+1-RLZmin,Ntilde,Egrid,rho_tot,rho2_tot
              close(13)
        endif
                deallocate(Egrid,rho_tot,rho2_tot,rho2_tot_collect,&
                rho_tot_collect)
        call MPI_FINALIZE(ierr)
              contains 
                      include "Chebyshev.f90"
      End program main



