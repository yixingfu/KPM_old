! Created =Tue 12 Dec 2017 03:10:19 PM STD
! Last Modified=Sun 13 May 2018 05:26:04 PM DST
! 

      ! This file prepares a few derived parameters from input file
      N = 2*(L**D)
      NNZ = (1+2*D)*N ! fwd & bwd each site per dim + disorder
!        if (BHZ) then --- for BHZ 
            NNZ = (1+4*D)*N ! cross spin term and same spin term
!        endif

      JNNZ = (2)*N ! fwd & bwd each site @ x
!      write(*,*)D,N,L,NNZ


      EigValTot = 0
      EigValLancTot = 0
