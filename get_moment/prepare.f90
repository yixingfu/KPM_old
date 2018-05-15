! Created =Tue 12 Dec 2017 03:10:19 PM STD
! Last Modified=Thu 10 May 2018 03:05:35 PM DST
! 

      ! This file prepares a few derived parameters from input file
      N = 2*(L**D)
        if (BHZ) then
        NNZ = (1+4*D)*N
        else
      NNZ = (1+2*D)*N ! fwd & bwd each site per dim + disorder
        endif

      JNNZ = (2)*N ! fwd & bwd each site @ x
!      write(*,*)D,N,L,NNZ


      EigValTot = 0
      EigValLancTot = 0
