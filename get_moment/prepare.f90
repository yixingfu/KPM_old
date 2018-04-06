! Created =Tue 12 Dec 2017 03:10:19 PM STD
! Last Modified=Wed 04 Apr 2018 03:44:24 PM EDT
! 

      ! This file prepares a few derived parameters from input file
      N = 2*(L**D)
      NNZ = (1+2*D)*N ! fwd & bwd each site per dim + disorder

      JNNZ = (2)*N ! fwd & bwd each site @ x
!      write(*,*)D,N,L,NNZ


        EigValTot = 0
        EigValLancTot = 0
