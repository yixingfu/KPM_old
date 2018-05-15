! Created=Tue 12 Dec 2017 03:19:48 PM STD
! Last Modified=Thu 10 May 2018 03:05:35 PM DST
        ! read inputfile name from command line
        call getarg(1,inputfile)
        ! read inputs D,L,Nc,W,QP?, from input file
        open(11, file=trim(inputfile))
        read(11,*) D,L,Nc,W,Rep,REALIZATION0
        read(11,*) QP,fixedTwist, slowOPTCOND
        read(11,*) task, RandType
        read(11,*) OrigTwist
        read(11,*) outputfile
        read(11,*) Inherit, SaveAll
        read(11,*) ExactSpectrum, ExactStates
        read(11,*) RandPhase, inputPhase
        read(11,*) fiboM,commC,setQ
        read(11,*) BHZ,BHZ_SPIN,BHZ_M
        close(11)

        if (QP) then
            fiboN = L
            L = fibonacci(fiboN)
        endif

!         check
!         write(*,*) D,L,Nc,W,QP,outputfile



