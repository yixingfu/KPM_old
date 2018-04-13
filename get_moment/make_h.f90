! Created=Tue 12 Dec 2017 03:28:22 PM STD
! Last Modified=Wed 04 Apr 2018 08:26:21 PM EDT
      !This file creates H
      !The matrix is stored as CSR(A,col,rp)
      !

        !call sleep(my_id*3)
      allocate(A(NNZ),col(NNZ),rp(N+1))
        if (SaveAll) then
        open (92,file=trim(outputfile_final)//".allout",&
                status="replace",access="stream",action="write")
        endif
        if (Inherit) then
                write(*,*)"inherit files"
                open(94,file=trim(outputfile_final)//".all",&
                        status="old",access="stream",action="read")
                read(94) A
                read(94) col
                read(94) rp
                read(94) norm_a
                read(94) norm_b
        
        else 
        if (fixedTwist) then
                Twist = OrigTwist*pi/real(L)
!                write(*,*) 'twist:',Twist
        else



        allocate(TwistAll(num_procs*3*seq_rep))
                call random_number(TwistAll)
!                write(*,*) 'proc',my_id
                Twist=TwistAll(my_id*3+1:my_id*3+3)
                write(*,*) 'twist:',Twist
        deallocate(TwistAll)
        
        

      Twist = (Twist)*pi/L! 0 to pi??
        endif

        expTwist = cdexp(III*Twist)

      txf = xf*cdexp(dcmplx(0d0,Twist(1)))
      txb = xb*cdexp(dcmplx(0d0,-Twist(1)))
      tyf = yf*cdexp(dcmplx(0d0,Twist(2)))
      tyb = yb*cdexp(dcmplx(0d0,-Twist(2)))
      tzf = zf*cdexp(dcmplx(0d0,Twist(3)))
      tzb = zb*cdexp(dcmplx(0d0,-Twist(3)))
!        write(*,*)'xf',xf(0)
!        write(*,*)'twist',Twist(1)
!        write(*,*)'exptwist',expTwist(1)
!        write(*,*)'mult',xf(0)*expTwist(1)
!        write(*,*)'txf',txf

      if (D.eq.2) then
              write(*,*) "2D"
              phase = 0
              if (QP) then
                      Wrnd = 0d0
                      WQP = W
                if (RandPhase) then
                      call random_number(phase)
                else 
                      phase = inputPhase
                endif

                      phase = phase*2d0*pi
                      P = 2.0d0*pi*fibonacci(fiboN-2)/fibonacci(fiboN)
                      Q = P
              else 
                      WQP = 0d0
                      Wrnd = W
              End if
                        TQP = 0d0
                        Trnd = 0d0
                        t0 = 1d0
                if (RandType.eq.RandHOP) then
                        WQP = 0d0
                        Wrnd = 0d0
                        t0 = dsqrt(1d0-W**2)
                        if (QP) then
                        TQP = W
                        Trnd = 0d0
                        else 
                        TQP = 0d0
                        Trnd = W
                        endif
                endif


              rp = 0
              col = 0
              col_ind = 1
              rp_ind = 1
               allocate(eps(L*L)) 
                eps_ind = 1
        if (RandType.eq.RandPOT) then
                do j=1,L
                do i=1,L

              eps(eps_ind) = WQP*quasiperiodic(i+0d0,j+0d0,P,Q,phase)&
                        + Wrnd*random2D(i+0d0,j+0d0,P,Q)
                eps_ind=eps_ind+1
                enddo
                enddo
                eps = eps - sum(eps)/real(L*L)
        else if (RandType.eq.RandHOP) then
                eps = 0d0
        endif

        eps_ind = 1
              do j=1,L! y
              do i=1,L! x
              do s = 0,1
                s_ = 1-s
                ! set row pointer and disorder term
                rp(rp_ind) = 5*rp_ind-4
                col(col_ind) = rp_ind
                rp_ind = rp_ind+1
                A(col_ind) = eps(eps_ind)
                col_ind = col_ind+1
                
                ! x forward
                t_tmp = t0 + Trnd*random2D(i-0.5d0,j+0d0,P,Q)&
                        + TQP*quasiperiodic(i-0.5d0,j+0d0,P,Q,phase)
                ind_r = xys2i(modulo(i-2,L)+1,j,s_,L)
                col(col_ind) = ind_r
                A(col_ind) = txf(s)*t_tmp
                col_ind = col_ind+1

                ! x backward
                t_tmp = t0 + Trnd*random2D(i+0.5d0,j+0d0,P,Q)&
                        + TQP*quasiperiodic(i+0.5d0,j+0d0,P,Q,phase)
                ind_r = xys2i(modulo(i,L)+1,j,s_,L)
                col(col_ind) = ind_r
                A(col_ind) = txb(s)*t_tmp
                col_ind = col_ind+1

                ! y forward
                t_tmp = t0 + Trnd*random2D(i+0.0d0,j-0.5d0,P,Q)&
                        + TQP*quasiperiodic(i+0.0d0,j-0.5d0,P,Q,phase)
                ind_r = xys2i(i,modulo(j-2,L)+1,s_,L)
                col(col_ind) = ind_r
                A(col_ind) = tyf(s)*t_tmp
                col_ind = col_ind+1

                ! y backward
                t_tmp = t0 + Trnd*random2D(i+0.0d0,j+0.5d0,P,Q)&
                        + TQP*quasiperiodic(i+0.0d0,j+0.5d0,P,Q,phase)
                ind_r = xys2i(i,modulo(j,L)+1,s_,L)
                col(col_ind) = ind_r
                A(col_ind) = tyb(s)*t_tmp
                col_ind = col_ind+1

              End do
                eps_ind = eps_ind+1
              End do
              End do
        rp(rp_ind) = NNZ+1

        deallocate(eps)
      else if (D.eq.3) then
              ! not differentiating QP or not
              rp = 0
              col = 0
              col_ind = 1
              rp_ind = 1
                allocate(eps(L*L*L))
                eps = 0
!                call ResetRandSeed(my_id*7)
!                call random_number(eps)
!                eps = W*2.0d0*(eps-0.5d0)
!                idum=-my_id*17
                call random3D(eps,W,my_id)
!                write(*,*)eps
                eps_ind = 1
              do k=1,L!z
              do j=1,L!y 
              do i=1,L!x
              do s = 0,1
                s_ = 1-s
                ! set row pointer and disorder term
                rp(rp_ind) = 7*rp_ind-6
                col(col_ind) = rp_ind
                rp_ind = rp_ind+1
                A(col_ind) = eps(eps_ind)
                col_ind = col_ind+1
!        write(*,*) rp_ind-1,',',col(col_ind-1),&
!        ',',real(A(col_ind-1)),',',imag(A(col_ind-1))
                
                ! x forward
                ind_r = xyzs2i(modulo(i-2,L)+1,j,k,s_,L)
                col(col_ind) = ind_r
                A(col_ind) = txf(s)
                col_ind = col_ind+1
!        write(*,*) rp_ind-1,',',col(col_ind-1),&
!        ',',real(A(col_ind-1)),',',imag(A(col_ind-1))

                ! x backward
                ind_r = xyzs2i(modulo(i,L)+1,j,k,s_,L)
                col(col_ind) = ind_r
                A(col_ind) = txb(s)
                col_ind = col_ind+1
!        write(*,*) rp_ind-1,',',col(col_ind-1),&
!        ',',real(A(col_ind-1)),',',imag(A(col_ind-1))

                ! y forward
                ind_r = xyzs2i(i,modulo(j-2,L)+1,k,s_,L)
                col(col_ind) = ind_r
                A(col_ind) = tyf(s)
                col_ind = col_ind+1
!        write(*,*) rp_ind-1,',',col(col_ind-1),&
!        ',',real(A(col_ind-1)),',',imag(A(col_ind-1))

                ! y backward
                ind_r = xyzs2i(i,modulo(j,L)+1,k,s_,L)
                col(col_ind) = ind_r
                A(col_ind) = tyb(s)
                col_ind = col_ind+1
!        write(*,*) rp_ind-1,',',col(col_ind-1),&
!        ',',real(A(col_ind-1)),',',imag(A(col_ind-1))

                ! caution! z has s->s, not s->s_
                ! z forward
                ind_r = xyzs2i(i,j,modulo(k-2,L)+1,s,L)
                col(col_ind) = ind_r
                A(col_ind) = tzf(s)
                col_ind = col_ind+1
!        write(*,*) rp_ind-1,',',col(col_ind-1),&
!        ',',real(A(col_ind-1)),',',imag(A(col_ind-1))

                ! z backward
                ind_r = xyzs2i(i,j,modulo(k,L)+1,s,L)
                col(col_ind) = ind_r
                A(col_ind) = tzb(s)
                col_ind = col_ind+1
!        write(*,*) rp_ind-1,',',col(col_ind-1),&
!        ',',real(A(col_ind-1)),',',imag(A(col_ind-1))

              End do
                eps_ind = eps_ind + 1
              End do
              End do
              End do
        rp(rp_ind) = NNZ+1
        deallocate(eps)

      endif


! -------------------------EIGENVALUE
        if (ExactSpectrum .or. ExactStates) then
                allocate(H_dense(N,N))
                call Sparse2Dense(N,NNZ,A,rp,col,H_dense)

        if (ExactStates) then
                JOBZ = 'V'
        else
                JOBZ = 'N'
        endif
              !Query
                allocate(work(2))
                allocate(rwork(3*N-2))
                allocate(EigVal(N))
                call zheev(JOBZ,'U',int(N),H_dense,int(N),EigVal,work,&
                        -1,rwork,info)
        if (info.ne.0) then
        write(*,*)"something wrong, exact diag query"
        endif
                lwork = work(1)
                write(*,*) "query result: ",lwork
              !End of query
                deallocate(work)
                allocate(work(lwork))
                call zheev(JOBZ,'U',int(N),H_dense,int(N),EigVal,work,&
                        lwork,rwork,info)
        if (info.ne.0) then
        write(*,*)"something wrong, exact diag "
        endif
        ! compare with new Lanzcos
!        call LanczosLowest(N,NNZ,A,rp,col,&
!                      1000*EIGVALCOUNT,EIGVALCOUNT,EigValLanc)
        ! empirical: 1.5*eigvalcount should be enough
        ! come back later. Not resolving the lowest ones
                deallocate(work,rwork)
 
        if (ExactStates) then
                write(*,*) "saving exact states NOT ready yet"
!                open(62,file=trim(outputfile_final)//".eigvec",&
!                        status="replace",access="stream",action="write")
!                write(62) N
!                write(62) H_dense
!                close(62)
        endif
        if (ExactSpectrum) then
                STARTPOINT=(N-EIGVALCOUNT)/2
                ENDPOINT = STARTPOINT+EIGVALCOUNT-1
                EigValTot = EigValTot + EigVal(STARTPOINT:ENDPOINT)
                EigValLancTot = EigValLancTot + EigValLanc
        do i_tmp=1,EIGVALCOUNT
                write(*,*) EigValTot(i_tmp),'vs',EigValLancTot(i_tmp),&
                        '@',my_id
        enddo
!                open(62,file=trim(outputfile_final)//".eigval",&
!                        status="replace",access="stream",action="write")
!                write(62) N
!                write(62) EigVal
!                close(62)
        endif
        deallocate(EigVal)

                deallocate(H_dense)
        endif        
! --------------------EIGEnVALUE end
      ! normalization same for both
        
      call LanczosBound(N,NNZ,A,rp,col,1000,Emax,Emin)
      norm_a = (Emax-Emin)/(2d0-0.2d0)
      norm_b = (Emax+Emin)/2
        call rescale_to_1(N,NNZ,A,rp,col,norm_a,norm_b)

        
        endif
        if (SaveAll) then 
                write(92) A
                write(92) col
                write(92) rp
                write(92) norm_a
                write(92) norm_b
        endif
