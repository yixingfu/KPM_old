! peeled from make_h.f90=Thu 10 May 2018 02:58:02 PM DST
! Last Modified=Thu 10 May 2018 02:58:24 PM DST
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
