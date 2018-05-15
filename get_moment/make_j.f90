! Created=Wed 13 Dec 2017 11:51:36 PM STD
! Last Modified=Thu 10 May 2018 03:02:41 PM DST
      ! This file makes J operator

      if (D.eq.2) then
          write(*,*) "J not implemented for 2D"
      else if (D.eq.3) then
          ! not differentiating QP or not
          write(*,*) "3D,J"
          Jrp = 0
          Jcol= 0
          col_ind = 1
          rp_ind = 1
          ! txf*(-i),txb*(i), that's it.
          Jtxf = -txf*III
          Jtxb = txb*III

          do k=1,L!z
          do j=1,L!y
          do i=1,L!x
          do s=0,1
          s_ = 1-s
          Jrp(rp_ind) = 2*rp_ind-1
          rp_ind = rp_ind+1

          ! x_forward
          ind_r = xyzs2i(modulo(i-2,L)+1,j,k,s_,L)
          Jcol(col_ind) = ind_r
          JA(col_ind) = Jtxf(s,s_)
          col_ind = col_ind+1
!        write(*,*) rp_ind-1,',',Jcol(col_ind-1),&
!        ',',real(JA(col_ind-1)),',',imag(JA(col_ind-1))

          ! x backward
          ind_r = xyzs2i(modulo(i,L)+1,j,k,s_,L)
          Jcol(col_ind) = ind_r
          JA(col_ind) = Jtxb(s,s_)
          col_ind = col_ind+1
!        write(*,*) rp_ind-1,',',Jcol(col_ind-1),&
!        ',',real(JA(col_ind-1)),',',imag(JA(col_ind-1))

          End do
          End do
          End do
          End do
          Jrp(rp_ind) = JNNZ+1
      endif

