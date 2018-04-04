!Created=Mon 12 Feb 2018 08:34:09 PM STD
!Last Modified=Tue 13 Feb 2018 11:05:16 AM EST
      ! This subroutine does 1d interpolation
      ! linear interpolation only for now

      subroutine interp(grid_in, arr_in, grid_out, arr_out)
              !grid_in: the x axis of the data
              !arr_in: data.
              !grid_out: target grid
              !arr_out: result. Should be ready allocated

              real*8,dimension(:),intent(in)::grid_in,arr_in,grid_out
              real*8,dimension(:),intent(inout)::arr_out
              integer::N_in,N_out,i,j
              integer,dimension(1)::temp_ind_l,temp_ind_r
              real*8::xi,yi,xl,xr,yl,yr

!                if (size(grid_in) .ne. size(arr_in)) then
!                        write(*,*) "Size not matching: interp, in"
!                endif
!                if (size(grid_out) .ne. size(arr_out)) then
!                        write(*,*) "Size not matching: interp, out"
!                endif
              N_in=size(grid_in)
              N_out=size(grid_out)
!              write(*,*)grid_in(1),grid_in(N_in)
!              write(*,*)grid_out(1),grid_out(N_out)
!              write(*,*)arr_in(1),arr_in(N_in)
!              write(*,*)arr_out(1),arr_out(N_out)

              do i=1,N_out
                      xi=grid_out(i)
                      if (xi .gt. maxval(grid_in)) then
                              yi=0
                      else if (xi .le. minval(grid_in)) then
                              yi=0
                      else
                        temp_ind_l=minloc(abs(xi-grid_in),&
                          (xi .gt. grid_in))
                        temp_ind_r=minloc(abs(xi-grid_in),&
                          (xi .le. grid_in))
                        xl=grid_in(temp_ind_l(1))
                        yl=arr_in(temp_ind_l(1))
                        xr=grid_in(temp_ind_r(1))
                        yr=arr_in(temp_ind_r(1))
                        yi=((xr-xi)*yl+(xi-xl)*yr)/(xr-xl)
                      endif
                      arr_out(i)=yi
              End do
              return
      end subroutine interp


      subroutine interp2d(gridx_in,gridy_in,arr_in,gridx_out,gridy_out&
                                        ,arr_out)
              !grid_in: the x axis of the data
              !arr_in: data.
              !grid_out: target grid
              !arr_out: result. Should be ready allocated

              real*8,dimension(:,:),intent(in)::arr_in
              real*8,dimension(:,:),intent(inout)::arr_out
              real*8,dimension(:),intent(in)::gridx_in,gridy_in
              real*8,dimension(:),intent(in)::gridx_out,gridy_out
              real*8,dimension(:),allocatable::arr_tmp
              real*8,dimension(:,:),allocatable::arr2d_tmp
              integer::i,j,Nx_in,Ny_in,Nx_out,Ny_out
              Nx_in=size(gridx_in)
              Ny_in=size(gridy_in)
              Nx_out=size(gridx_out)
              Ny_out=size(gridy_out)
              allocate(arr2d_tmp(Nx_in,Ny_out))
              allocate(arr_tmp(Ny_out))
              do i=1,Nx_in
              arr_tmp=0d0
              call interp(gridy_in,arr_in(i,:),gridy_out,arr_tmp)
              
              arr2d_tmp(i,:)=arr_tmp(:)
              enddo
              deallocate(arr_tmp)
                
              allocate(arr_tmp(Nx_out))
              do i=1,Ny_out
                arr_tmp=0d0
           call interp(gridx_in,arr2d_tmp(:,i),gridx_out,arr_tmp)
                arr_out(:,i)=arr_tmp(:)
              enddo
              deallocate(arr_tmp,arr2d_tmp)
              return 

      End subroutine interp2d

