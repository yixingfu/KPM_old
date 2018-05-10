! created=Tue 12 Dec 2017 03:55:22 PM STD
! Last Modified=Thu 10 May 2018 03:05:34 PM DST
      ! changing indices between x,y,z,s to single index
      ! s=0 for spin up, s=1 for spin down
      integer*8 function xyzs2i(x,y,z,s,L)
          integer,intent(in)::x,y,z
          integer,intent(in)::L,s
          xyzs2i = 2_8*((z-1)*L**2_8+(y-1)*L+(x-1))+1+s
          return
      end function xyzs2i
      ! changing indices between x,y,s to single index
      ! s=0 for spin up, s=1 for spin down
      integer*8 function xys2i(x,y,s,L)
          integer,intent(in)::x,y
          integer,intent(in)::L,s
          xys2i = 2_8*((y-1)*L+(x-1))+1+s
          return
      end function xys2i
