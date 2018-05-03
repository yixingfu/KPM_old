! Created=Wed 13 Dec 2017 02:27:33 PM STD
! Last Modified=Fri 02 Mar 2018 03:00:50 PM EST
      ! This file makes kernel gn for n=0,Nc-1

      ! Jackson kernel
      ! gn = [(N-n+1)cos(pi n/(N+1)sin(pi n/(N+1))cot(pi/(N+1))]/(N+1)
        if (ForceNc .eq. 0) then
                ForceNc=Nc
        endif
      allocate(gJ(0:Nc-1))
        a_ = pi/(ForceNc+1)
        gJ=0
      do i = 0,ForceNc-1
              gJ(i) = ((ForceNc-i+1)*dcos(a_*i)+dsin(a_*i)/dtan(a_)) &
                      /(ForceNc+1)
      End do

