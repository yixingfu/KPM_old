! Created=Wed 13 Dec 2017 02:27:33 PM STD
! Last Modified=Tue 06 Feb 2018 08:44:44 PM EST
      ! This file makes kernel gn for n=0,Nc-1

      ! Jackson kernel
      ! gn = [(N-n+1)cos(pi n/(N+1)sin(pi n/(N+1))cot(pi/(N+1))]/(N+1)
      allocate(gJ(0:Nc-1),hm(0:Nc-1))
        a_ = pi/(Nc+1)
      do i = 0,Nc-1
              gJ(i) = ((Nc-i+1)*dcos(a_*i)+dsin(a_*i)/dtan(a_)) &
                      /(Nc+1)
      End do
        hm=2d0
        hm(0)=1d0
        write(*,*) "kernel prepared"
