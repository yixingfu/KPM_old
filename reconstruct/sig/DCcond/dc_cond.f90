! Created=Tue 16 Jan 2018 05:13:10 AM STD
! Last Modified=Wed 24 Jan 2018 05:03:12 PM STD
      ! Computes DC conductivity
      ! Data to be prepared: mu_mn(0:Ntilde-1,0:Ntilde-1)
      ! Need to reserve a few loop variables: m,n
      ! prepare phase
        allocate(intgammamn(0:Ntilde-1,0:Ntilde-1))
!      do m=0,Ntilde-1
!      write(*,*)m
!      do n=0,Ntilde-1
!        intgammamn(m,n) = int_gamma_mn(m,n,0d0,10000000d0,100000)
!      End do
!      End do 
      open(18,file="intgammamnT1000.dat",&
                form="unformatted",access="stream")
      read(18)intgammamn
      close(18)

      ! working phase
      dc_conductivity = 0d0
      do m=0,Ntilde-1
      do n=0,Ntilde-1
        dc_conductivity = dc_conductivity+mu2d(m,n)*intgammamn(m,n)
      End do
      End do 
      deallocate(intgammamn)
