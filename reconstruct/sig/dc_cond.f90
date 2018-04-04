! Created=Tue 16 Jan 2018 05:13:10 AM STD
! Last Modified=Mon 26 Mar 2018 06:17:00 PM EDT
      ! Computes DC conductivity
      ! Data to be prepared: mu_mn(0:Ntilde-1,0:Ntilde-1)
      ! Need to reserve a few loop variables: m,n
      ! prepare phase
!      do m=0,Ntilde-1
!      write(*,*)m
!      do n=0,Ntilde-1
!        intgammamn(m,n) = int_gamma_mn(m,n,0d0,10000000d0,100000)
!      End do
!      End do 

      ! working phase
      dc_conductivity = 0d0
        !no, not in use.
!      do m=0,Nc-1
!      do n=0,Nc-1
!        dc_conductivity = dc_conductivity+mu_tilde(m,n)*intgammamn(m,n)
!      End do
!      End do 
