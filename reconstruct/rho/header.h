!Created=Wed 13 Dec 2017 01:48:29 PM STD
!Last Modified=Sun 25 Mar 2018 01:08:16 PM EDT
	character(200)::inputfile 
	character(200)::outputfile
	character(200)::arg_tmp
	logical::useFFT
	real*8::outputEmax
	integer::RLZmax,RLZmin
	integer::Nc,ForceNc,SetNtilde
	integer::Der
	real*8::norm_a,norm_b,tempTsum
	real*8,dimension(:),allocatable::mu_avg,mu2_avg,orig_mu
	integer::iT

	
	
	! rho calculation
	integer::Ntilde
	real*8::Ei,rhoi
	real*8,dimension(:),allocatable::Egrid,rho,rho_tot,Egrid_t

	
	! Jackson Kernel
	real*8,dimension(:),allocatable::gJ
	real*8::a_



	! parameters
	real*8,parameter::pi=3.1415926535897932384626433832795d0

	! index
	integer*8::i8,j8,k8
	integer::i,j,k
