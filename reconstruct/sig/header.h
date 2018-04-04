!Created=Wed 13 Dec 2017 01:48:29 PM STD
!Last Modified=Tue 13 Feb 2018 11:29:37 PM EST
	character(200)::inputfile 
	character(200)::outputfile
	character(200)::arg_tmp
	logical::useFFT
	real*8::EmaxNONFFT
	real*8::outputEmax
	integer::RLZmax,RLZmin
	integer::Nc
	real*8::norm_a,norm_b
	real*8,dimension(:,:),allocatable::mu_avg,mu2_avg
	
	! sig calculation
	integer::Ntilde,Noutput
	real*8::Ei,rhoi
	real*8,dimension(:),allocatable::xygrid,xygrid_ref
	real*8,dimension(:,:),allocatable::Jxy,Jxy_tot,Jxy_ref
	real*8::Jij,xi,yj
	integer,parameter::ref_lev=4

	
	! Jackson Kernel
	real*8,dimension(:),allocatable::gJ,hm
	real*8::a_



	! parameters
	real*8,parameter::pi=3.1415926535897932384626433832795d0

	! index & stats
	integer*8::i8,j8,k8
	integer::i,j,k,m,n
	integer::badfiles
	integer::ierr

	! useFFT version
	real*8,dimension(:,:),allocatable::mu_tilde
	real*8,dimension(:),allocatable::Jtemp_in,Jtemp_out
	
	! DC
	real*8,dimension(:,:),allocatable::intgammamn
	real*8::DCtot, dc_conductivity
	
	
