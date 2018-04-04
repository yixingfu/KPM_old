!File=StatMech.f95
!Author=yxfu
!Created=Sat 30 Sep 2017 11:03:24 PM DST
!Last Modified=Tue 16 Jan 2018 05:06:42 AM STD
        
        real*8 function FermiFunction(E,mu,beta)
                real*8,intent(in)::E,mu,beta
                if (beta==0d0) then 
                        if (E.lt.mu) then
                                FermiFunction=1d0
                        else
                                FermiFunction=0d0
                        endif
                else
                        FermiFunction = 1d0/(exp((E-mu)*beta)+1d0)
                endif
                return
        End function FermiFunction
        real*8 function BoseFunction(E,mu,beta)
                real*8,intent(in)::E,mu,beta
                
                BoseFunction = 1d0/(exp((E-mu)*beta)-1d0)
                return
        End function BoseFunction

