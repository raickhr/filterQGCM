module constants
    use kinds
    implicit none
    
    REAL(kind =r8) :: f0  !! coriolis paramater
                      

    REAL(kind=r8),PARAMETER :: rho_o = 1000, & !! ocean density
                               rho_a = 1  !! atmosphere density 

end module
