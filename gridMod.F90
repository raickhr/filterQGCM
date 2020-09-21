module gridMod
    ! This module has the information about the grid sizes
    
    use kinds
    use configMod
    implicit none

    INTEGER(kind = i4) :: atm_npx, &     !! number of atm p points for x dim
                          atm_npy, &     !! number of atm p points for y dim
                          atm_ntx, &     !! number of atm t points for x dim
                          atm_nty, &     !! number of atm t points for y dim
                          atm_nl         !! number of atm layers 


    INTEGER(kind = i4) :: ocn_npx, &     !! number of ocn p points for x dim
                          ocn_npy, &     !! number of ocn p points for y dim
                          ocn_ntx, &     !! number of ocn t points for x dim
                          ocn_nty, &     !! number of ocn t points for y dim
                          ocn_nl         !! number of ocn layers 

    real(kind=r4), ALLOCATABLE, DIMENSION(:)  :: xpo,  &  !! x values for ocn p points
                                                 ypo,  &  !! y values for ocn p points
                                                 xto,  &  !! x values for ocn t points
                                                 yto,  &  !! y values for ocn t points 
                                                 zo,   &  !! ocn midlayer depths 
                                                 zio,  &  !! ocn interface depths
                                                 xpa,  &  !! x values for atm p points
                                                 ypa,  &  !! y values for atm p points
                                                 xta,  &  !! x values for atm t points
                                                 yta, &    !! y values for atm t points
                                                 za, &    !! atm midlayer heights  
                                                 zia      !! atm interaface heights

    INTEGER(kind=i4) :: nxpo, nypo, &   !! ocn grid size for p points
                        nxto, nyto, &   !! ocn grid size for t points (number of t points are one less than number of p points) 
                        nzo,         &   !! ocn number of midlayer depths
                        nzio,        &   !! ocn number of interfaces
                        nxpa, nypa, &   !! atm gird size for p points
                        nxta, nyta, &    !! atm grid size for t points
                        nza,         &   !! atm number of midlayer depths
                        nzia           !! atm number of interfaces

    contains
    subroutine init_grid()
        !! This subroutine initializes grid sizes from the config file given at terminal
        print *,''
        print *, 'Initializing grid ...' 
        nxto = nxto_c 
        nyto = nyto_c
        nzo = nzo_c
        nxta = nxta_c
        nyta = nyta_c
        nza = nza_c

        nxpo = nxto+1 
        nypo = nyto+1
        nzio = nzo-1
        nxpa = nxta+1
        nypa = nyta+1
        nzia = nza-1

        allocate(xpo(nxpo),  & 
                 ypo(nypo),  &
                 xto(nxto),  &
                 yto(nyto),  &
                 zo(nzo),    &  
                 zio(nzio),  &
                 xpa(nxpa),  &
                 ypa(nypa),  &
                 xta(nxta),  &
                 yta(nyta),  &
                 za(nza),    &  
                 zia(nzia))

    end subroutine


    

end module gridMod