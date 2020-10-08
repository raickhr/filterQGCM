module gridMod
    ! This module has the information about the grid sizes
    
    use kinds
    use configMod
    use mpiMod
    
    implicit none
    private 
    save


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

    REAL(kind=r8):: prevTimeVal, timeVal
    CHARACTER(len = char_len_short), PUBLIC :: timeUnits 

    PUBLIC :: init_grid, &
              setOcnPgridXYsizeto, &
              setOcnTgridXYsizeto, &
              setAtmPgridXYsizeto, &
              setAtmTgridXYsizeto, &
              saveReadXpoYpo, &
              getDxpo, getDypo, &
              getPrevTimeVal, getCurrTimeVal, &
              setPrevTimeVal, setCurrTimeVal


    contains
    subroutine init_grid()
        !! This subroutine initializes grid sizes from the config file given at terminal

        INTEGER(kind = i4) :: err

        if (thisProc() == MASTER) then 
            print *,''
            print *, 'Initializing grid ...' 

            call setAtmXYto(nxta, nyta,  nza)
            call setOcnXYto(nxto, nyto,  nzo)

        endif

        prevTimeVal = -111.00
        timeVal = -111.00

        call MPI_BCAST(nxto, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, err)
        call MPI_BCAST(nyto, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, err)
        call MPI_BCAST(nzo, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, err)
        call MPI_BCAST(nxta, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, err)
        call MPI_BCAST(nyta, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, err)
        call MPI_BCAST(nza, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, err)
        
        call MPI_Barrier(MPI_COMM_WORLD, err)

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

    subroutine setOcnPgridXYsizeto(ocnX, ocnY )
        integer(kind = i4), INTENT(OUT) :: ocnX, ocnY
        ocnX = nxpo
        ocnY = nypo
    end subroutine 

    subroutine setOcnTgridXYsizeto(ocnX, ocnY )
        integer(kind = i4), INTENT(OUT) :: ocnX, ocnY
        ocnX = nxto
        ocnY = nyto
    end subroutine 

    subroutine setAtmPgridXYsizeto(AtmX, AtmY )
        integer(kind = i4), INTENT(OUT) :: AtmX, AtmY
        AtmX = nxpa
        AtmY = nypa
    end subroutine 

    subroutine setAtmTgridXYsizeto(AtmX, AtmY )
        integer(kind = i4), INTENT(OUT) :: AtmX, AtmY
        AtmX = nxta
        AtmY = nyta
    end subroutine 

    subroutine saveReadXpoYpo(inxpo, inypo)

        REAL(kind = r8), INTENT(IN) :: inxpo(nxpo), inypo(nypo)
        xpo(:) = inxpo(:)
        ypo(:) = inypo(:)

    end subroutine

    function getDxpo() result(returnVal)
        REAL(kind = r8) :: returnVal
        returnVal = xpo(2) - xpo(1)
    end function

    function getDypo() result(returnVal)
        REAL(kind = r8) :: returnVal
        returnVal = ypo(2) - ypo(1)
    end function

    function getPrevTimeVal() result(returnVal)
        REAL(kind = r8) :: returnVal
        returnVal = prevTimeVal
    end function

    function getCurrTimeVal() result(returnVal)
        REAL(kind = r8) :: returnVal
        returnVal = timeVal
    end function

    subroutine setPrevTimeVal( inTimeVal) 
        REAL(kind = r8), INTENT(IN) :: inTimeVal
        INTEGER(kind = i4) :: err
        prevTimeVal = inTimeVal

        call MPI_BCAST(prevTimeVal, 1, MPI_REAL , MASTER, MPI_COMM_WORLD, err)
        call MPI_Barrier(MPI_COMM_WORLD, errorCode)
    end subroutine

    subroutine setCurrTimeVal( inTimeVal) 
        REAL(kind = r8), INTENT(IN) :: inTimeVal
        INTEGER(kind = i4) :: err
        timeVal = inTimeVal

        call MPI_BCAST(prevTimeVal, 1, MPI_REAL , MASTER, MPI_COMM_WORLD, err)
        call MPI_Barrier(MPI_COMM_WORLD, errorCode)
    end subroutine



    

end module gridMod