program filterQGCM
    !use fields
    use configMod
    !use mpiMod
    !use ioMod
    use gridMod
    implicit none

    !! call startMPI()   ! initialize MPI

    call makeConfig() !!! Reads the i/o location and variables to read 
    call init_grid()
end program filterQGCM
