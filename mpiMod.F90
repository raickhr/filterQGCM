module mpiMod
include 'mpif.h'
    integer :: i_err, taskid, numtasks, numworkers, mtype,  &
        & MASTER, source, dest, FROM_MASTER, FROM_WORKER, ERRORCODE

    integer status(MPI_STATUS_SIZE)

    contains

    subroutine startMPI()
        call MPI_INIT( i_err )
        call MPI_COMM_RANK( MPI_COMM_WORLD, taskid, i_err )
        call MPI_COMM_SIZE( MPI_COMM_WORLD, numtasks, i_err )
        !print *, 'FROM INITIALIZATIIN FUNCTION numtasks =', numtasks

        MASTER = 0
        mtype = 0
        numworkers = numtasks-1
        FROM_MASTER =1
        FROM_WORKER =2

    end subroutine

end module mpiMod
