module mpiMod
    use kinds   
    implicit none

    include 'mpif.h'

    integer :: i_err, taskid, numtasks, numworkers, msg_tag,  &
        & MASTER, source, dest, FROM_MASTER, FROM_WORKER, ERRORCODE, &
        & startCol_tag, endCol_tag, dataSize_tag, OL_UVEL_tag, OL_VVEL_tag, &
        & OL_TAUX_tag, OL_TAUY_tag, OL_PPA_tag

    integer status(MPI_STATUS_SIZE)

    contains

    subroutine startMPI()
        call MPI_INIT( i_err )
        call MPI_COMM_RANK( MPI_COMM_WORLD, taskid, i_err )
        call MPI_COMM_SIZE( MPI_COMM_WORLD, numtasks, i_err )
        !print *, 'FROM INITIALIZATIIN FUNCTION numtasks =', numtasks

        MASTER = 0
        msg_tag = 0
        numworkers = numtasks-1
        FROM_MASTER =1
        FROM_WORKER =2

        startCol_tag = 1000
        endCol_tag = 1100
        dataSize_tag = 2000
        OL_UVEL_tag = 10000
        OL_VVEL_tag = 20000
        OL_TAUX_tag = 30000
        OL_TAUY_tag = 40000
        OL_PPA_tag = 50000

    end subroutine

    subroutine mpiAbort()
        call MPI_ABORT(MPI_COMM_WORLD, i_err)
    end subroutine

end module mpiMod
