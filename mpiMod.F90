module mpiMod
    use kinds   
    implicit none
    private 
    save

    include 'mpif.h'

    INTEGER(kind = i4) :: i_err, &
                        taskid, &
                        numtasks, &
                        numworkers
   
    INTEGER(kind = i4), PUBLIC :: MASTER 

    integer status(MPI_STATUS_SIZE)

    PUBLIC :: startMPI, &
              mpiAbort, &
              thisProc, &
              totalWorkers

    contains

    subroutine startMPI()
        call MPI_INIT( i_err )
        call MPI_COMM_RANK( MPI_COMM_WORLD, taskid, i_err )
        call MPI_COMM_SIZE( MPI_COMM_WORLD, numtasks, i_err )
        
        MASTER = 0
        numworkers = numtasks-1
    end subroutine

    subroutine mpiAbort()
        call MPI_ABORT(MPI_COMM_WORLD, i_err)
    end subroutine

    function thisProc() result(procNo)
        INTEGER(kind = i4 ) :: procNo
        procNo = taskid
    end function

    function totalWorkers() result(nWorkers)
        INTEGER(kind = i4 ) :: nWorkers
        nWorkers = numworkers
    end function

end module mpiMod
