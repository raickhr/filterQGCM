module filter
    use kinds
    use gridMod
    use fields
    use workDiv

    implicit none


    real(kind = r8), ALLOCATABLE, DIMENSION(:,:) :: kernel

    real(kind = r8) :: dArea


    contains

    subroutine makeKernel(filterLength, gridType, ocnORatmGrid )
        REAL(kind =r4), INTENT(IN) :: filterLength
        CHARACTER, INTENT(IN) :: gridType
        CHARACTER, OPTIONAL :: ocnORatmGrid
        REAL(kind = r8) :: dx, dy

        INTEGER(kind=i4) :: knx, kny , &
                            i, j, &  !! iterators
                            n  !! dummy number
        real(kind = r8), ALLOCATABLE, DIMENSION(:,:) :: kX2D, kY2D, &
                                                        distFromCenter, &
                                                        ell_by2, weight

        if (ALLOCATED(kernel)) then
            DEALLOCATE(kernel)
        endif

        if (PRESENT(ocnORatmGrid)) then
            if ((ocnORatmGrid .EQ. 'A') .OR. &
                (ocnORatmGrid .EQ. 'a')) then
                if ( (gridType .EQ. 'T') .OR. &
                   (gridType .EQ. 't')) then 
                    dx = xta(2) - xta(1) 
                    dy = yta(2) - yta(1)
                else
                    dx = xpa(2) - xpa(1) 
                    dy = ypa(2) - ypa(1)
                endif
            else
                if ( (gridType .EQ. 'T') .OR. &
                   (gridType .EQ. 't')) then 
                    dx = xto(2) - xto(1) 
                    dy = yto(2) - yto(1)
                else
                    dx = xpo(2) - xpo(1) 
                    dy = ypo(2) - ypo(1)
                endif
            endif
        else    
            print*,'default dx dy'
            dx = xpo(2) - xpo(1) 
            dy = ypo(2) - ypo(1)
        endif

        dArea = dx*dy

        !print *,'in make kernel dx, dy= ', dx, dy

        knx = NINT(((filterLength/2.0)/dx)+1) * 2 + 1
        kny = NINT(((filterLength/2.0)/dy)+1) * 2 + 1

        !print *,'in make kernel knx, kny= ', knx, kny
        

        allocate(distFromCenter(knx, kny), &
                 kernel(knx, kny), &
                 kX2D(knx, kny), kY2D(knx, kny), &
                 ell_by2(knx, kny), weight(knx, kny))

        n = -knx/2
        do i = 1, knx
            kX2D(i,:) = n * dx
            n = n + 1
        enddo

        n = -kny/2
        do j = 1, kny
            kY2D(:,j) = n * dy
            n = n + 1
        enddo

        distFromCenter(:,:) = SQRT(kX2D*kX2D + kY2D *kY2D)

        weight(:,:) = 0.5
        ell_by2(:,:) = filterLength/2

        weight = weight-0.5*tanh((distFromCenter-ell_by2)/10.0) 
        weight = weight * dArea
        
        kernel = weight/sum(weight)

        ! where (distFromCenter > filterLength/2)
        !     kernel = 0
        ! endwhere
        
        
    end subroutine


    subroutine getFilterFields(inField, outField)
        REAL(kind = r8), INTENT(IN) :: inField(nxpo,nypo)
        REAL(kind = r8), INTENT(OUT) :: outField(nxpo,nypo)

        INTEGER (kind=i4) :: strt_col, ed_col, data_size, row_size, &
                             i, j, &
                             knx, kny, &
                             is, ie, js, je

        REAL(kind = r8) , ALLOCATABLE, DIMENSION(:,:):: UN_FILT_FEILD, FILT_FEILD

        knx = size(kernel,1)
        kny = size(kernel,2)

        is = knx/2 +1
        ie = nxpo + knx/2+1
        js = kny/2 +1
        je = nypo + kny/2+1

        allocate(UN_FILT_FEILD(nxpo + knx-1, nypo+kny-1), &
                FILT_FEILD(nxpo + knx-1, nypo+kny-1))

        UN_FILT_FEILD(:,:) = 0
        FILT_FEILD(:,:) = 0

        UN_FILT_FEILD(is:ie, js:je) = inField(:,:)

        if (taskid .NE. MASTER) then   
            strt_col = start_col +kny/2
            ed_col = end_col +kny/2

            !print *, 'before filtering at point'
            do j = strt_col, ed_col
                !write(*,'(A13 I3 I4 A3 I4)') 'at processor',taskid, j-strt_col+1,'of',ed_col-strt_col+1  
                do i = knx/2 + 1, knx/2+1 + nxpo
                    FILT_FEILD(i,j) = filterAtPoint(i,j, UN_FILT_FEILD(:,:))
                enddo
            enddo

            row_size = nxpo + knx-1
            data_size = row_size * (ed_col - strt_col +1)

            dest = MASTER

            call MPI_SEND( strt_col, 1, MPI_INTEGER, MASTER, startCol_tag + taskid, &
            &                 MPI_COMM_WORLD, i_err )

            call MPI_SEND( data_size, 1, MPI_INTEGER, MASTER, dataSize_tag + taskid, &
            &                 MPI_COMM_WORLD, i_err )

            call MPI_SEND( FILT_FEILD(1,strt_col), data_size, MPI_REAL, dest, &
            & OL_UVEL_tag + taskid, MPI_COMM_WORLD, status, i_err )

        
        else
            do source = 1, numworkers

                call MPI_RECV( strt_col, 1, MPI_INTEGER, source, startCol_tag + source, &
                &                 MPI_COMM_WORLD, status, i_err )

                if (i_err .NE. MPI_SUCCESS ) then 
                    print *, 'from taskid',taskid,' recv start_col'
                    call MPI_ABORT(MPI_COMM_WORLD, i_err)
                endif

                call MPI_RECV( data_size, 1, MPI_INTEGER, source, dataSize_tag + source, &
                &                 MPI_COMM_WORLD, status, i_err )

                if (i_err .NE. MPI_SUCCESS ) then 
                    print *, 'from taskid',taskid,' recv data_size'
                    call MPI_ABORT(MPI_COMM_WORLD, i_err)
                endif

                call MPI_RECV( FILT_FEILD(1,strt_col), data_size, MPI_REAL, source, &
                & OL_UVEL_tag + source, MPI_COMM_WORLD, status, i_err )

                if (i_err .NE. MPI_SUCCESS ) then 
                    print *, 'from taskid',taskid,' recv OL_UVEL'
                    call MPI_ABORT(MPI_COMM_WORLD, i_err)
                endif

                !print *, 'RECEIVED from Source', source
            
            enddo

        endif

        call MPI_BARRIER(MPI_COMM_WORLD, i_err)

        outField = FILT_FEILD(is:ie, js:je)
        DEALLOCATE(UN_FILT_FEILD, FILT_FEILD)


    end subroutine

    function filterAtPoint(index_i,index_j, UN_FILT_FEILD) result(fltPoint)
        REAL(kind = r8), INTENT(IN), DIMENSION(:,:) :: UN_FILT_FEILD
        INTEGER(kind=i4), INTENT(IN) :: index_i, index_j
        REAL(kind = r8) :: fltPoint

        INTEGER(kind=i4) :: knx, kny, is, ie, js, je
        REAL(kind = r8), ALLOCATABLE, DIMENSION(:,:) :: uf_WORK

        knx = size(kernel,1)
        kny = size(kernel,2)

        allocate(uf_WORK(knx, kny))

        is = index_i - knx/2
        ie = index_i + knx/2

        js = index_j - kny/2
        je = index_j + kny/2

        uf_WORK = UN_FILT_FEILD(is:ie,js:je) 
        fltPoint = sum(uf_WORK * kernel * dArea)/sum(kernel * dArea)

        fltPoint = UN_FILT_FEILD(index_i, index_j)

        DEALLOCATE(uf_WORK)
    
    end function

end module filter