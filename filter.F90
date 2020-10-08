module filter
    use kinds
    use gridMod
    use fields
    use workDiv

    implicit none
    PRIVATE
    save


    real(kind = r8), ALLOCATABLE, DIMENSION(:,:) :: kernel

    REAL(kind = r8) , ALLOCATABLE, DIMENSION(:,:):: inUVEL, &
                                                    outUVEL, &
                                                    inVVEL, &
                                                    outVVEL, &
                                                    inTAUX, &
                                                    outTAUX, &
                                                    inTAUY, &
                                                    outTAUY, &
                                                    inPPA, &
                                                    outPPA

    REAL(kind = r8), ALLOCATABLE, DIMENSION(:,:) :: uf_UVEL, &
                                                    uf_VVEL, &
                                                    uf_TAUX, &
                                                    uf_TAUY, &
                                                    uf_PPA

    real(kind = r8) :: dArea

    INTEGER(kind = i4) :: knx, kny


    contains

    subroutine makeKernel(filterLength, gridType, ocnORatmGrid )
        REAL(kind =r4), INTENT(IN) :: filterLength
        CHARACTER, INTENT(IN) :: gridType
        CHARACTER, OPTIONAL :: ocnORatmGrid
        REAL(kind = r8) :: dx, dy

        INTEGER(kind=i4) :: i, j, &  !! iterators
                            n  !! dummy number

        real(kind = r8), ALLOCATABLE, DIMENSION(:,:) :: kX2D, kY2D, &
                                                        distFromCenter, &
                                                        ell_by2, weight

        if (ALLOCATED(kernel)) then
            DEALLOCATE(kernel)
        endif

        dx = getDxpo()
        dy = getDypo()

        call MPI_BCAST(filterLength, 1, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
        call MPI_BCAST(dx, 1, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
        call MPI_BCAST(dy, 1, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)

        call MPI_BARRIER(MPI_COMM_WORLD, i_err)
                

        ! if (PRESENT(ocnORatmGrid)) then
        !     if ((ocnORatmGrid .EQ. 'A') .OR. &
        !         (ocnORatmGrid .EQ. 'a')) then
        !         if ( (gridType .EQ. 'T') .OR. &
        !            (gridType .EQ. 't')) then 
        !             dx = xta(2) - xta(1) 
        !             dy = yta(2) - yta(1)
        !         else
        !             dx = xpa(2) - xpa(1) 
        !             dy = ypa(2) - ypa(1)
        !         endif
        !     else
        !         if ( (gridType .EQ. 'T') .OR. &
        !            (gridType .EQ. 't')) then 
        !             dx = xto(2) - xto(1) 
        !             dy = yto(2) - yto(1)
        !         else
        !             dx = xpo(2) - xpo(1) 
        !             dy = ypo(2) - ypo(1)
        !         endif
        !     endif
        ! else    
        !     print*,'default dx dy'
        !     dx = xpo(2) - xpo(1) 
        !     dy = ypo(2) - ypo(1)
        ! endif

        dArea = dx*dy

        !print *,'in make kernel dx, dy= ', dx, dy

        knx = NINT(((filterLength/2.0)/dx)+1) * 2 + 1
        kny = NINT(((filterLength/2.0)/dy)+1) * 2 + 1

        print *,'in proc',thisProc(),'kernel knx, kny= ', knx, kny
        

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

        if (thisProc() .EQ. MASTER ) then !MASTER
            write(*,'(A15,I4,I4)') 'kernel size', knx, kny
            ! do dummmy =1, knx
            !     write(*,'(13F11.5)') kernel(dummmy, :)
            ! enddo
            print *, 'sum Kernel=', sum(kernel)
        endif
        
        
    end subroutine

    subroutine filterFields()
        INTEGER (kind=i4) :: strt_col, ed_col, data_size, row_size, &
                             i, j, &
                             nxpo, nypo, &
                             is, ie, js, je, &
                             startCol_tag, &
                             dataSize_tag, &
                             OL_UVEL_tag,  &
                             OL_VVEL_tag,  &
                             OL_TAUX_tag,  & 
                             OL_TAUY_tag,  &
                             OL_PPA_tag, &
                             i_err


        REAL(kind = r8), ALLOCATABLE, DIMENSION(:,:) :: UVEL, VVEL, TAUX, TAUY, PPA

        call setOcnPgridXYsizeto(nxpo, nypo)

        is = knx/2 +1
        ie = nxpo + knx/2
        js = kny/2 +1
        je = nypo + kny/2

        ALLOCATE(UVEL(nxpo, nypo), &
                 VVEL(nxpo, nypo), &
                 TAUX(nxpo, nypo), &
                 TAUY(nxpo, nypo), &
                 PPA(nxpo, nypo)) 

        allocate(inUVEL(nxpo + knx-1, nypo+kny-1), &
                 outUVEL(nxpo + knx-1, nypo+kny-1), &
                 inVVEL(nxpo + knx-1, nypo+kny-1), &
                 outVVEL(nxpo + knx-1, nypo+kny-1), &
                 inTAUX(nxpo + knx-1, nypo+kny-1), &
                 outTAUX(nxpo + knx-1, nypo+kny-1), &
                 inTAUY(nxpo + knx-1, nypo+kny-1), &
                 outTAUY(nxpo + knx-1, nypo+kny-1), &
                 inPPA(nxpo + knx-1, nypo+kny-1), &
                 outPPA(nxpo + knx-1, nypo+kny-1))

        allocate(uf_UVEL(knx, kny), &
                 uf_VVEL(knx, kny), &
                 uf_TAUX(knx, kny), &
                 uf_TAUY(knx, kny), &
                 uf_PPA(knx, kny))

        call getInputFields(nxpo, nypo, UVEL, VVEL, &
                            TAUX, TAUY, PPA)


        startCol_ta = 1000
        dataSize_ta = 2000
        OL_UVEL_tag = 3000
        OL_VVEL_tag = 4000
        OL_TAUX_tag = 5000
        OL_TAUY_tag = 6000
        OL_PPA_tag = 7000

        inUVEL(:,:) = 0
        outUVEL(:,:) = 0
        inVVEL(:,:) = 0
        outVVEL(:,:) = 0
        inTAUX(:,:) = 0
        outTAUX(:,:) = 0
        inTAUY(:,:) = 0
        outTAUY(:,:) = 0
        inPPA(:,:) = 0
        outPPA(:,:) = 0

        inUVEL(is:ie, js:je) = UVEL(:,:)
        inVVEL(is:ie, js:je) = VVEL(:,:)
        inTAUX(is:ie, js:je) = TAUX(:,:)
        inTAUY(is:ie, js:je) = TAUY(:,:)
        inPPA(is:ie, js:je) = PPA(:,:)

        if (thisProc() .NE. MASTER) then   
            strt_col = start_col +kny/2
            ed_col = end_col +kny/2

            !print *, 'before filtering at point'
            do j = strt_col, ed_col
                !write(*,'(A13 I3 I4 A3 I4)') 'at processor',thisProc(), j-strt_col+1,'of',ed_col-strt_col+1  
                do i = knx/2 + 1, knx/2+1 + nxpo
                    call filterAllFeildAtPoint(i,j)
                enddo
            enddo

            row_size = nxpo + knx-1
            data_size = row_size * (ed_col - strt_col +1)

            dest = MASTER

            call MPI_SEND( strt_col, 1, MPI_INTEGER, MASTER, startCol_tag + thisProc(), &
            &                 MPI_COMM_WORLD, i_err )
            !print *, 'from thisProc()',thisProc(),' send start_col'

            call MPI_SEND( data_size, 1, MPI_INTEGER, MASTER, dataSize_tag + thisProc(), &
            &                 MPI_COMM_WORLD, i_err )
            !print *, 'from thisProc()',thisProc(),' send data size'

            call MPI_SEND( outUVEL(1,strt_col), data_size, MPI_REAL, dest, &
            & OL_UVEL_tag + thisProc(), MPI_COMM_WORLD, status, i_err )
            !print *, 'from thisProc()',thisProc(),' send outUVEL'

            call MPI_SEND( outVVEL(1,strt_col), data_size, MPI_REAL, dest, &
            & OL_VVEL_tag + thisProc(), MPI_COMM_WORLD, status, i_err )
            !print *, 'from thisProc()',thisProc(),' send outVVEL'

            call MPI_SEND( outTAUX(1,strt_col), data_size, MPI_REAL, dest, &
            & OL_TAUX_tag + thisProc(), MPI_COMM_WORLD, status, i_err )
            !print *, 'from thisProc()',thisProc(),' send outTAUX'

            call MPI_SEND( outTAUY(1,strt_col), data_size, MPI_REAL, dest, &
            & OL_TAUY_tag + thisProc(), MPI_COMM_WORLD, status, i_err )
            !print *, 'from thisProc()',thisProc(),' send outTAUY'

            call MPI_SEND( outPPA(1,strt_col), data_size, MPI_REAL, dest, &
            & OL_PPA_tag + thisProc(), MPI_COMM_WORLD, status, i_err )
            !print *, 'from thisProc()',thisProc(),' send outPPA'
       
        else
            do source = 1, numworkers

                call MPI_RECV( strt_col, 1, MPI_INTEGER, source, startCol_tag + source, &
                &                 MPI_COMM_WORLD, status, i_err )
                !print *, 'from thisProc()', source,' recv start_col'

                call MPI_RECV( data_size, 1, MPI_INTEGER, source, dataSize_tag + source, &
                &                 MPI_COMM_WORLD, status, i_err )
                !print *, 'from thisProc()', source,' recv datasize'

                call MPI_RECV( outUVEL(1,strt_col), data_size, MPI_REAL, source, &
                & OL_UVEL_tag + source, MPI_COMM_WORLD, status, i_err )
                !print *, 'from thisProc()', source,' recv outUVEL'

                call MPI_RECV( outVVEL(1,strt_col), data_size, MPI_REAL, source, &
                & OL_VVEL_tag + source, MPI_COMM_WORLD, status, i_err )
                !print *, 'from thisProc()', source,' recv outVVEL'

                call MPI_RECV( outTAUX(1,strt_col), data_size, MPI_REAL, source, &
                & OL_TAUX_tag + source, MPI_COMM_WORLD, status, i_err )
                !print *, 'from thisProc()', source,' recv outTAUX'

                call MPI_RECV( outTAUY(1,strt_col), data_size, MPI_REAL, source, &
                & OL_TAUY_tag + source, MPI_COMM_WORLD, status, i_err )
                !print *, 'from thisProc()', source,' recv outTAUY'

                call MPI_RECV( outPPA(1,strt_col), data_size, MPI_REAL, source, &
                & OL_PPA_tag + source, MPI_COMM_WORLD, status, i_err )
                !print *, 'from thisProc()', source,' recv outPPA'

                !print *, 'RECEIVED from Source', source
            
            enddo

            is = knx/2 +1
            ie = nxpo + knx/2
            js = kny/2 +1
            je = nypo + kny/2

            call reset_FilteredFields()
            
            call saveOutputFields(nxpo, nypo, outUVEL(is:ie, js:je), &
                                  outVVEL(is:ie, js:je), &
                                  outTAUX(is:ie, js:je), &
                                  outTAUY(is:ie, js:je), &
                                  outPPA(is:ie, js:je))

        endif

        call MPI_BARRIER(MPI_COMM_WORLD, i_err)

        DEALLOCATE(inUVEL, outUVEL, inVVEL, outVVEL, inTAUX, &
                 & outTAUX, inTAUY, outTAUY, inPPA, outPPA)

        DEALLOCATE(uf_UVEL, uf_VVEL, uf_TAUX, uf_TAUY, uf_PPA)

        DEALLOCATE(UVEL, VVEL, TAUX, TAUY, PPA) 

        !PRINT *, "\n\n DEALLOCATED"

    end subroutine

    subroutine filterAllFeildAtPoint(index_i,index_j)
        INTEGER(kind=i4), INTENT(IN) :: index_i, index_j
        
        INTEGER(kind=i4) :: is, ie, js, je

        is = index_i - knx/2
        ie = index_i + knx/2

        js = index_j - kny/2
        je = index_j + kny/2

        uf_UVEL(:,:) = inUVEL(is:ie,js:je)
        uf_VVEL(:,:) = inVVEL(is:ie,js:je)
        uf_TAUX(:,:) = inTAUX(is:ie,js:je)
        uf_TAUY(:,:) = inTAUY(is:ie,js:je)
        uf_PPA(:,:) = inPPA(is:ie,js:je)

        ! outUVEL(index_i, index_j) = sum(uf_UVEL* kernel * dArea)/sum(kernel * dArea)
        ! outVVEL(index_i, index_j) = sum(uf_VVEL* kernel * dArea)/sum(kernel * dArea)
        ! outTAUX(index_i, index_j) = sum(uf_TAUX* kernel * dArea)/sum(kernel * dArea)
        ! outTAUY(index_i, index_j) = sum(uf_TAUY* kernel * dArea)/sum(kernel * dArea)
        ! outPPA(index_i, index_j) = sum(uf_PPA* kernel * dArea)/sum(kernel * dArea)

        outUVEL(index_i, index_j) = inUVEL(index_i, index_j)
        outVVEL(index_i, index_j) = inVVEL(index_i, index_j)
        outTAUX(index_i, index_j) = inTAUX(index_i, index_j)
        outTAUY(index_i, index_j) = inTAUY(index_i, index_j)
        outPPA(index_i, index_j) = inPPA(index_i, index_j)

    end subroutine

end module filter