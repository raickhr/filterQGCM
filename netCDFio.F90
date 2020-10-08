module netCDFio
    use kinds
    use configMod
    use gridMod
    use fields
    use ncdf_wrapper
    implicit none
    
    contains

    subroutine readInputOcnVars2D_P(fileCounter, recDim_count)
        INTEGER(kind = i4) , INTENT(IN) :: fileCounter, recDim_count

        CHARACTER (len = char_len_long) :: inputpath
        CHARACTER (len = char_len) :: filename

        character(len=char_len):: path_n_file

        CHARACTER(len=char_len_short) :: varName, &
                                         longName, &
                                         units

        integer(kind = i4):: k_level, nxpo, nypo
        
        integer(kind=i4):: num_of_var, &
                  var_counter, &
                  var_id, file_id, field_id, &
                  ierr, &
                  endcount(4), &
                  startcount(4), &
                  endcount2(3), &
                  startcount2(3)


        real(kind =r8), ALLOCATABLE :: tempField2D(:,:,:), &
                                       xpo(:), &
                                       ypo(:)

        real(kind = r8) :: timeVal

        call setOcnPgridXYsizeto(nxpo, nypo)
        k_level = 1

        ALLOCATE(tempField2D(nxpo, nypo, 3), &
                 xpo(nxpo), &
                 ypo(nypo))


        inputpath = trim(adjustl(getInputLoc()))
        filename = trim(adjustl(getInFileNo(fileCounter)))

        path_n_file = trim(adjustl(inputpath))//trim(adjustl(filename))

        print *,''
        print *, 'Opening file  ', trim(path_n_file), ' ... '

        ierr = nf_open (trim(path_n_file), nf_nowrite, file_id)
        if ( ierr /= nf_noerr )  call handle_err(ierr, 'nf_open')

        print *,''
        print *, 'Opened  ', trim(path_n_file)

        varName = 'xp'
        ierr = nf_inq_varid(file_id,trim(adjustl(varName)),var_id)
        if ( ierr /= nf_noerr )  call handle_err(ierr, 'nf_inq_varid: xpo')

        ierr = nf_get_vara_real(file_id,var_id,(/1/),(/nxpo/),xpo)
        if ( ierr /= nf_noerr )  call handle_err(ierr, 'nf_get_vara_real: time Val')

        varName = 'yp'
        ierr = nf_inq_varid(file_id,trim(adjustl(varName)),var_id)
        if ( ierr /= nf_noerr )  call handle_err(ierr, 'nf_inq_varid: ypo')

        ierr = nf_get_vara_real(file_id,var_id,(/1/),(/nypo/),ypo)
        if ( ierr /= nf_noerr )  call handle_err(ierr, 'nf_get_vara_real: time Val')

        varName = 'time'
        ierr = nf_inq_varid(file_id,trim(adjustl(varName)),var_id)
        if ( ierr /= nf_noerr )  call handle_err(ierr, 'nf_inq_varid: time')

        ierr = NF_GET_ATT_TEXT (file_id, var_id, 'units', timeUnits)
        if ( ierr /= nf_noerr )  call handle_err(ierr, 'nf_get_att: time units')

        ierr = nf_get_vara_real(file_id,var_id,(/recDim_count/),(/1/),timeVal)
        if ( ierr /= nf_noerr )  call handle_err(ierr, 'nf_get_vara_real: time Val')

        call setCurrTimeVal(timeVal)
        
        num_of_var = numVarsToRead()

        do var_counter = 1,num_of_var
            varName = trim(adjustl(getVarNameNo(var_counter)))
            ierr = nf_inq_varid(file_id,trim(adjustl(varName)),var_id)

            if ( ierr /= nf_noerr )  call handle_err(ierr, 'nf_inq_varid')
            startcount = (/1, 1, k_level, recDim_count/)
            endcount = (/nxpo, nypo, 1 , 1/)

            startcount2 = (/1, 1, recDim_count/)
            endcount2 = (/nxpo, nypo, 1/)

            ierr = NF_GET_ATT_TEXT (file_id, var_id, 'long_name', longName)
            if ( ierr /= nf_noerr )  call handle_err(ierr, 'nf_get_att')
            input2DOcnFields(var_counter)%longName = longName

            ierr = NF_GET_ATT_TEXT (file_id, var_id, 'units', units)
            if ( ierr /= nf_noerr )  call handle_err(ierr, 'nf_get_att')
            input2DOcnFields(var_counter)%longName = longName

            if (var_counter .EQ. 1) then
              ierr = nf_get_vara_real(file_id, var_id, startcount, endcount, tempField2D(:,:,var_counter))
            else
              ierr = nf_get_vara_real(file_id, var_id, startcount2, endcount2, tempField2D(:,:,var_counter))
            endif

            if ( ierr /= nf_noerr )  call handle_err(ierr, 'nf_get_field')

            !input2DOcnFields(var_counter)%fieldVal(:,:) = tempField2D(:,:)

            WRITE(*,'(A,A)')  '  Variable name :',trim(varName)
            WRITE(*,'(A,A)')  '  Long name     :',trim(longName)
            WRITE(*,'(A,A)')  '  Units         :',trim(units)
            WRITE(*,'(A,I3)') '  Variable id   :',var_id
            print *,''

        end do

        print *,''
        print *, 'Closing file  ', trim(path_n_file), ' ... '

        ierr = nf_close(file_id)

        ! print *, 'after reading TAUX(500,600)', tempField2D(500,600,2) * 1000  !! * density

        call saveReadInputFields(nxpo, nypo, tempField2D(:,:,1), &
                                 tempField2D(:,:,2), &
                                 tempField2D(:,:,3))

        call saveReadXpoYpo(xpo,ypo)

        DEALLOCATE(tempField2D, xpo, ypo)


    end subroutine

    subroutine writeOcnOutPut_Pfields(recDim, filterLen)
      ! outPath,  &
      !                              outFileName)
      INTEGER(kind = i4), INTENT(IN) :: recDim
      REAL(kind = r4), INTENT(IN) :: filterLen 

      CHARACTER(len = 4 ) :: str_recDim, str_filterLen

      CHARACTER (len = char_len_long) :: outPath, fileWithPath

      CHARACTER(len=char_len) :: filename, &
                                 varName, &
                                 varLongName, &
                                 varUnits

      integer(kind = i4) :: counter, nxpo, nypo, nf_status, i_err

      integer(kind = i4) :: f_id, &
                            x_dim_id, &
                            y_dim_id, &
                            time_dim_id, &
                            coord_ids(3), &
                            xpo_id, &
                            ypo_id, &
                            time_id

      integer ::field_id(7)
      REAL(kind=r8), ALLOCATABLE, DIMENSION(:,:,:):: wfield

      REAL(kind=r8), ALLOCATABLE, DIMENSION(:):: xpo, ypo

      REAL(kind=r8) :: timeVal

      call setOcnPgridXYsizeto(nxpo,nypo)

      ALLOCATE(wfield(nxpo, nypo, 7), xpo(nxpo), ypo(nypo))

      call getXpoYpo(xpo(:), ypo(:))

      call getOutputFields(nxpo, nypo, wfield(:,:,1), wfield(:,:,2), &
                                     wfield(:,:,3), wfield(:,:,4), &   
                                     wfield(:,:,5))

      wfield(:,:,6) = wfield(:,:,1)*wfield(:,:,3) + wfield(:,:,2) * wfield(:,:,4)
      wfield(:,:,7) = wfield(:,:,5) - wfield(:,:,6)

      timeVal = getCurrTimeVal()

      !if (thisProc() == MASTER) then
      !        print *,'before writing xpo(500), ypo(500)', xpo(500), ypo(500)
      !        print *,'before writing TAUX(500,600)', wfield(500,600,3)
      !endif

      write(str_recDim,'(i4.4)') recDim
      write(str_filterLen,'(i4.4)') int(filterLen)

      filename = 'recNo_'//trim(str_recDim)//'.filterLen_'//trim(str_filterLen)//'.nc'

      outPath = getOutputLoc()

      fileWithPath = trim(adjustl(outPath))//'/'//&
                 trim(adjustl(filename))

      fileWithPath = trim(adjustl(fileWithPath))

      !-------------------------------------------------------------------
      !  open netcdf file
      !-------------------------------------------------------------------

      ! print *,'Opening to write ...',filename
   	  nf_status = nf_create(fileWithPath, nf_clobber, f_id)
   	  if (nf_status /= nf_noerr) stop 'at create file'

      !-------------------------------------------------------------------
      !  define dimensions
      !-------------------------------------------------------------------
      ! print *, 'Defining dimensions ...'

      nf_status = nf_def_dim(f_id, 'xpo', nxpo, x_dim_id)
      if (nf_status /= nf_noerr) stop 'at def x_dim'

      nf_status = nf_def_dim(f_id, 'ypo', nypo, y_dim_id)
      if (nf_status /= nf_noerr) stop 'at def y_dim'

      nf_status = nf_def_dim(f_id, 'time', 1, time_dim_id)
      if (nf_status /= nf_noerr) stop 'at def y_dim'

      ! print *, 'Dimensions defined'

      coord_ids(1)=x_dim_id
      coord_ids(2)=y_dim_id
      coord_ids(3)=time_dim_id


      ! print *,''
      ! print *, 'Defining field variable'

      nf_status = nf_def_var(f_id, &
                          'xpo', &
                          nf_float, 1, coord_ids(1), &
                          xpo_id)
      if (nf_status /= nf_noerr) stop 'at field var def: xpo'

      nf_status = nf_put_att_text(f_id, &
                               xpo_id, &
                               "units",len("km") , &
                                "km" )
      if (nf_status /= nf_noerr) stop 'at put field attribute units: xpo'

      nf_status = nf_def_var(f_id, &
                          'ypo', &
                          nf_float, 1, coord_ids(2), &
                          ypo_id)
      if (nf_status /= nf_noerr) stop 'at field var def: ypo'

      nf_status = nf_put_att_text(f_id, &
                               ypo_id, &
                               "units",len("km") , &
                                "km" )
      if (nf_status /= nf_noerr) stop 'at put field attribute units: ypo'

      nf_status = nf_def_var(f_id, &
                          'time', &
                          nf_float, 1, coord_ids(3), &
                          time_id)
      if (nf_status /= nf_noerr) stop 'at field var def: time'

      nf_status = nf_put_att_text(f_id, &
                               time_id, &
                               "units",len(trim(adjustl(timeUnits))) , &
                                timeUnits )
      if (nf_status /= nf_noerr) stop 'at put field attribute units: time'
        

      do counter = 1, 7

        varName = trim(adjustl(output2DOcnFields(counter)%fieldName))
        varLongName = trim(adjustl(output2DOcnFields(counter)%longName))
        varUnits = trim(adjustl(output2DOcnFields(counter)%units))

        nf_status = nf_def_var(f_id, &
                            varName, &
                            nf_float, 3, coord_ids(1:3), &
                            field_id(counter))
        if (nf_status /= nf_noerr) stop 'at field var def'

        nf_status = nf_put_att_text(f_id, &
                                 field_id(counter), &
                                 "units", len(trim(adjustl(varUnits))), &
                                  varUnits )
        if (nf_status /= nf_noerr) stop 'at put field attribute units'

        nf_status = nf_put_att_text(f_id, &
                                 field_id(counter), &
                                 "long_name", len(trim(adjustl(varLongName))), &
                                 varLongName)                         
        if (nf_status /= nf_noerr) stop 'at put field attribute long_name'



      enddo ! define fields

      nf_status = nf_enddef(f_id)
      if (nf_status /= nf_noerr) stop 'at enddef'

      !print *, 'Defining variables SUCCESS...'

      !-------------------------------------------------------------------
      !  start writing the file
      !-------------------------------------------------------------------
      nf_status = nf_put_var(f_id, &
                          xpo_id, &
                          xpo)

      if (nf_status /= nf_noerr) stop 'at writing xpo'

      nf_status = nf_put_var(f_id, &
                          ypo_id, &
                          ypo)
        if (nf_status /= nf_noerr) stop 'at writing ypo'

      nf_status = nf_put_var(f_id, &
                          time_id, &
                          timeVal )
        if (nf_status /= nf_noerr) stop 'at writing timeVal'

    ! print *, 'xpo ypo and time written ... '
    
      do counter = 1, 7
        nf_status = nf_put_var(f_id, &
                            field_id(counter), &
                            ! (/1,1,1/), (/nxpo, nypo, 1/), &
                             wfield(:,:,counter))
        if (nf_status /= nf_noerr) stop 'at put var'
      enddo ! write feilds


      print *, 'Written file ', trim(fileWithPath) 
      nf_status = nf_close(f_id)
        if (nf_status /= nf_noerr) stop 'at close'

    end subroutine

end module netCDFio