module netCDFio
    use kinds
    use gridMod
    use fields
    use ncdf_wrapper
    implicit none
    
    contains

    subroutine readInputOcnVars2D_P(inputpath, &
                                    filename, & 
                                    k_level, &
                                    recDim_count)

        character(len=*), intent(in) :: inputpath,filename
        integer, intent(in) :: k_level, recDim_count

        character(len=char_len):: path_n_file

        CHARACTER(len=char_len_short) :: varName, &
                                         longName, &
                                         units
        
        integer(kind=i4):: num_of_var, &
                  var_counter, &
                  var_id, file_id, field_id, &
                  ierr, &
                  endcount(4), &
                  startcount(4), &
                  endcount2(3), &
                  startcount2(3)


        real(kind =r8) :: tempField2D(nxpo,nypo,3)

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

        
        num_of_var = size(input2DOcnFields)

        do var_counter = 1,num_of_var
            varName = trim(adjustl(input2DOcnFields(var_counter)%info%fieldName))
            ierr = nf_inq_varid(file_id,trim(adjustl(varName)),var_id)

            if ( ierr /= nf_noerr )  call handle_err(ierr, 'nf_inq_varid')
            startcount = (/1, 1, k_level, recDim_count/)
            endcount = (/nxpo, nypo, 1 , 1/)

            startcount2 = (/1, 1, recDim_count/)
            endcount2 = (/nxpo, nypo, 1/)

            ierr = NF_GET_ATT_TEXT (file_id, var_id, 'long_name', longName)
            if ( ierr /= nf_noerr )  call handle_err(ierr, 'nf_get_att')
            input2DOcnFields(var_counter)%info%longName = longName

            ierr = NF_GET_ATT_TEXT (file_id, var_id, 'units', units)
            if ( ierr /= nf_noerr )  call handle_err(ierr, 'nf_get_att')
            input2DOcnFields(var_counter)%info%longName = longName

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

        Pressure(:,:) = tempField2D(:,:,1)
        TAUX(:,:) = tempField2D(:,:,2)
        TAUY(:,:) = tempField2D(:,:,3)
        

        

    end subroutine

    subroutine writeOcnOutPut_Pfields(outPath,  &
                                   outFileName)

      integer:: counter
      character(len=*),intent(in)::outPath, outFileName

      CHARACTER(len=char_len) :: filename, &
                                 varName, &
                                 varLongName, &
                                 varUnits

      integer :: status
      integer :: f_id, &
                 x_dim_id, &
                 y_dim_id, &
                 time_dim_id, &
                 coord_ids(3), &
                 xpo_id, &
                 ypo_id, &
                 time_id

      integer ::field_id(7)
      REAL(kind=r8):: wfield(nxpo, nypo, 7)

      filename = trim(adjustl(outPath))//'/'//&
                 trim(adjustl(outFileName))

      filename = trim(adjustl(filename))

      !-------------------------------------------------------------------
      !  open netcdf file
      !-------------------------------------------------------------------

      ! print *,'Opening to write ...',filename
   	  status = nf_create(filename, nf_clobber, f_id)
   	  if (status /= nf_noerr) stop 'at create file'

      !-------------------------------------------------------------------
      !  define dimensions
      !-------------------------------------------------------------------
      ! print *, 'Defining dimensions ...'

      status = nf_def_dim(f_id, 'xpo', nxpo, x_dim_id)
      if (status /= nf_noerr) stop 'at def x_dim'

      status = nf_def_dim(f_id, 'ypo', nypo, y_dim_id)
      if (status /= nf_noerr) stop 'at def y_dim'

      status = nf_def_dim(f_id, 'time', 1, time_dim_id)
      if (status /= nf_noerr) stop 'at def y_dim'

      ! print *, 'Dimensions defined'

      coord_ids(1)=x_dim_id
      coord_ids(2)=y_dim_id
      coord_ids(3)=time_dim_id


      ! print *,''
      ! print *, 'Defining field variable'

      status = nf_def_var(f_id, &
                          'xpo', &
                          nf_float, 1, coord_ids(1), &
                          xpo_id)
      if (status /= nf_noerr) stop 'at field var def: xpo'

      status = nf_put_att_text(f_id, &
                               xpo_id, &
                               "units",len("km") , &
                                "km" )
      if (status /= nf_noerr) stop 'at put field attribute units: xpo'

      status = nf_def_var(f_id, &
                          'ypo', &
                          nf_float, 1, coord_ids(2), &
                          ypo_id)
      if (status /= nf_noerr) stop 'at field var def: ypo'

      status = nf_put_att_text(f_id, &
                               ypo_id, &
                               "units",len("km") , &
                                "km" )
      if (status /= nf_noerr) stop 'at put field attribute units: ypo'

      status = nf_def_var(f_id, &
                          'time', &
                          nf_float, 1, coord_ids(3), &
                          time_id)
      if (status /= nf_noerr) stop 'at field var def: time'

      status = nf_put_att_text(f_id, &
                               time_id, &
                               "units",len(trim(adjustl(timeUnits))) , &
                                timeUnits )
      if (status /= nf_noerr) stop 'at put field attribute units: time'
        

      do counter = 1, 7

        varName = trim(adjustl(output2DOcnFields(counter)%info%fieldName))
        varLongName = trim(adjustl(output2DOcnFields(counter)%info%longName))
        varUnits = trim(adjustl(output2DOcnFields(counter)%info%units))

        status = nf_def_var(f_id, &
                            varName, &
                            nf_float, 3, coord_ids(1:3), &
                            field_id(counter))
        if (status /= nf_noerr) stop 'at field var def'

        status = nf_put_att_text(f_id, &
                                 field_id(counter), &
                                 "units", len(trim(adjustl(varUnits))), &
                                  varUnits )
        if (status /= nf_noerr) stop 'at put field attribute units'

        status = nf_put_att_text(f_id, &
                                 field_id(counter), &
                                 "long_name", len(trim(adjustl(varLongName))), &
                                 varLongName)                         
        if (status /= nf_noerr) stop 'at put field attribute long_name'



      enddo ! define fields

      status = nf_enddef(f_id)
      if (status /= nf_noerr) stop 'at enddef'

      !print *, 'Defining variables SUCCESS...'

      !-------------------------------------------------------------------
      !  start writing the file
      !-------------------------------------------------------------------
      status = nf_put_var(f_id, &
                          xpo_id, &
                          xpo)

      if (status /= nf_noerr) stop 'at writing xpo'

      status = nf_put_var(f_id, &
                          ypo_id, &
                          ypo)
        if (status /= nf_noerr) stop 'at writing ypo'

      status = nf_put_var(f_id, &
                          time_id, &
                          timeVal )
        if (status /= nf_noerr) stop 'at writing timeVal'

    ! print *, 'xpo ypo and time written ... '

    if (taskid == MASTER) then
      wfield(:,:,1) = OL_UVEL
      wfield(:,:,2) = OL_VVEL
      wfield(:,:,3) = OL_TAUX
      wfield(:,:,4) = OL_TAUY
      wfield(:,:,5) = OL_PowerPerArea
      wfield(:,:,6) = OL_TAUX * OL_UVEL + OL_TAUY * OL_VVEL
      wfield(:,:,7) = wfield(:,:,5) - wfield(:,:,6)

    endif
    

      do counter = 1, 7
        status = nf_put_var(f_id, &
                            field_id(counter), &
                            ! (/1,1,1/), (/nxpo, nypo, 1/), &
                             wfield(:,:,counter))
        if (status /= nf_noerr) stop 'at put var'
      enddo ! write feilds


      print *, 'Written file ', filename 
      status = nf_close(f_id)
        if (status /= nf_noerr) stop 'at close'

    end subroutine

end module netCDFio