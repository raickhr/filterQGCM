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
                  startcount(4)

        real(kind =r8) :: tempField2D(nxpo,nypo)

        path_n_file = trim(adjustl(inputpath))//trim(adjustl(filename))

        print *,''
        print *, 'Opening file  ', trim(path_n_file), ' ... '

        ierr = nf_open (trim(path_n_file), nf_nowrite, file_id)
        if ( ierr /= nf_noerr )  call handle_err(ierr, 'nf_open')

        print *,''
        print *, 'Opened  ', trim(path_n_file)

        num_of_var = size(input2DOcnFields)

        do var_counter = 1,num_of_var
            varName = trim(adjustl(input2DOcnFields(var_counter)%info%fieldName))
            ierr = nf_inq_varid(file_id,trim(adjustl(varName)),var_id)

            if ( ierr /= nf_noerr )  call handle_err(ierr, 'nf_inq_varid')
            startcount = (/1, 1, k_level, recDim_count/)
            endcount = (/nxpo, nypo, 1 , 1/)

            ierr = NF_GET_ATT_TEXT (file_id, var_id, 'long_name', longName)
            if ( ierr /= nf_noerr )  call handle_err(ierr, 'nf_get_att')
            input2DOcnFields(var_counter)%info%longName = longName

            ierr = NF_GET_ATT_TEXT (file_id, var_id, 'units', units)
            if ( ierr /= nf_noerr )  call handle_err(ierr, 'nf_get_att')
            input2DOcnFields(var_counter)%info%longName = longName

            ierr = nf_get_vara_real(file_id,var_id,startcount,endcount,tempField2D)
            if ( ierr /= nf_noerr )  call handle_err(ierr, 'nf_get_field')

            input2DOcnFields(var_counter)%fieldVal(:,:) = tempField2D(:,:)

            WRITE(*,'(A,A)')  '  Variable name :',trim(varName)
            WRITE(*,'(A,A)')  '  Long name     :',trim(longName)
            WRITE(*,'(A,A)')  '  Units         :',trim(units)
            WRITE(*,'(A,I3)') '  Variable id   :',var_id
            print *,''

        end do

        print *,''
        print *, 'Closing file  ', trim(path_n_file), ' ... '

        ierr = nf_close(file_id)
        

    end subroutine
end module netCDFio