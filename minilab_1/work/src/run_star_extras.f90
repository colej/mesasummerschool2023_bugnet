! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful,
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************

      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      use gyre_lib


      use num_lib
      use utils_lib

      implicit none

      real(dp):: delta_omega_g, omega_max, Delta_Pi1
      real(dp), allocatable, save :: frequencies(:,:), inertias(:,:)
      logical, save :: gyre_has_run
      integer, save :: nmax, nmax_prev
      real(dp), save :: chi2_old, chi2


      ! these routines are called by the standard run_star check_model
      contains

        subroutine extras_controls(id, ierr)
           integer, intent(in) :: id
           integer, intent(out) :: ierr
           type (star_info), pointer :: s
           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return

           ! this is the place to set any procedure pointers you want to change
           ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)


           ! the extras functions in this file will not be called
           ! unless you set their function pointers as done below.
           ! otherwise we use a null_ version which does nothing (except warn).

           s% extras_startup => extras_startup
           s% extras_start_step => extras_start_step
           s% extras_check_model => extras_check_model
           s% extras_finish_step => extras_finish_step
           s% extras_after_evolve => extras_after_evolve
           s% how_many_extra_history_columns => how_many_extra_history_columns
           s% data_for_extra_history_columns => data_for_extra_history_columns
           s% how_many_extra_profile_columns => how_many_extra_profile_columns
           s% data_for_extra_profile_columns => data_for_extra_profile_columns

           s% how_many_extra_history_header_items => how_many_extra_history_header_items
           s% data_for_extra_history_header_items => data_for_extra_history_header_items
           s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
           s% data_for_extra_profile_header_items => data_for_extra_profile_header_items


        end subroutine extras_controls


        subroutine extras_startup(id, restart, ierr)
           integer, intent(in) :: id
           logical, intent(in) :: restart
           integer, intent(out) :: ierr
           type (star_info), pointer :: s
           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return


    ! >>> Insert allocation code below

           allocate(frequencies(2,300), inertias(2,300))
           nmax = 0
           !chi2_old = 1d99


        end subroutine extras_startup


        integer function extras_start_step(id)
           integer, intent(in) :: id
           integer :: ierr
           type (star_info), pointer :: s
           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return
           extras_start_step = 0

           frequencies = 0._dp
           inertias = 0._dp
           nmax_prev = nmax
           if (s% center_h1 < 0.7) s% write_profiles_flag = .true.
        end function extras_start_step


        ! returns either keep_going, retry, or terminate.
        integer function extras_check_model(id)
           integer, intent(in) :: id
           integer :: ierr
           type (star_info), pointer :: s
           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return
           extras_check_model = keep_going
           !if (safe_log10(s% Teff) < 3.7) then
              ! stop when star hydrogen mass drops to specified level
            !  extras_check_model = terminate
             ! write(*, *) 'log Teff < 3.7'
              !return
           !end if


           ! if you want to check multiple conditions, it can be useful
           ! to set a different termination code depending on which
           ! condition was triggered.  MESA provides 9 customizeable
           ! termination codes, named t_xtra1 .. t_xtra9.  You can
           ! customize the messages that will be printed upon exit by
           ! setting the corresponding termination_code_str value.
           ! termination_code_str(t_xtra1) = 'my termination condition'

           ! by default, indicate where (in the code) MESA terminated
           if (extras_check_model == terminate) s% termination_code = t_extras_check_model
        end function extras_check_model


        integer function how_many_extra_history_columns(id)
           integer, intent(in) :: id
           integer :: ierr
           type (star_info), pointer :: s
           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return
           how_many_extra_history_columns = 0
        end function how_many_extra_history_columns


        subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
           integer, intent(in) :: id, n
           character (len=maxlen_history_column_name) :: names(n)
           real(dp) :: vals(n), integral_N, integral_N3, I, mu_0, Br_mean
           integer, intent(out) :: ierr
           type (star_info), pointer :: s
           double precision, allocatable :: brunt_N(:)
           integer :: k
           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return
           mu_0 = 4d-8*pi

           ! note: do NOT add the extras names to history_columns.list
           ! the history_columns.list is only for the built-in history column options.
           ! it must not include the new column names you are adding here.

        end subroutine data_for_extra_history_columns


        integer function how_many_extra_profile_columns(id)
           integer, intent(in) :: id
           integer :: ierr
           type (star_info), pointer :: s
           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return
           how_many_extra_profile_columns = 4
        end function how_many_extra_profile_columns


        subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
           integer, intent(in) :: id, n, nz
           character (len=maxlen_profile_column_name) :: names(n)
           real(dp) :: vals(nz,n)
           integer, intent(out) :: ierr
           type (star_info), pointer :: s
           integer :: k, ki, km
           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return

           ! note: do NOT add the extra names to profile_columns.list
           ! the profile_columns.list is only for the built-in profile column options.
           ! it must not include the new column names you are adding here.

           ! here is an example for adding a profile column
           !if (n /= 1) stop 'data_for_extra_profile_columns'
           !names(1) = 'beta'
           !do k = 1, nz
           !   vals(k,1) = s% Pgas(k)/s% P(k)
           !end do

              write (names(1),    '(A,I0)') 'freq_l', 0
              write (names(2), '(A,I0)') 'freq_l', 1
              write (names(3),    '(A,I0)') 'Enorm_l', 0
              write (names(4), '(A,I0)') 'Enorm_l', 1

              do k =1, 300
                ki = k
                if (frequencies(2, k) > 0) exit
              end do
              do k =300, 1, -1
                km = k
                if (frequencies(2, k) > 0) exit
              end do
              do k =km, s% nz
                vals(k,1) = frequencies(1, km)
                vals(k,2) = frequencies(2, km)
                vals(k,3) = safe_log10(inertias(1, km))
                vals(k,4) = safe_log10(inertias(2, km))
              end do
              do k =1, ki
                vals(k,1) = frequencies(1, ki)
                vals(k,2) = frequencies(2, ki)
                vals(k,3) = safe_log10(inertias(1, ki))
                vals(k,4) = safe_log10(inertias(2, ki))
              end do
              ! save the frequencies of the radial and dipole modes
              do k = ki, km
                if (frequencies(1, k) > 0 ) then
                  vals(k,1) = frequencies(1, k)
                  vals(k,3) = safe_log10(inertias(1, k))
                end if
                if (frequencies(2, k) > 0 ) then
                    vals(k,2) = frequencies(2, k)
                    vals(k,4) = safe_log10(inertias(2, k))
                end if
              end do

        end subroutine data_for_extra_profile_columns


        integer function how_many_extra_history_header_items(id)
           integer, intent(in) :: id
           integer :: ierr
           type (star_info), pointer :: s
           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return
           how_many_extra_history_header_items = 0
        end function how_many_extra_history_header_items


        subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
           integer, intent(in) :: id, n
           character (len=maxlen_history_column_name) :: names(n)
           real(dp) :: vals(n)
           type(star_info), pointer :: s
           integer, intent(out) :: ierr
           ierr = 0
           call star_ptr(id,s,ierr)
           if(ierr/=0) return

           ! here is an example for adding an extra history header item
           ! also set how_many_extra_history_header_items
           ! names(1) = 'mixing_length_alpha'
           ! vals(1) = s% mixing_length_alpha

        end subroutine data_for_extra_history_header_items


        integer function how_many_extra_profile_header_items(id)
           integer, intent(in) :: id
           integer :: ierr
           type (star_info), pointer :: s
           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return
           how_many_extra_profile_header_items = 0
        end function how_many_extra_profile_header_items


        subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
           integer, intent(in) :: id, n
           character (len=maxlen_profile_column_name) :: names(n)
           real(dp) :: vals(n)
           type(star_info), pointer :: s
           integer, intent(out) :: ierr
           ierr = 0
           call star_ptr(id,s,ierr)
           if(ierr/=0) return

           ! here is an example for adding an extra profile header item
           ! also set how_many_extra_profile_header_items
           ! names(1) = 'mixing_length_alpha'
           ! vals(1) = s% mixing_length_alpha

        end subroutine data_for_extra_profile_header_items


        ! returns either keep_going or terminate.
        ! note: cannot request retry; extras_check_model can do that.
        integer function extras_finish_step(id)
           integer, intent(in) :: id
           integer :: ierr
           type (star_info), pointer :: s
           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return
           extras_finish_step = keep_going

           ! to save a profile,
              ! s% need_to_save_profiles_now = .true.
           ! to update the star log,
              ! s% need_to_update_history_now = .true.

           if (safe_log10(s% Teff) < 3.7 .and. s%x_logical_ctrl(1)) call run_gyre(id, ierr)
           !chi2_old = chi2


           ! see extras_check_model for information about custom termination codes
           ! by default, indicate where (in the code) MESA terminated
           if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
        end function extras_finish_step


        subroutine extras_after_evolve(id, ierr)
           integer, intent(in) :: id
           integer, intent(out) :: ierr
           type (star_info), pointer :: s
           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return
        end subroutine extras_after_evolve

        subroutine run_gyre (id, ierr)

           integer, intent(in)  :: id
           integer, intent(out) :: ierr

           real(dp), allocatable :: global_data(:)
           real(dp), allocatable :: point_data(:,:)
           integer               :: ipar(0)
           real(dp)              :: rpar(0)

           ! Pass model data to GYRE

           call star_get_pulse_data(id, 'GYRE', .FALSE., .TRUE., .FALSE., &
                global_data, point_data, ierr)
           if (ierr /= 0) then
              print *,'Failed when calling star_get_pulse_data'
              return
           end if

           call gyre_set_model(global_data, point_data, 101)

           ! Run GYRE to get modes

           call gyre_get_modes(1, process_mode, ipar, rpar)

           gyre_has_run = .true.

        contains

           subroutine process_mode (md, ipar, rpar, retcode)

              type(mode_t), intent(in) :: md
              integer, intent(inout)   :: ipar(:)
              real(dp), intent(inout)  :: rpar(:)
              integer, intent(out)     :: retcode
              integer :: k

              type (star_info), pointer :: s
              ierr = 0
              call star_ptr(id, s, ierr)
              if (ierr /= 0) return

              ! Print out degree, radial order, mode inertia, and frequency
              print *, 'Found mode: l, m, n_p, n_g, zeta, nu = ', &
                md%id-nmax_prev, md%l, md%m, md%n_p, md%n_g, md%n_pg, REAL(md%E_norm()),REAL(md%freq('UHZ'))


              frequencies(md%l+1, md%id-nmax_prev) = REAL(md%freq('UHZ'))
              inertias(md%l+1, md%id-nmax_prev) = REAL(md%E_norm())
              nmax = md%id 
              retcode = 0
           end subroutine process_mode


        end subroutine run_gyre


      end module run_star_extras
