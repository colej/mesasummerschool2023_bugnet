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
      use kap_lib, only: kap_get_elect_cond_opacity


      use num_lib
      use utils_lib

      implicit none

      real(dp):: delta_omega_g, omega_max, Delta_Pi1
      real(dp), allocatable, save :: frequencies(:,:), inertias(:,:)
      logical, save :: gyre_has_run
      integer, save :: nmax, nmax_prev


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
           s% other_am_mixing => TSF_dynamo

           print *
            print "('Insert delta_omega_g in nHz')"
            !read(*,*) delta_omega_g
            delta_omega_g = 126
            write(*,*)'delta_omega_g= ',delta_omega_g


        end subroutine extras_controls


        subroutine extras_startup(id, restart, ierr)
           integer, intent(in) :: id
           logical, intent(in) :: restart
           integer, intent(out) :: ierr
           type (star_info), pointer :: s
           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return

           ! Initialize GYRE

           call gyre_init('gyre_mix.in')

           ! Set constants

           call gyre_set_constant('G_GRAVITY', standard_cgrav)
           call gyre_set_constant('C_LIGHT', clight)
           call gyre_set_constant('A_RADIATION', crad)

           call gyre_set_constant('M_SUN', Msun)
           call gyre_set_constant('R_SUN', Rsun)
           call gyre_set_constant('L_SUN', Lsun)

           call gyre_set_constant('GYRE_DIR', TRIM(mesa_dir)//'/gyre/gyre')
           !call gyre_set_constant('GYRE_DIR', '/Users/joey/Documents/Software/gyre-7.0')

    ! >>> Insert allocation code below

           allocate(frequencies(2,300), inertias(2,300))
           nmax = 0


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
           if (.false. .and. s% star_mass_h1 < 0.35d0) then
              ! stop when star hydrogen mass drops to specified level
              extras_check_model = terminate
              write(*, *) 'have reached desired hydrogen mass'
              return
           end if

           if (safe_log10(s% Teff) < 3.7 .and. .false.) call run_gyre(id, ierr)

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
           how_many_extra_history_columns = 3
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

           allocate(brunt_N(s% nz))
           names(1) = 'I'
           names(2) = 'Br_mean'
           names(3) = 'Delta_Pi1'
           brunt_N = sqrt(max(0._dp,s% brunt_N2))
           integral_N3 = 0.0_dp
           integral_N = 0.0_dp
           do k = 1, s%nz - 1
             integral_N3 = integral_N3 + 0.5_dp*(brunt_N(k)**3/(s% rho(k)) + brunt_N(k+1)**3/(s% rho(k+1)))*abs(s% r(k+1) - s% r(k)) / (s% r(k))**3
             integral_N  = integral_N + 0.5_dp*(brunt_N(k) + brunt_N(k+1))*abs(s% r(k+1) - s% r(k)) / s% r(k)
           end do
           I = integral_N3 / integral_N
           vals(1) = I
           omega_max = 2 * pi * s% nu_max * 1d-6
           Br_mean = sqrt(mu_0 * (2*pi*delta_omega_g*1d-9) * omega_max**3 / I)*10 ! In kG.
           vals(2) = Br_mean
           Delta_Pi1 = (2._dp*pi**2)/integral_N / (sqrt(2._dp))
           vals(3) = Delta_Pi1
           write(*,*) 'Br_mean [kG] = ', Br_mean, 'Delta_Pi1 [s] = ', Delta_Pi1, 'nu_max [uHz] = ', s% nu_max, 'delta_nu [uHz]', s% delta_nu,   'I = ', I
           deallocate(brunt_N)




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

              !km = 90
              !do k =1, s% nz
              !  vals(k,1) = 1d6/frequencies(1, km)/86400!-k*1d-10 !MODULO(1d6/maxval(frequencies(1, :)), Delta_Pi1)
              !  vals(k,2) = 1d6/frequencies(2, km)/86400!-k*1d-10 !MODULO(1d6/maxval(frequencies(2, :)), Delta_Pi1)
              !  vals(k,3) = (1d6/frequencies(1, km) - 1d6/frequencies(1, km-1)) !maxval(frequencies(1, :))
              !  vals(k,4) = (1d6/frequencies(2, km) - 1d6/frequencies(2, km-1)) !maxval(frequencies(2, :))
              !end do
              !! save the frequencies of the radial and dipole modes
              !do k = 2, km
              !  if (frequencies(1, k) > 0 .and. frequencies(1, k-1) > 0) then
              !    vals(k-1,1) = 1d6/frequencies(1, k)/86400 !MODULO(1d6/frequencies(1, k), Delta_Pi1)
              !    vals(k-1,3)  = (1d6/frequencies(1, k) - 1d6/frequencies(1, k-1))!  = frequencies(1, k)
              !  end if
              !  if (frequencies(2, k) > 0 .and. frequencies(2, k-1) > 0) then
              !      vals(k-1,2) = 1d6/frequencies(2, k)/86400 !MODULO(1d6/frequencies(2, k), Delta_Pi1)
              !      vals(k-1,4) = (1d6/frequencies(2, k) - 1d6/frequencies(2, k-1)) !frequencies(2, k)
              !      write(*,*) 'point', k, 1d6/frequencies(2, k)/86400, (1d6/frequencies(2, k) - 1d6/frequencies(2, k-1))
              !
              !  end if
              !end do

              do k =1, 300
                ki = k
                if (frequencies(2, k) > 0) exit
              end do
              do k =300, 1, -1
                km = k
                if (frequencies(2, k) > 0) exit
              end do
              !write(*,*) 'k', ki, km
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
                    !write(*,*) 'point', k, frequencies(2, k), safe_log10(inertias(2, k))
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

           !call gyre_get_modes(0, process_mode, ipar, rpar)
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
              !    MODULO(REAL(1d6/md%freq('UHZ')), Delta_Pi1)

              !    print *, 'Found mode: l, m, n_g, omega  = ', &
              !        md%l, md%m, md%n_g, REAL(md%freq('CYC_PER_DAY'))

              frequencies(md%l+1, md%id-nmax_prev) = REAL(md%freq('UHZ'))
              inertias(md%l+1, md%id-nmax_prev) = REAL(md%E_norm())
              nmax = md%id !The mode ID in GYRE is unique throughout the run. That means if in the first call N modes where found, the mode ID of the first mode of the second call has ID = N+1.
              ! We thus save the number of modes found such that we always store them in entries 1,...,N.
              retcode = 0
           end subroutine process_mode


        end subroutine run_gyre

        subroutine TSF_dynamo(id, ierr)
          integer, intent(in) :: id
          integer, intent(out) :: ierr
          type (star_info), pointer :: s
          real(dp) :: alpha, kap, dlnkap_dlnRho, dlnkap_dlnT, sig, eta, Kth, nu_TSF
          !real(dp), allocatable :: Kth(:), nu_TSF(:), eta(:)
          integer :: i
          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return

          !allocate(Kth(s% nz), nu_TSF(s% nz), eta(s% nz))

          alpha = 1d0
          do i = 1, s% nz
             call kap_get_elect_cond_opacity( &
                s% kap_handle, s% zbar(i), log10(s% rho(i)), log10(s% T(i)),  &
                kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
             if (ierr /= 0) return
             sig = calc_sige(s% abar(i), s% zbar(i), s% rho(i), s% T(i), s% Cp(i), kap, s% opacity(i))
             eta = calc_eta(sig)
             Kth = 16._dp/3._dp*boltz_sigma * (s% T(i))**3 / (s% opacity(i) * (s% rho(i))**2 * s% cp(i))
             nu_TSF = alpha**3 * (s% r(i))**2 * (s% Omega(i))**3 / (eta/Kth*(s% brunt_N2(i) - s% brunt_n2_composition_term(i)) &
              + s% brunt_n2_composition_term(i))
             s% am_nu_omega(i) = s% am_nu_omega(i) + nu_TSF
          end do


          !deallocate(Kth, eta)
        end subroutine TSF_dynamo

        !> Compute the magnetic diffusivity from the electric conductivity.
        !! @param sig The electrical conductivity (1/s).
        !! @param eta The magnetic diffusivity (output, cm^2/s).
        real(dp) function calc_eta(sig) result(eta)
           real(dp), intent(in) :: sig

           eta = (clight * clight / (4d0 * pi)) /sig
        end function calc_eta

        !> Computes the electrical conductivity following
        !! S.-C. YOON Oct. 10, 2003.
        !!
        !! @param abar The mean atomic mass number.
        !! @param zbar The mean atomic charge.
        !! @param rho The density (g/cm^3).
        !! @param T The temperature (K).
        !! @param Cp The specific heat at constant pressure (erg/g/K).
        !! @param kap_cond The electronic thermal opacity (cm^2/g).
        !! @param opacity The opacity (cm^2/g).
        !! @param sig The electrical conductivity (output, 1/s).
        real(dp) function calc_sige(abar, zbar, rho, T, Cp, kap_cond, opacity) result(sig)
           real(dp), intent(in) :: abar, zbar, rho, T, Cp, kap_cond, opacity
           real(dp) :: gamma, xlg, alpha, ffff, xxx, xsig1, xsig2, xsig3

           gamma = 0.2275d0*pow2(zbar) * pow(rho * 1.d-6 / abar, one_third)*1.d8/T
           xlg = log10(gamma)

           alpha = 16d0 * boltz_sigma * pow3(T) / (3d0 * opacity * pow2(rho) * Cp)

           if (xlg < -1.5d0) then
              sig = sige1(zbar,T,gamma)
           else if (xlg >= -1.5d0 .and. xlg <= 0d0) then
              xxx = (xlg + 0.75d0)*4d0/3d0
              ffff = 0.25d0*(2d0-3d0*xxx + pow3(xxx))
              xsig1 = sige1(zbar,T,gamma)
              xsig2 = sige2(T,rho,kap_cond)
              sig = (1d0-ffff)*xsig2 + ffff*xsig1
           else if (xlg > 0d0 .and. xlg < 0.5d0) then
              xsig2 = sige2(T,rho,kap_cond)
              sig = xsig2
           else if (xlg >= 0.5d0 .and. xlg < 1d0) then
              xxx = (xlg-0.75d0)*4d0
              ffff = 0.25d0*(2d0-3d0*xxx + pow3(xxx))
              xsig2 = sige2(T,rho,kap_cond)
              xsig3 = sige3(zbar,T,gamma)
              sig = (1d0-ffff)*xsig3 + ffff*xsig2
           else
              sig = sige3(zbar,T,gamma)
           endif

        end function calc_sige

        !> Computes one regime of the electrical conductivity.
        !! Written by S.-C. Yoon, Oct. 10, 2003
        !! See also Spitzer 1962 and Wendell et al. 1987, ApJ 313:284
        !! @param Z species charge
        !! @param T Temperature (K)
        !! @param xgamma The ion coupling strength (dimensionless).
        !! @param sige1 The electrical conductivity (1/s).
        real(dp) function sige1(z,t,xgamma)
           real(dp), intent(in) :: z, t, xgamma
           real(dp) :: etan, xlambda,f
           if (t >= 4.2d5) then
              f = sqrt(4.2d5/t)
           else
              f = 1.d0
           end if
           xlambda = sqrt(3d0*z*z*z)*pow(xgamma,-1.5d0)*f + 1d0
           etan = 3.d11*z*log(xlambda)*pow(t,-1.5d0)             ! magnetic diffusivity
           etan = etan/(1.d0-1.20487d0*exp(-1.0576d0*pow(z,0.347044d0))) ! correction: gammae
           sige1 = clight*clight/(pi4*etan)                    ! sigma = c^2/(4pi*eta)
        end function sige1

        !> Computes one regime of the electrical conductivity using the conductive opacity.
        !! Written by S.-C. Yoon, Oct. 10, 2003
        !! See Wendell et al. 1987, ApJ 313:284
        !! @param T Temperature (K)
        !! @param rho Temperature (g/cm^3)
        !! @param kap_cond The electronic thermal opacity (cm^2/g).
        !! @param sige2 The electrical conductivity (1/s).
        real(dp) function sige2(T,rho,kap_cond)
           real(dp), intent(in) :: t,rho,kap_cond
           sige2 = 1.11d9*T*T/(rho*kap_cond)
        end function sige2

        !> Computes the electrical conductivity in degenerate matter.
        !! Written by S.-C. Yoon, Oct. 10, 2003
        !! See Nandkumar & Pethick (1984)
        !! @param Z species charge
        !! @param T Temperature (K)
        !! @param xgamma The ion coupling strength (dimensionless).
        !! @param sige3 The electrical conductivity (1/s).
        real(dp) function sige3(z,t,xgamma)
           real(dp), intent(in) :: z, t, xgamma
           real(dp) :: rme, rm23, ctmp, xi
           rme = 8.5646d-23*t*t*t*xgamma*xgamma*xgamma/pow5(z)  ! rme = rho6/mue
           rm23 = pow(rme,2d0/3d0)
           ctmp = 1d0 + 1.018d0*rm23
           xi= sqrt(3.14159d0/3.)*log(z)/3.d0 + 2.d0*log(1.32d0+2.33d0/sqrt(xgamma))/3.d0-0.484d0*rm23/ctmp
           sige3 = 8.630d21*rme/(z*ctmp*xi)
        end function sige3


      end module run_star_extras
