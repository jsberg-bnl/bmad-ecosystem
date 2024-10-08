!+
! Subroutine track1_bunch_hook (bunch, ele, err, centroid, direction, finished, bunch_track)
!
! Prototype routine that can be customized for tracking a bunch through a single element.
!
! Input:
!   bunch          -- Bunch_struct: Starting bunch position.
!   ele           -- Ele_struct: Element to track through.
!   centroid(0:)  -- coord_struct, optional: Approximate centroid orbit. Only needed if CSR is on.
!                      Hint: Calculate this before bunch tracking by tracking a single particle.
!   direction     -- integer, optional: +1 (default) -> Track forward, -1 -> Track backwards.
!   bunch_track   -- bunch_track_struct, optional: Existing tracks. If bunch_track%n_pt = -1 then
!                        Overwrite any existing track.
!
! Output:
!   bunch       -- bunch_struct: Ending bunch position.
!   err         -- Logical: Set true if there is an error. 
!                    EG: Too many particles lost for a CSR calc.
!   finished    -- logical: When set True, the standard track1_bunch code will not be called.
!   bunch_track -- bunch_track_struct, optional: track information if the tracking method does
!                        tracking step-by-step. When tracking through multiple elements, the 
!                        trajectory in an element is appended to the existing trajectory. 
!-

subroutine track1_bunch_hook (bunch, ele, err, centroid, direction, finished, bunch_track)

use e_cooling_mod
use beam_utils

implicit none

type (bunch_struct), target :: bunch
type (ele_struct), target :: ele
type (coord_struct), optional :: centroid(0:)
type (bunch_track_struct), optional :: bunch_track
type (coord_struct), pointer :: p_init, p_final
type (bunch_params_struct), target :: bunch_params
integer ib, ia
real (rp) :: delay
real (rp) :: rand1
real (rp) :: A_xm, A_ym, A_xk, A_yk
real (rp) :: k_xm, k_ym, k_xk, k_yk
real (rp) :: lambda_xm, lambda_ym, lambda_xk, lambda_yk
real (rp) :: D_h_xk, D_h_yk
real (rp) :: D_e11_xk, D_e11_yk
real (rp) :: D_e12_xk, D_e12_yk
real (rp) :: D_e22_xk, D_e22_yk
real (rp) :: A_m, A_k
real (rp) :: k_m, k_k
real (rp) :: lambda_m, lambda_k
real (rp) :: sigma_z, sigma_delta
real (rp) :: emit_x_init, emit_y_init, sigma_delta_init
real (rp) :: emit_x_final, emit_y_final, sigma_delta_final
real (rp) :: A, k, lambda
real (rp) :: D_h, D_e11, D_e12, D_e22, D_e
real (rp) :: sigma_sup, local_density_e_rel
real (rp) :: densities(100), bin_edges(101)
real (rp) :: local_h_density, bin_width, local_e_density
real (rp) :: num_p_per_mark
real (rp) :: beta_x, alpha_x, gamma_x, D_x, Dp_x
real (rp) :: beta_y, alpha_y, gamma_y, D_y, Dp_y
real (rp) :: diff_rate, diff_rate_x, diff_rate_y
real (rp) :: gamma_rel
real (rp) :: num_sigmas, frac_sigmas
integer :: num_bins

integer, optional :: direction
logical err, is_ok, finished, error
character(*), parameter :: r_name = 'track1_bunch_hook'


num_sigmas = ec%wd%num_sigmas  ! For binning h+ longitudinally
frac_sigmas = ec%wd%frac_sigmas
num_bins = nint(2*num_sigmas*frac_sigmas)

! Number of protons represented by single marker.
num_p_per_mark = ec%wd%p_per_bunch/size(bunch%particle)

! Convert e- bunch length into parameter for use in supergaussian
sigma_sup = ec%wd%supergaussian_order*sqrt(pi)/(gamma(1./(ec%wd%supergaussian_order*2.)))*ec%wd%sigma_ze

finished = .false.

! For the modulator element, simply save the initial hadron 6D distribution at the center.
! Track through as a drift.
if (associated(ec_com%input_ele, ele)) then
  call out_io (s_info$, r_name, 'This is the input element: ' // ele_full_name(ele))

  ! Track through the first half of the modulator
  do ib = 1, size(bunch%particle)
    call track_a_drift(bunch%particle(ib), ele%value(l$)/2.)
  enddo

  ! Get the bunch parameters at the center of the modulator.
  ec_com%input_bunch = bunch

  ! Track through the second half of the modulator
  do ib = 1, size(bunch%particle)
    call track_a_drift(bunch%particle(ib), ele%value(l$)/2.)
  enddo

  finished = .true.

  return

! For the kicker, run the interpolation and apply energy kicks.
elseif (associated(ec_com%output_ele, ele)) then
  call out_io (s_info$, r_name, 'This is the output element: ' // ele_full_name(ele))


  ! Get relativistic gamma
  gamma_rel = ele%value(E_tot$) / mass_of(ele%ref_species)

  ! Get bunch length and momentum spread at kicker.
  ! Also, track through first half of the kicker.
  sigma_z = 0
  sigma_delta = 0

  do ib = 1, size(bunch%particle)
    p_final => bunch%particle(ib)
    call track_a_drift(p_final, ec_com%output_ele%value(l$)/2.)
    sigma_z = sigma_z + p_final%vec(5)*p_final%vec(5)/size(bunch%particle)
    sigma_delta = sigma_delta + p_final%vec(6)*p_final%vec(6)/size(bunch%particle)
  enddo
  sigma_z = sqrt(sigma_z)
  sigma_delta = sqrt(sigma_delta)
  sigma_delta_init = sigma_delta

  ! Get initial emittance and Courant-Snyder parameters.
  call calc_bunch_params (bunch, bunch_params, error)
  emit_x_init = bunch_params%a%emit
  emit_y_init = bunch_params%b%emit

  beta_x = bunch_params%a%beta
  alpha_x = bunch_params%a%alpha
  gamma_x = bunch_params%a%gamma
  D_x = bunch_params%x%eta
  Dp_x = bunch_params%x%etap

  beta_y = bunch_params%b%beta
  alpha_y = bunch_params%b%alpha
  gamma_y = bunch_params%b%gamma
  D_y = bunch_params%y%eta
  Dp_y = bunch_params%y%etap


  ! Get densities of hadrons longitudinally.

  bin_width = sigma_z/frac_sigmas

  ! Create longitudinal bins.
  do ia = 1, num_bins
   densities(ia) = 0
   bin_edges(ia) = bin_width*(ia-1.) - bin_width*num_bins/2.
  enddo
  ! Define end of final bin.
  bin_edges(num_bins+1) = bin_width*num_bins/2.

  ! Fill said bins.
  do ib = 1, size(ec_com%input_bunch%particle)
    p_final => bunch%particle(ib)
    do ia = 1, num_bins
      if(p_final%vec(5) < bin_edges(ia+1) .and. p_final%vec(5) >= bin_edges(ia)) then
        densities(ia) = densities(ia) + num_p_per_mark/bin_width
      endif
    enddo
  enddo


  ! Apply cooling and diffusion.
  diff_rate = 0. ! Analytic diffusion rate.
  do ib = 1, size(ec_com%input_bunch%particle)
    p_init => ec_com%input_bunch%particle(ib)
    p_final => bunch%particle(ib)
    ! Get delay of each particle.
    ! Positive values mean shifts toward the head of the bunch.
    ! Delay evaluated between centers of modulator and kicker elements.
    delay = p_final%vec(5) - p_init%vec(5)

    A_m = 0
    A_k = 0
    k_m = 0
    k_k = 0
    lambda_m = 0
    lambda_k = 0
    D_h = 0
    D_e11 = 0
    D_e12 = 0
    D_e22 = 0

    ! Get fit parameters from splines.
    ! Average these over the length of the modulator and kicker as appropriate.
    do ia = 1, nint(ec_com%input_ele%value(num_steps$))
      call track_a_drift(p_init, -ec_com%input_ele%value(l$)/2. + ec_com%input_ele%value(l$)/ec_com%input_ele%value(num_steps$)*(ia-0.5))

      call spline_evaluate(ec%wd%xm%A_spl, abs(p_init%vec(1)), is_ok, A_xm)
      call spline_evaluate(ec%wd%xm%k_spl, abs(p_init%vec(1)), is_ok, k_xm)
      call spline_evaluate(ec%wd%xm%lambda_spl, abs(p_init%vec(1)), is_ok, lambda_xm)
 
      call spline_evaluate(ec%wd%ym%A_spl, abs(p_init%vec(3)), is_ok, A_ym)
      call spline_evaluate(ec%wd%ym%k_spl, abs(p_init%vec(3)), is_ok, k_ym)
      call spline_evaluate(ec%wd%ym%lambda_spl, abs(p_init%vec(3)), is_ok, lambda_ym)

      A_m = A_m + A_xm*A_ym/ec_com%input_ele%value(num_steps$)
      k_m = k_m + k_xm*k_ym/ec_com%input_ele%value(num_steps$)
      lambda_m = lambda_m + lambda_xm*lambda_ym/ec_com%input_ele%value(num_steps$)

      call track_a_drift(p_init, ec_com%input_ele%value(l$)/2. - ec_com%input_ele%value(l$)/ec_com%input_ele%value(num_steps$)*(ia-0.5))
    enddo

    ! Same for kicker.
    ! Also include diffusion info here.
    do ia = 1, nint(ec_com%output_ele%value(num_steps$))
      call track_a_drift(p_final, -ec_com%output_ele%value(l$)/2. + ec_com%output_ele%value(l$)/ec_com%output_ele%value(num_steps$)*(ia-0.5))

      call spline_evaluate(ec%wd%xk%A_spl, abs(p_final%vec(1)), is_ok, A_xk)
      call spline_evaluate(ec%wd%xk%k_spl, abs(p_final%vec(1)), is_ok, k_xk)
      call spline_evaluate(ec%wd%xk%lambda_spl, abs(p_final%vec(1)), is_ok, lambda_xk)
      call spline_evaluate(ec%wd%xk%D_h_spl, abs(p_final%vec(1)), is_ok, D_h_xk)
      call spline_evaluate(ec%wd%xk%D_e11_spl, abs(p_final%vec(1)), is_ok, D_e11_xk)
      call spline_evaluate(ec%wd%xk%D_e12_spl, abs(p_final%vec(1)), is_ok, D_e12_xk)
      call spline_evaluate(ec%wd%xk%D_e22_spl, abs(p_final%vec(1)), is_ok, D_e22_xk)
 
      call spline_evaluate(ec%wd%yk%A_spl, abs(p_final%vec(3)), is_ok, A_yk)
      call spline_evaluate(ec%wd%yk%k_spl, abs(p_final%vec(3)), is_ok, k_yk)
      call spline_evaluate(ec%wd%yk%lambda_spl, abs(p_final%vec(3)), is_ok, lambda_yk)
      call spline_evaluate(ec%wd%yk%D_h_spl, abs(p_final%vec(3)), is_ok, D_h_yk)
      call spline_evaluate(ec%wd%yk%D_e11_spl, abs(p_final%vec(3)), is_ok, D_e11_yk)
      call spline_evaluate(ec%wd%yk%D_e12_spl, abs(p_final%vec(3)), is_ok, D_e12_yk)
      call spline_evaluate(ec%wd%yk%D_e22_spl, abs(p_final%vec(3)), is_ok, D_e22_yk)

      A_k = A_k + A_xk*A_yk/ec_com%output_ele%value(num_steps$)
      k_k = k_k + k_xk*k_yk/ec_com%output_ele%value(num_steps$)
      lambda_k = lambda_k + lambda_xk*lambda_yk/ec_com%output_ele%value(num_steps$)
      D_h = D_h + D_h_xk*D_h_yk/ec_com%output_ele%value(num_steps$)
      D_e11 = D_e11 + D_e11_xk*D_e11_yk/ec_com%output_ele%value(num_steps$)
      D_e12 = D_e12 + D_e12_xk*D_e12_yk/ec_com%output_ele%value(num_steps$)
      D_e22 = D_e22 + D_e22_xk*D_e22_yk/ec_com%output_ele%value(num_steps$)

      call track_a_drift(p_final, ec_com%output_ele%value(l$)/2. - ec_com%output_ele%value(l$)/ec_com%output_ele%value(num_steps$)*(ia-0.5))
    enddo

    ! Get out the coefficients for the cooling and diffusion functions.
    A = A_m*A_k
    k = k_m*k_k
    lambda = lambda_m*lambda_k

    ! Get local h+ density near a given hadron.
    ! 0 if outside the longitudinal binning
    ! (won't see the e- anyway.)
    local_h_density = 0
    do ia = 1, num_bins
      if(p_final%vec(5) < bin_edges(ia+1) .and. p_final%vec(5) >= bin_edges(ia)) then
        local_h_density = densities(ia)
      endif
    enddo

    ! Get electron density assuming a supergaussian distribution.
    local_density_e_rel = exp(-((p_final%vec(5)/sigma_sup)**2 / 2.)**(ec%wd%supergaussian_order))
    local_e_density = local_density_e_rel * ec%wd%Ie_peak/c_light/e_charge

    ! Correct for wake reduction if proton not at center of e- bunch.
    A = A*local_density_e_rel**2*(sin(ec%wd%phi_avg*sqrt(local_density_e_rel))/sin(ec%wd%phi_avg))**2

    D_h = D_h*local_density_e_rel**4*(sin(ec%wd%phi_avg*sqrt(local_density_e_rel))/sin(ec%wd%phi_avg))**4
    D_e11 = D_e11*local_density_e_rel**4*(sin(ec%wd%phi_avg*sqrt(local_density_e_rel))/sin(ec%wd%phi_avg))**4
    D_e12 = D_e12*local_density_e_rel**3*(sin(ec%wd%phi_avg*sqrt(local_density_e_rel))/sin(ec%wd%phi_avg))**4
    D_e22 = D_e22*local_density_e_rel**2*(sin(ec%wd%phi_avg*sqrt(local_density_e_rel))/sin(ec%wd%phi_avg))**4

    ! Correct for off-momentum proton slipping relative to electrons in modulator and kicker.
    A = A*(sin(ec%wd%off_E_reduction*p_final%vec(6)*ec_com%input_ele%value(l$)/2./gamma_rel**2.)\
         /(ec%wd%off_E_reduction*p_final%vec(6)*ec_com%input_ele%value(l$)/2./gamma_rel**2.))\
         *(sin(ec%wd%off_E_reduction*p_final%vec(6)*ec_com%output_ele%value(l$)/2./gamma_rel**2.)\
         /(ec%wd%off_E_reduction*p_final%vec(6)*ec_com%output_ele%value(l$)/2./gamma_rel**2.))

    ! Analytic diffusion rate.
    diff_rate = diff_rate + (D_h*local_h_density + (D_e11+D_e12+D_e22)*local_e_density)/size(bunch%particle)

    ! Diffusion scales with local particle density, and turns-per-timestep.
    D_h = sqrt(D_h*local_h_density*ec%wd%turns_per_step)
    D_e = sqrt(max(D_e11 + D_e12 + D_e22, 0.)*local_e_density*ec%wd%turns_per_step)

    ! Apply coherent kick, scaling with turns per timestep.
    p_final%vec(6) = p_final%vec(6) + A*sin(k*delay)*exp(-delay*delay/2./lambda/lambda)*ec%wd%turns_per_step

    ! Apply hadron diffusion.
    call ran_gauss (rand1)
    p_final%vec(6) = p_final%vec(6) + rand1*D_h

    ! Apply electron diffusion.
    call ran_gauss (rand1)
    p_final%vec(6) = p_final%vec(6) + rand1*D_e

  enddo


  ! Get final momentum spread.
  sigma_delta = 0
  do ib = 1, size(bunch%particle)
    p_final => bunch%particle(ib)
    sigma_delta = sigma_delta + p_final%vec(6)*p_final%vec(6)/size(bunch%particle)
  enddo
  sigma_delta = sqrt(sigma_delta)
  sigma_delta_final = sigma_delta

  ! Also get final emittance, and Courant-Snyder parameters.
  call calc_bunch_params (bunch, bunch_params, error)
  emit_x_final = bunch_params%a%emit
  emit_y_final = bunch_params%b%emit

  beta_x = bunch_params%a%beta
  alpha_x = bunch_params%a%alpha
  gamma_x = bunch_params%a%gamma
  D_x = bunch_params%x%eta
  Dp_x = bunch_params%x%etap

  beta_y = bunch_params%b%beta
  alpha_y = bunch_params%b%alpha
  gamma_y = bunch_params%b%gamma
  D_y = bunch_params%y%eta
  Dp_y = bunch_params%y%etap

  print *, "Final Twiss parameters:"
  ! Also show analytic diffusion rates alone.
  print '(5es12.4)', beta_x, alpha_x, gamma_x, D_x, Dp_x
  print '(5es12.4)', beta_y, alpha_y, gamma_y, D_y, Dp_y
  print *, "H Invariants:"
  print '(1es12.4)', 0.5*(beta_x*Dp_x**2. + 2.*alpha_x*D_x*Dp_x + gamma_x*D_x**2.)
  print '(1es12.4)', 0.5*(beta_y*Dp_y**2. + 2.*alpha_y*D_y*Dp_y + gamma_y*D_y**2.)

  ! Get heating rates in transverse planes.
  diff_rate_x = diff_rate*0.5*(beta_x*Dp_x**2. + 2.*alpha_x*D_x*Dp_x + gamma_x*D_x**2.)/emit_x_init
  diff_rate_y = diff_rate*0.5*(beta_y*Dp_y**2. + 2.*alpha_y*D_y*Dp_y + gamma_y*D_y**2.)/emit_y_init

  ! Show the initial and final emittances.
  print *, "Initial (top) and final (bottom) a and b mode emittances:"
  print '(2es16.8)', emit_x_init, emit_x_final
  print '(2es16.8)', emit_y_init, emit_y_final

  ! Print fractional change in emittance for one turn.
  ! Assumes longitudinal emittance is 2x squared energy spread, to account for half of the emittance in z^2 term.
  print *, "Fractional effective cooling per turn in x, y, z:"
  print '(3es12.4)', (emit_x_final - emit_x_init)/emit_x_init/ec%wd%turns_per_step,\
      (emit_y_final -  emit_y_init)/emit_y_init/ec%wd%turns_per_step,\
      (sigma_delta_final**2. - sigma_delta_init**2.)/sigma_delta_init**2/2./ec%wd%turns_per_step

  print *
  print *, "Fractional pure heating per turn in x, y, z:"
  ! Also show analytic diffusion rates alone.
  print '(3es12.4)', diff_rate_x, diff_rate_y, diff_rate/sigma_delta**2./2.

  ! Get to end of the kicker.
  do ib = 1, size(bunch%particle)
    p_final => bunch%particle(ib)
    call track_a_drift(p_final, ec_com%output_ele%value(l$)/2.)
  enddo

  finished = .true.

  return

else
  call out_io (s_info$, r_name, 'This element is not part of the cooling circut: ', ele_full_name(ele))
endif

end subroutine
