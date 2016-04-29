! Main math subroutines

module math
	use iso_c_binding

	! Constants
	real(kind=c_double), parameter :: pi = 3.14159265358979323846_c_double
	
	! List of beta coefficiens
	real(kind=c_double), dimension(:), allocatable :: beta_array

contains
	subroutine print_array (a)
		real(c_double), dimension(:), intent(in) :: a

		integer(c_int) :: i

		do i=1, size(a)
			write (*, "(e20.10, a1)"), a(i), ","
		end do		
	end subroutine print_array	



	subroutine calc_periodic_spline_coefs (f, u)
		! Arrays of the same size
		real(c_double), dimension(:), intent(in) :: f ! Function values
		real(c_double), dimension(:), intent(out) :: u ! Spline coefficients
	
		real(c_double), dimension(:), allocatable :: d, v, l
		integer(c_int) :: N, i, M
		integer(c_int) :: status_var
		real(c_double) :: FF
		
		real(c_double) :: ac, apN, tmp, c

		! MUST BE EVEN NUMBER
		N = size(f)

		allocate(d(N), v(N), l(N), STAT = status_var)
		if (status_var /= 0) stop "*** Cannot allocate memory ***"

		! Progonka for matrix (N-1)x(N-1)
		! (We assume there that u(N) = 0 and solve first N-1 equations)
		! (The result is stored in the d-array)
		! (Богачев Практикум на ЭВМ (v_i = 1 and u_i = 1 / l_i in those notations))
		!
		! | 4 1 0 0 0 ... 0 |
		! | 1 4 1 0 0 ... 0 |
		! | 0 1 4 1 0 ... 0 |
		! | ............... |
		! | 0 0 0 0 ... 1 4 |
	
		l(1) = 4.0d0
		d(1) = 6.0d0 * (f(N) - 2.0d0 * f(1) + f(2)) / l(1)

		do i=2, N-1
			l(i) = 4.0d0 - 1.0d0 / l(i-1)	
            FF = 6.0d0 * (f(i-1) - 2.0d0 * f(i) + f(i+1))
			d(i) = (FF - d(i-1)) / l(i)
		end do

		do i=N-2,1,-1
			d(i) = d(i) - d(i+1)/l(i)
		end do

		! Don't forget to set d(N) 
		d(N)=0
	
		! Let us return to the notations of Petrov
		! Solution for homogeneous system
		ac = 2.0d0 + sqrt(3.0d0)

		! Init v(i) as ac^i
		! We want to avoid overflow/underflow here for large N

		if (N <= 512) then
			tmp = 1.0d0
			
			do i=1,N/2
				tmp = tmp * ac
				v(i) = tmp
			end do	

			apN = tmp * tmp
		
			tmp = -1.0d0 / (1.0d0 + 1.0d0 / apN)
			do i=1, N/2
				v(i) = (v(i)/apN + 1.0d0 / v(i)) * tmp
				v(N-i) = v(i)
				tmp = -tmp
			end do
			v(N) = 1.0d0

			! Linear combination of two solutions
			c = (apN + 1.0d0) / (apN - 1.0d0) * (6.0d0 * (f(N-1) - 2.0d0 * f(N) +f(1)) - &
				& d(1) - d(N-1)) / (ac - 1.0d0 / ac)
		else
			! v(i) ~ a^(-i)
			tmp = 1.0d0

			M = min(N/2, 538)
			do i=1,M
				tmp = tmp / ac
				v(i) = tmp
				v(N-i) = tmp
			end do	
			
			! Avoid underflow. v(i) ~ 0
			do i=M+1,N/2
				v(i) = 0.0d0
				v(N-i) = 0.0d0
			end do

			v(N) = 1.0d0
			
			! Linear combination of two solutions
			c = (6.0d0 * (f(N-1) - 2.0d0 * f(N) +f(1)) - d(1) - d(N-1)) / (ac - 1.0d0 / ac)

		end if

	
		do i=1,N
			u(i) = d(i) + c * v(i)
		end do

		deallocate( d, v, l, STAT = status_var)
		if (status_var /= 0) stop "*** Cannot deallocate memory ***"
	end subroutine calc_periodic_spline_coefs



	! Calculate f'(i) through spline coefficient
	function diff_spline (fp1, fm1, up1, um1, N) result(s)
		! f_(i+1), f_(i-1)
		real(c_double) :: fp1, fm1
		
		! Corresponding spline coefficients
		real(c_double) :: up1, um1

		! Number of points
		integer(c_int) :: N

		! Result
		real(c_double) :: s

		s = N / 2.0d0 * (fp1 - fm1 - (up1 - um1) / 6.0d0)
	end function diff_spline



	! Integrate f(i) on [dzeta_(i-1), dzeta_i] through spline
	function int_spline_im1_i (f, fm1, u, um1, N) result(s)
		! f_(i), f_(i-1)
		real(c_double) :: f, fm1
		
		! Corresponding spline coefficients
		real(c_double) :: u, um1

		! Number of points
		integer(c_int) :: N

		! Result
		real(c_double) :: s

		s = (0.5d0 * (fm1 + f) - (um1+u) /24.0d0) / N
	end function int_spline_im1_i 



	function alpha (m, N) result(s)
		integer(c_int), intent(in) :: m
		integer(c_int), intent(in) :: N ! Must be even

		integer(c_int) :: i
		real(c_double) :: s
	
		if (mod(m, 2) == 0) then
			s = log(2.0d0) + 1.0d0 / N
		else
			s = log(2.0d0) - 1.0d0 / N
		end if
			
		do i=1, N/2-1
			s = s + cos(2.0d0 * pi * i * m / N) / i
		end do

		s = -s
	end function alpha



	function beta (m, N) result (s)
		integer(c_int), intent(in) :: m
		integer(c_int), intent(in) :: N ! Must be even

		real(c_double) :: s

		if (m == 0) then
			s = alpha(m, N)
		else
			s = -log(abs(sin(pi * m / N))) + alpha(m, N)
		end if	
	end function beta



	subroutine init_beta_array (N)
		integer(c_int), intent(in) :: N
		integer(c_int) :: i, status_var

		allocate(beta_array(N), STAT = status_var)
		if (status_var /= 0) stop "*** Cannot allocate memory ***"

		do i=1,N
			beta_array(i) = beta(i-1, N)
		end do

	end subroutine init_beta_array 



	subroutine deallocate_beta_array
		integer(c_int)    :: status_var

		deallocate( beta_array, STAT = status_var)
		if (status_var /= 0) stop "*** Cannot deallocate memory ***"
	end subroutine deallocate_beta_array

end module math 
