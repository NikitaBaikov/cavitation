! General drawing subroutines

module handlers
	use iso_c_binding
	use gtk, only: gtk_init, gtk_window_new, GTK_WINDOW_TOPLEVEL, gtk_window_set_default_size, & 
	& gtk_window_set_title, gtk_container_set_border_width, gtk_table_new, gtk_container_add, &
	& TRUE, FALSE, gtk_drawing_area_new, gtk_table_attach_defaults, gtk_widget_show_all, gtk_main, &
	& gtk_widget_get_window, g_signal_connect, gtk_main_quit, gtk_table_set_row_spacings, &
	& gtk_table_set_col_spacings, c_null_ptr, GDK_SCROLL_DOWN, GDK_SCROLL_UP, gtk_widget_size_request, &
	& gtk_text_view_new, gtk_text_view_get_buffer, c_null_char, c_new_line

	use g, only: g_timeout_add

	use gdk, only: gdk_cairo_create

	use gdk_events, only: gdkeventscroll

	use cairo, only: cairo_set_line_width, cairo_arc, cairo_set_source_rgb, cairo_rectangle, &
	& cairo_new_path, cairo_fill, cairo_move_to, cairo_line_to, cairo_stroke, cairo_set_dash

	use gtk_hl
	use gtk_draw_hl
	use gdk_events

	use math

	! Let's define user type for point
	type point
		! Coordinates for drawing area
		real(c_double) :: dx, dy
			
		! Cartesian coordinates
		real(c_double) :: x, y

		! Functions in this point
		real(c_double) :: phi, psi, V, U, K, theta

		! Arrays for Adams method
		real(c_double), dimension(4) :: K_rhs, phi_rhs

		! Normal vector (FOR TESTS)
		real(c_double) :: nx, ny
	end type point	


	! My widgets
	type(c_ptr) :: main_window
	type(c_ptr) :: main_drawing_area
	type(c_ptr) :: button_pause, button_continue
	type(c_ptr) :: table
	type(c_ptr) :: view, text_buffer
	
	! Global variables
	! Initial window size
	integer(c_int) :: width = 1000
	integer(c_int) :: height = 500

	! Size of domain for drawing in cartesian coordinates
	real(c_double) :: x_min = -1.5d0
	real(c_double) :: y_min = -3.4d0
	real(c_double) :: x_max, y_max
	
	! Size of drawing area widget
	real(c_double) :: drawing_area_width, drawing_area_height

	! Zoom	
	real(c_double) :: ZOOM_VAR = 100.0d0

	integer(c_int) :: PAUSE_VAR = 0

	! tmp
	integer(c_int) :: cx = 200
	integer(c_int) :: cy = 350

	! Constants
	
	! User data	
	! Number of points (Must be even)
	integer(c_int), parameter :: N = 256

	! In the initial state bubble has ellipse shape
	! Ellipse parameters
	real(kind=c_double), parameter :: a = 1.0d0
	real(kind=c_double), parameter :: b = 1.0d0
	
	! Time step size
	real(kind=c_double), parameter :: time_delta = 5.0d-4

	! Bernoulli const
	real(kind=c_double), parameter :: bernoulli_const = 0.0d0

	! List of all points
	type(point), dimension(N) :: bubble_points_array

	! Total length of the bubble boundary
	real(c_double) :: l = 2.0d0 * pi

	! Array for Adams method
	real(c_double), dimension(4) :: l_rhs
	
	! We will use it to derive array index in Adams method
	integer(c_int) :: adams_counter = 0

	! Density
	
	! First case
	!real(c_double), dimension(N) :: density  = (/ (1.0d0, I=1,N) /)
	!real(c_double), dimension(N) :: idensity = (/ (I * 1.0d0/N, I=1,N) /)

	! Second case
	real(c_double), parameter :: dd = 0.11d0
	real(c_double), dimension(N) :: density  = &
		& (/ ((2.0d0 * sin(2.0d0 * pi * I /N)**2  + dd) / (1.0d0 + dd), I=1,N) /)
	real(c_double), dimension(N) :: idensity = &
		& (/ (I * 1.0d0 / N - sin(4.0d0 * pi * I / N) / (4.0d0 * pi * (1.0d0 + dd)), I=1,N) /)

	! Temporary arrays for splines coefficients
	real(c_double), dimension(N) :: spline_ux
	real(c_double), dimension(N) :: spline_uy
	real(c_double), dimension(N) :: spline_phi
	real(c_double), dimension(N) :: spline_psi
	real(c_double), dimension(N) :: spline_densityKV
	real(c_double), dimension(N) :: spline_V
	real(c_double), dimension(N) :: spline_UKmdVof
	real(c_double), dimension(N) :: spline_densityK
	real(c_double), dimension(N) :: spline_densityCosTheta
	real(c_double), dimension(N) :: spline_densitySinTheta

	! Temporary arrays for different purposes
	real(c_double), dimension(N) :: tmp_array, tmp_array2
	real(c_double), dimension(N,N) :: tmp_matrix
	
	! Time elapsed
	real(kind=c_double) :: time = 0.0d0

contains
	function quit (widget, event, gdata) result(ret)  bind(c)
		integer(c_int)    :: ret
		type(c_ptr), value :: widget, event, gdata

		call deallocate_beta_array()
		call gtk_main_quit()

		ret = FALSE
  	end function quit



	function A_operator (i, j) result(s)
		integer(c_int), intent(in) :: i, j ! Points indexes
		real(c_double) :: s, g_ij

		if ( i /= j ) then
			g_ij = log((bubble_points_array(i)%x - bubble_points_array(j)%x) * &
	            & (bubble_points_array(i)%x - bubble_points_array(j)%x) + &
    	        & (bubble_points_array(i)%y - bubble_points_array(j)%y) * &
        	    & (bubble_points_array(i)%y - bubble_points_array(j)%y)) * 0.5d0
		else
			g_ij = log(l * density(i) / pi)
		end if

		! List indexes starts from 1
		s = - 1.0d0 / N * (beta_array(abs(i-j)+1) + g_ij)	
	end function A_operator



	function B_operator (i, j, ux, uy) result(s)
		! Points indexes
		integer(c_int), intent(in) :: i, j 		
		
		! Splines coefficients arrays for x' and y' 
		real(c_double), dimension(:), intent(in) :: ux, uy  

		! j+1 and j-1 indexes, taking into account that function is periodic
		integer(c_int) :: indp1, indm1
		
		integer(c_int) :: k
		
		real(c_double) :: s, gn, rs, dx, dy

		if ( i /= j ) then
			rs = (bubble_points_array(i)%x - bubble_points_array(j)%x) * &
            	& (bubble_points_array(i)%x - bubble_points_array(j)%x) + &
	            & (bubble_points_array(i)%y - bubble_points_array(j)%y) * &
    	        & (bubble_points_array(i)%y - bubble_points_array(j)%y)

			indp1 = mod(j+1, N)
			if (indp1 == 0) indp1 = N
 
			indm1 = mod(j-1, N)
			if (indm1 == 0) indm1 = N

			dx = diff_spline (bubble_points_array(indp1)%x, bubble_points_array(indm1)%x, &
                & ux(indp1), ux(indm1), N)

			dy = diff_spline (bubble_points_array(indp1)%y, bubble_points_array(indm1)%y, &
                & uy(indp1), uy(indm1), N)

			gn = (bubble_points_array(j)%x - bubble_points_array(i)%x) / rs * dy - &
				& (bubble_points_array(j)%y - bubble_points_array(i)%y) / rs * dx
			
			s = gn / N
		else
			s = 0.0d0

			do k=1,N
				! Sum by k=1,N instead of fixed j
				if (k /= i) then	
					rs = (bubble_points_array(i)%x - bubble_points_array(k)%x) * &
		            	& (bubble_points_array(i)%x - bubble_points_array(k)%x) + &
			            & (bubble_points_array(i)%y - bubble_points_array(k)%y) * &
    			        & (bubble_points_array(i)%y - bubble_points_array(k)%y)

					indp1 = mod(k+1, N)
					if (indp1 == 0) indp1 = N
 	
					indm1 = mod(k-1, N)
					if (indm1 == 0) indm1 = N
	
					dx = diff_spline (bubble_points_array(indp1)%x, bubble_points_array(indm1)%x, &
        		        & ux(indp1), ux(indm1), N)

					dy = diff_spline (bubble_points_array(indp1)%y, bubble_points_array(indm1)%y, &
        		        & uy(indp1), uy(indm1), N)

					gn = (bubble_points_array(k)%x - bubble_points_array(i)%x) / rs * dy - &
						& (bubble_points_array(k)%y - bubble_points_array(i)%y) / rs * dx

					s = s + gn
				end if
			end do
			
			s = - s / N
		end if

	end function B_operator



	subroutine initial_bubble_state ()	
		integer(c_int) :: i

		do i=1,N
			bubble_points_array(i)%x = a * cos(2. * pi * idensity(i))
			bubble_points_array(i)%y = b * sin(2. * pi * idensity(i))

			! TODO
			bubble_points_array(i)%phi = cos(2. * pi * idensity(i))
			bubble_points_array(i)%K = 2.0d0 * pi
		end do
	end subroutine initial_bubble_state

	

	subroutine initial_domain ()
		type(gtkallocation), target:: alloc

		! c_loc - create c pointer for target	
		call gtk_widget_get_allocation(main_drawing_area,c_loc(alloc))
	    drawing_area_width = alloc%width
    	drawing_area_height = alloc%height
		
		x_max = x_min + drawing_area_width / ZOOM_VAR
		y_max = y_min + drawing_area_height / ZOOM_VAR	
	end subroutine initial_domain

	
	
	subroutine init_all ()
		call init_beta_array (N)
		call initial_domain ()
		call initial_bubble_state ()
	end subroutine init_all
	


	subroutine draw_current_bubble_state (cr)	
		type(c_ptr), intent(in) :: cr

		integer(c_int) :: i

		do i=1,N
			call get_coordinates_for_drawing (bubble_points_array(i)%x, bubble_points_array(i)%y, &
				& bubble_points_array(i)%dx, bubble_points_array(i)%dy)
		end do

		call cairo_set_source_rgb(cr, 1d0, 1d0, 1d0)
		call cairo_new_path(cr)
	
		call cairo_move_to (cr, bubble_points_array(1)%dx, bubble_points_array(1)%dy)
		
		do i=2,N
			call cairo_line_to (cr, bubble_points_array(i)%dx, &
    	        & bubble_points_array(i)%dy)
		end do

		call cairo_line_to (cr, bubble_points_array(1)%dx, bubble_points_array(1)%dy)
		
		call cairo_fill(cr)
	end subroutine draw_current_bubble_state



	subroutine draw_axis (cr)	
		type(c_ptr), intent(in) :: cr
		integer(c_int) :: i, i_min, i_max, i_d
		real(c_double) :: tmp_x, tmp_y, tmp_dx, tmp_dy
		real(c_double), dimension(2), target:: dashed = (/14.0d0, 6.0d0/)

		! Set color and width
		call cairo_set_source_rgb(cr, 0.93_c_double, 0.8_c_double, 0.53_c_double)
		call cairo_set_line_width(cr, 0.3d0)
		call cairo_set_dash(cr, c_loc(dashed), 2, 0.0d0)	

		! Horisontal axis

		i_min = ceiling(y_min)
		i_max = floor(y_max)
		i_d = i_max - i_min

		do i = 0, i_d
			call cairo_new_path(cr)

			tmp_x = x_min
			tmp_y = real(i_min + i, c_double)
			call get_coordinates_for_drawing (tmp_x, tmp_y, tmp_dx, tmp_dy)
			call cairo_move_to (cr, tmp_dx, tmp_dy)

			tmp_x = x_max
			tmp_y = real(i_min + i, c_double)
			call get_coordinates_for_drawing (tmp_x, tmp_y, tmp_dx, tmp_dy)
			call cairo_line_to (cr, tmp_dx, tmp_dy)
			call cairo_stroke (cr)
		end do

		! Vertical axis

		i_min = ceiling(x_min)
		i_max = floor(x_max)
		i_d = i_max - i_min

		do i = 0, i_d
			call cairo_new_path(cr)

			tmp_x = real(i_min + i, c_double)
			tmp_y = y_min
			call get_coordinates_for_drawing (tmp_x, tmp_y, tmp_dx, tmp_dy)
			call cairo_move_to (cr, tmp_dx, tmp_dy)

			tmp_x = real(i_min + i, c_double)
			tmp_y = y_max
			call get_coordinates_for_drawing (tmp_x, tmp_y, tmp_dx, tmp_dy)
			call cairo_line_to (cr, tmp_dx, tmp_dy)
			call cairo_stroke (cr)
		end do

	end subroutine draw_axis



	function draw_current_state (area, event, gdata) result(ret) bind(c)
		type(c_ptr), value :: area, event, gdata
		integer(c_int) :: ret
		type(gtkallocation), target:: alloc
		type(c_ptr) :: cr
		character(len=1024) :: text_full
	
	    ! Update drawing area size	
		call gtk_widget_get_allocation(main_drawing_area,c_loc(alloc))
	    drawing_area_width = alloc%width
    	drawing_area_height = alloc%height

		! Update text in the text area
		write (text_full, '(a10, i5, a1, a10, f10.3, a1, a10, f10.5)'), "Pause:", PAUSE_VAR, c_new_line, &
			& "Zoom:", ZOOM_VAR, c_new_line, "Time:", time
		
		call gtk_text_buffer_set_text (text_buffer, trim(text_full) // c_null_char, -1_c_int)

		cr = hl_gtk_drawing_area_cairo_new(area)
				
		! Fill background 
		call cairo_set_source_rgb(cr, 0.0_c_double, 0.086_c_double, 0.11_c_double)
		call cairo_rectangle(cr, 0._c_double, 0._c_double, real(width, c_double), real(height, c_double))
		call cairo_paint(cr)

		! Draw axis
		call draw_axis(cr)

		! Draw bubble
		call draw_current_bubble_state (cr)

		call hl_gtk_drawing_area_cairo_destroy(cr)
		call gtk_widget_queue_draw(area)

		if (PAUSE_VAR == 0) then
			call calc_next_state_test_1 ()
		end if 

		ret = TRUE
	end function draw_current_state



	function rescale_domain_for_drawing (area, event, gdata) result(ret) bind(c)
		type(c_ptr), value :: area, event, gdata
    	type(gdkeventscroll), pointer :: fevent
		integer(c_int) :: rescale_flag
		
		! Mouse coordinates
		real(c_double) :: mdx, mdy ! Drawing
		real(c_double) :: mx, my   ! Cartesian
		
		real(c_double) :: new_max, new_min

		if (.not. c_associated(event)) return ! shouldn't happen

		! Assign fortran pointer to the target of c pointer
    	call c_f_pointer(event, fevent)

		rescale_flag = 0

		if (fevent%direction == GDK_SCROLL_DOWN) then
			ZOOM_VAR = ZOOM_VAR / 1.1d0
			rescale_flag = 1
		end if

		if (fevent%direction == GDK_SCROLL_UP) then
			ZOOM_VAR = ZOOM_VAR * 1.1d0
			rescale_flag = 1
		end if

		if (rescale_flag == 1) then
			mdx = fevent%x
			mdy = fevent%y

			mx = x_min + (x_max - x_min) * mdx / drawing_area_width
			my = y_max - (y_max - y_min) * mdy / drawing_area_height
			
			
			new_max = mx + (x_max - mx) / (x_max - x_min) * drawing_area_width / ZOOM_VAR
			new_min = mx + (x_min - mx) / (x_max - x_min) * drawing_area_width / ZOOM_VAR
			x_max = new_max
			x_min = new_min

			new_max = my + (y_max - my) / (y_max - y_min) * drawing_area_height / ZOOM_VAR
			new_min = my + (y_min - my) / (y_max - y_min) * drawing_area_height / ZOOM_VAR
			y_max = new_max
			y_min = new_min

		end if
	
		ret = TRUE
	end function rescale_domain_for_drawing



	function button_pause_handler() result(ret) bind(c)
		integer(c_int) :: ret
		
		PAUSE_VAR = 1
		ret = TRUE
	end function button_pause_handler



	function button_continue_handler() result(ret) bind(c)
		integer(c_int) :: ret
		
		PAUSE_VAR = 0
		ret = TRUE
	end function button_continue_handler



	subroutine calc_next_state ()
		integer(c_int) :: ret, i, j, indp1, indm1, info
		integer(c_int) :: adams_index, adams_index_m1, adams_index_m2, adams_index_m3
		real(c_double) :: r, x, y, s, dlol 

		! Init all splines
		do i=1,N
			tmp_array(i) = bubble_points_array(i)%phi
		end do	
		call calc_periodic_spline_coefs (tmp_array, spline_phi)

		do i=1,N
			tmp_array(i) = bubble_points_array(i)%x
		end do	
		call calc_periodic_spline_coefs (tmp_array, spline_ux)

		do i=1,N
			tmp_array(i) = bubble_points_array(i)%y
		end do	
		call calc_periodic_spline_coefs (tmp_array, spline_uy)


		! Init matrix for equation B psi + 2pi psi = - A d(phi)/dn
		do i=1,N
			tmp_array(i) = 0.0d0
			do j=1,N
				indp1 = mod(j+1, N)
				if (indp1 == 0) indp1 = N
 
				indm1 = mod(j-1, N)
				if (indm1 == 0) indm1 = N

				tmp_matrix(i,j) = B_operator(i,j,spline_ux,spline_uy)
				
				tmp_array(i) = tmp_array(i) - A_operator(i,j) * &
					& diff_spline (bubble_points_array(indp1)%phi, &
					& bubble_points_array(indm1)%phi, spline_phi(indp1), &
					& spline_phi(indm1), N)
			end do
			tmp_matrix(i,i) = tmp_matrix(i,i) + 2.0d0 * pi
		end do

		! Solve it (Solution is stored in tmp_array)
		call dgesv (N, 1, tmp_matrix, N, tmp_array2, tmp_array, N, info)
		
		if (info /= 0) then
			stop "Cannot solve linear system" 
		end if

		do i=1,N
			 bubble_points_array(i)%psi = tmp_array(i)
		end do

		! Let v and u be normal and tangent speed respectively
		! Calculate V = v/l as V = d(psi)/d(dzeta) * 1/(l^2 * f)
		! We still have psi values in tmp_array

		call calc_periodic_spline_coefs (tmp_array, spline_psi)

		do i=1,N
			indp1 = mod(i+1, N)
			if (indp1 == 0) indp1 = N
 
			indm1 = mod(i-1, N)
			if (indm1 == 0) indm1 = N

			bubble_points_array(i)%V = diff_spline (tmp_array(indp1), &
				& tmp_array(indm1), spline_psi(indp1), spline_psi(indm1), N) / &
				& (l * l * density(i))
		end do

		! Now we want to calculate l'/l = int_0^1 (density K V) d(dzeta)
		! and U(dzeta_i) = U(0) - int_0^(dzeta_i) (density K V) d(dzeta) + 
		! 		+ idensity(dzeta_i) * l'/l
		! Note that U(0) = 0

		dlol = 0.0d0

		! Note that we don't need spline there
		do i=1,N
			tmp_array(i) = bubble_points_array(i)%K * density(i) * &
				& bubble_points_array(i)%V
			dlol = dlol + tmp_array(i) 
		end do		

		! l' / l
		dlol = dlol / N

		! Init spline for (density K V)
		call calc_periodic_spline_coefs (tmp_array, spline_densityKV)

		! Calculate U(dzeta_i)
		s = - int_spline_im1_i (tmp_array(1), tmp_array(N), &
			& spline_densityKV(1), spline_densityKV(N), N)
		bubble_points_array(1)%U = s + idensity(1) * dlol 
		
		do i=2,N
			s = s - int_spline_im1_i (tmp_array(i), tmp_array(i-1), &
				& spline_densityKV(i), spline_densityKV(i-1), N)
			bubble_points_array(i)%U = s + idensity(i) * dlol 
		end do

		! Now we can calculate new K from equation
		! K' = d (UK - dV/ddzeta / density) / ddzeta / density
		! Adams method

		adams_index = mod(adams_counter, 4)

		do i=1,N
			tmp_array(i) = bubble_points_array(i)%V
		end do

		! Init spline for V
		call calc_periodic_spline_coefs (tmp_array, spline_V)

		! Calculate (UK - dV/ddzeta / density)
		do i=1,N
			indp1 = mod(i+1, N)
			if (indp1 == 0) indp1 = N
 
			indm1 = mod(i-1, N)
			if (indm1 == 0) indm1 = N

			tmp_array2(i) = bubble_points_array(i)%U * bubble_points_array(i)%K - & 
				& diff_spline (tmp_array(indp1), tmp_array(indm1), spline_V(indp1), &
				& spline_V(indm1), N) / density(i)
		end do

		! Init spline for (UK - dV/ddzeta / density)
		call calc_periodic_spline_coefs (tmp_array2, spline_UKmdVof)

		! Calculate K_rhs

		! Increase counter of iterations first
		
		adams_counter = adams_counter + 1
		if (adams_counter == 9) adams_counter = 5

		! Calculate all indexes
		adams_index = mod(adams_counter, 4)
		if (adams_index == 0)  adams_index = 4

		adams_index_m1 = mod(adams_counter - 1, 4)
		if (adams_index_m1 == 0)  adams_index_m1 = 4

		adams_index_m2 = mod(adams_counter - 2, 4)
		if (adams_index_m2 == 0)  adams_index_m2 = 4

		adams_index_m3 = mod(adams_counter - 3, 4)
		if (adams_index_m3 == 0)  adams_index_m3 = 4

		do i=1,N
			indp1 = mod(i+1, N)
			if (indp1 == 0) indp1 = N
 
			indm1 = mod(i-1, N)
			if (indm1 == 0) indm1 = N

			bubble_points_array(i)%K_rhs(adams_index) = diff_spline(tmp_array2(indp1), &
				& tmp_array2(indm1), spline_UKmdVof(indp1), spline_UKmdVof(indm1), N) / &
				& density(i)
		end do

		! Adams method
		if (adams_counter == 1) then
			do i=1,N
				bubble_points_array(i)%K = bubble_points_array(i)%K + time_delta * &
					& bubble_points_array(i)%K_rhs(adams_index)
			end do
		else if (adams_counter == 2) then
			do i=1,N
				bubble_points_array(i)%K = bubble_points_array(i)%K + time_delta * &
					& (1.5d0 * bubble_points_array(i)%K_rhs(adams_index) - &
					& 0.5d0 * bubble_points_array(i)%K_rhs(adams_index_m1))
			end do
		else if (adams_counter == 3) then
			do i=1,N
				bubble_points_array(i)%K = bubble_points_array(i)%K + time_delta * &
					& (23.0d0 / 12.0d0  * bubble_points_array(i)%K_rhs(adams_index) - &
					& 4.0d0 / 3.0d0 * bubble_points_array(i)%K_rhs(adams_index_m1) + &
					& 5.0d0 / 12.0d0 * bubble_points_array(i)%K_rhs(adams_index_m2))
			end do
		else
			do i=1,N
				bubble_points_array(i)%K = bubble_points_array(i)%K + time_delta * &
					& (55.0d0 / 24.0d0 * bubble_points_array(i)%K_rhs(adams_index) - &
					& 59.0d0 / 24.0d0 * bubble_points_array(i)%K_rhs(adams_index_m1) + &
					& 37.0d0 / 24.0d0 * bubble_points_array(i)%K_rhs(adams_index_m2) - &
					& 3.0d0 / 8.0d0 * bubble_points_array(i)%K_rhs(adams_index_m3))

			end do
		end if
	
		! Now we can calculate teta as theta = pi/2 + int_0^dzeta (K density) d(dzeta)

		! Init spline for (K density)
		do i=1,N
			tmp_array(i) = bubble_points_array(i)%K * density(i)
		end do	
		
		call calc_periodic_spline_coefs (tmp_array, spline_densityK)

		s = pi * 0.5d0 + int_spline_im1_i (tmp_array(1), tmp_array(N), &
			& spline_densityK(1), spline_densityK(N), N)
		bubble_points_array(1)%theta = s

		do i=2,N
			s = s + int_spline_im1_i (tmp_array(i), tmp_array(i-1), &
	            & spline_densityK(i), spline_densityK(i-1), N)
			bubble_points_array(i)%theta = s
		end do

		! Now we can calculate x as 
		! x = x_0 + l * int_0^dzeta (cos(theta) * density) d(dzeta)
		
		! Init spline for (cos(theta) * density)
		do i=1,N
			tmp_array(i) = cos(bubble_points_array(i)%theta) * density(i)
		end do	
		
		call calc_periodic_spline_coefs (tmp_array, spline_densityCosTheta)

		! Update x_0 as x_0(t+dt) = x_0(t) + dt * v_x = x_0(t) + dt * V * l
		
		s = (bubble_points_array(N)%x + time_delta * bubble_points_array(N)%V * l) + &
			& int_spline_im1_i (tmp_array(1), tmp_array(N), &
			& spline_densityCosTheta(1), spline_densityCosTheta(N), N) * l
		bubble_points_array(1)%x = s

		do i=2,N
			s = s + int_spline_im1_i (tmp_array(i), tmp_array(i-1), &
	            & spline_densityCosTheta(i), spline_densityCosTheta(i-1), N) * l
			bubble_points_array(i)%x = s
		end do

		! Now we can calculate y as 
		! y = y_0 + l * int_0^dzeta (sin(theta) * density) d(dzeta)
		
		! Init spline for (cos(theta) * density)
		do i=1,N
			tmp_array(i) = sin(bubble_points_array(i)%theta) * density(i)
		end do	
		
		call calc_periodic_spline_coefs (tmp_array, spline_densitySinTheta)

		! y_0 = 0
		
		s = int_spline_im1_i (tmp_array(1), tmp_array(N), &
			& spline_densitySinTheta(1), spline_densitySinTheta(N), N) * l
		bubble_points_array(1)%y = s

		do i=2,N
			s = s + int_spline_im1_i (tmp_array(i), tmp_array(i-1), &
	            & spline_densitySinTheta(i), spline_densitySinTheta(i-1), N) * l
			bubble_points_array(i)%y = s
		end do

		! Update phi (Bernoulli integral)
		! We already have spline coefficients for phi
		! And we have adams indexes

		do i=1,N
			indp1 = mod(i+1, N)
			if (indp1 == 0) indp1 = N
 
			indm1 = mod(i-1, N)
			if (indm1 == 0) indm1 = N

			s = diff_spline(bubble_points_array(indp1)%phi, &
                & bubble_points_array(indm1)%phi, spline_phi(indp1), &
				& spline_phi(indm1), N)

			bubble_points_array(i)%phi_rhs(adams_index) = &
				& (bubble_points_array(i)%U * s / density(i) + s * s * 0.5d0 / (l * l * &
				& density(i) * density(i)) - 0.5d0 * bubble_points_array(i)%V * &
				& bubble_points_array(i)%V * l * l + bernoulli_const)
		end do
	
		! Adams method
		if (adams_counter == 1) then
			do i=1,N
				bubble_points_array(i)%phi = bubble_points_array(i)%phi + time_delta * &
					& bubble_points_array(i)%phi_rhs(adams_index)
			end do
		else if (adams_counter == 2) then
			do i=1,N
				bubble_points_array(i)%phi = bubble_points_array(i)%phi + time_delta * &
					& (1.5d0 * bubble_points_array(i)%phi_rhs(adams_index) - &
					& 0.5d0 * bubble_points_array(i)%phi_rhs(adams_index_m1))
			end do
		else if (adams_counter == 3) then
			do i=1,N
				bubble_points_array(i)%phi = bubble_points_array(i)%phi + time_delta * &
					& (23.0d0 / 12.0d0  * bubble_points_array(i)%phi_rhs(adams_index) - &
					& 4.0d0 / 3.0d0 * bubble_points_array(i)%phi_rhs(adams_index_m1) + &
					& 5.0d0 / 12.0d0 * bubble_points_array(i)%phi_rhs(adams_index_m2))
			end do
		else
			do i=1,N
				bubble_points_array(i)%phi = bubble_points_array(i)%phi + time_delta * &
					& (55.0d0 / 24.0d0 * bubble_points_array(i)%phi_rhs(adams_index) - &
					& 59.0d0 / 24.0d0 * bubble_points_array(i)%phi_rhs(adams_index_m1) + &
					& 37.0d0 / 24.0d0 * bubble_points_array(i)%phi_rhs(adams_index_m2) - &
					& 3.0d0 / 8.0d0 * bubble_points_array(i)%phi_rhs(adams_index_m3))
			end do
		end if
	
		! Update l as l(t+dt) = l(t) + dt * (l'/l)(t) * l(t)
		! And we have adams indexes
		
		l_rhs(adams_index) = l * dlol	

		! Adams method
		if (adams_counter == 1) then
			l = l + time_delta * l_rhs(adams_index)
		else if (adams_counter == 2) then
			l = l + time_delta * &
				& (1.5d0 * l_rhs(adams_index) - &
				& 0.5d0 * l_rhs(adams_index_m1))
		else if (adams_counter == 3) then
			l = l + time_delta * &
				& (23.0d0 / 12.0d0  * l_rhs(adams_index) - &
				& 4.0d0 / 3.0d0 * l_rhs(adams_index_m1) + &
				& 5.0d0 / 12.0d0 * l_rhs(adams_index_m2))
		else
			l = l + time_delta * &
				& (55.0d0 / 24.0d0 * l_rhs(adams_index) - &
				& 59.0d0 / 24.0d0 * l_rhs(adams_index_m1) + &
				& 37.0d0 / 24.0d0 * l_rhs(adams_index_m2) - &
				& 3.0d0 / 8.0d0 * l_rhs(adams_index_m3))
		end if

		! Update time
		time = time + time_delta

		do i=N/2 - 5, N/2 + 5
			print *, "K", i, bubble_points_array(i)%K
			print *, "V", i, bubble_points_array(i)%V
			print *, "U", i, bubble_points_array(i)%U
			print *, "psi", i, bubble_points_array(i)%psi
			print *, "phi", i, bubble_points_array(i)%phi
		end do

	end subroutine calc_next_state


	subroutine calc_next_state_test_1
		integer(c_int) :: ret, i, j, indp1, indm1, info
		integer(c_int) :: adams_index, adams_index_m1, adams_index_m2, adams_index_m3
		real(c_double) :: r, x, y, s, dlol 

		! Init all splines
		do i=1,N
			tmp_array(i) = bubble_points_array(i)%x
		end do	
		call calc_periodic_spline_coefs (tmp_array, spline_ux)

		do i=1,N
			tmp_array(i) = bubble_points_array(i)%y
		end do	
		call calc_periodic_spline_coefs (tmp_array, spline_uy)

		do i=1,N
			indp1 = mod(i+1, N)
			if (indp1 == 0) indp1 = N
 
			indm1 = mod(i-1, N)
			if (indm1 == 0) indm1 = N

			bubble_points_array(i)%nx = diff_spline (bubble_points_array(indp1)%y, &
				& bubble_points_array(indm1)%y, spline_uy(indp1), spline_uy(indm1), N) / &
				& (l * density(i))

			bubble_points_array(i)%ny = - diff_spline (bubble_points_array(indp1)%x, &
				& bubble_points_array(indm1)%x, spline_ux(indp1), spline_ux(indm1), N) / &
				& (l * density(i))
		end do

		! Init V=v/l from v_x = x, v_y = 0 and (nx, ny)
		do i=1,N
			bubble_points_array(i)%V = bubble_points_array(i)%nx * &
				& bubble_points_array(i)%x / l
		end do

! COPY calc_next_state (except phi)

		! Now we want to calculate l'/l = int_0^1 (density K V) d(dzeta)
		! and U(dzeta_i) = U(0) - int_0^(dzeta_i) (density K V) d(dzeta) + 
		! 		+ idensity(dzeta_i) * l'/l
		! Note that U(0) = 0

		dlol = 0.0d0

		! Note that we don't need spline there
		do i=1,N
			tmp_array(i) = bubble_points_array(i)%K * density(i) * &
				& bubble_points_array(i)%V
			dlol = dlol + tmp_array(i) 
		end do		

		! l' / l
		dlol = dlol / N

		! Init spline for (density K V)
		call calc_periodic_spline_coefs (tmp_array, spline_densityKV)

		! Calculate U(dzeta_i)
		s = - int_spline_im1_i (tmp_array(1), tmp_array(N), &
			& spline_densityKV(1), spline_densityKV(N), N)
		bubble_points_array(1)%U = s + idensity(1) * dlol 
		
		do i=2,N
			s = s - int_spline_im1_i (tmp_array(i), tmp_array(i-1), &
				& spline_densityKV(i), spline_densityKV(i-1), N)
			bubble_points_array(i)%U = s + idensity(i) * dlol 
		end do

		! Now we can calculate new K from equation
		! K' = d (UK - dV/ddzeta / density) / ddzeta / density
		! Adams method

		adams_index = mod(adams_counter, 4)

		do i=1,N
			tmp_array(i) = bubble_points_array(i)%V
		end do

		! Init spline for V
		call calc_periodic_spline_coefs (tmp_array, spline_V)

		! Calculate (UK - dV/ddzeta / density)
		do i=1,N
			indp1 = mod(i+1, N)
			if (indp1 == 0) indp1 = N
 
			indm1 = mod(i-1, N)
			if (indm1 == 0) indm1 = N

			tmp_array2(i) = bubble_points_array(i)%U * bubble_points_array(i)%K - & 
				& diff_spline (tmp_array(indp1), tmp_array(indm1), spline_V(indp1), &
				& spline_V(indm1), N) / density(i)
		end do

		! Init spline for (UK - dV/ddzeta / density)
		call calc_periodic_spline_coefs (tmp_array2, spline_UKmdVof)

		! Calculate K_rhs

		! Increase counter of iterations first
		
		adams_counter = adams_counter + 1
		if (adams_counter == 9) adams_counter = 5

		! Calculate all indexes
		adams_index = mod(adams_counter, 4)
		if (adams_index == 0)  adams_index = 4

		adams_index_m1 = mod(adams_counter - 1, 4)
		if (adams_index_m1 == 0)  adams_index_m1 = 4

		adams_index_m2 = mod(adams_counter - 2, 4)
		if (adams_index_m2 == 0)  adams_index_m2 = 4

		adams_index_m3 = mod(adams_counter - 3, 4)
		if (adams_index_m3 == 0)  adams_index_m3 = 4

		do i=1,N
			indp1 = mod(i+1, N)
			if (indp1 == 0) indp1 = N
 
			indm1 = mod(i-1, N)
			if (indm1 == 0) indm1 = N

			bubble_points_array(i)%K_rhs(adams_index) = diff_spline(tmp_array2(indp1), &
				& tmp_array2(indm1), spline_UKmdVof(indp1), spline_UKmdVof(indm1), N) / &
				& density(i)
		end do

		! Adams method
		if (adams_counter == 1) then
			do i=1,N
				bubble_points_array(i)%K = bubble_points_array(i)%K + time_delta * &
					& bubble_points_array(i)%K_rhs(adams_index)
			end do
		else if (adams_counter == 2) then
			do i=1,N
				bubble_points_array(i)%K = bubble_points_array(i)%K + time_delta * &
					& (1.5d0 * bubble_points_array(i)%K_rhs(adams_index) - &
					& 0.5d0 * bubble_points_array(i)%K_rhs(adams_index_m1))
			end do
		else if (adams_counter == 3) then
			do i=1,N
				bubble_points_array(i)%K = bubble_points_array(i)%K + time_delta * &
					& (23.0d0 / 12.0d0  * bubble_points_array(i)%K_rhs(adams_index) - &
					& 4.0d0 / 3.0d0 * bubble_points_array(i)%K_rhs(adams_index_m1) + &
					& 5.0d0 / 12.0d0 * bubble_points_array(i)%K_rhs(adams_index_m2))
			end do
		else
			do i=1,N
				bubble_points_array(i)%K = bubble_points_array(i)%K + time_delta * &
					& (55.0d0 / 24.0d0 * bubble_points_array(i)%K_rhs(adams_index) - &
					& 59.0d0 / 24.0d0 * bubble_points_array(i)%K_rhs(adams_index_m1) + &
					& 37.0d0 / 24.0d0 * bubble_points_array(i)%K_rhs(adams_index_m2) - &
					& 3.0d0 / 8.0d0 * bubble_points_array(i)%K_rhs(adams_index_m3))

			end do
		end if
	
		! Now we can calculate teta as theta = pi/2 + int_0^dzeta (K density) d(dzeta)

		! Init spline for (K density)
		do i=1,N
			tmp_array(i) = bubble_points_array(i)%K * density(i)
		end do	
		
		call calc_periodic_spline_coefs (tmp_array, spline_densityK)

		s = pi * 0.5d0 + int_spline_im1_i (tmp_array(1), tmp_array(N), &
			& spline_densityK(1), spline_densityK(N), N)
		bubble_points_array(1)%theta = s

		do i=2,N
			s = s + int_spline_im1_i (tmp_array(i), tmp_array(i-1), &
	            & spline_densityK(i), spline_densityK(i-1), N)
			bubble_points_array(i)%theta = s
		end do

		! Now we can calculate x as 
		! x = x_0 + l * int_0^dzeta (cos(theta) * density) d(dzeta)
		
		! Init spline for (cos(theta) * density)
		do i=1,N
			tmp_array(i) = cos(bubble_points_array(i)%theta) * density(i)
		end do	
		
		call calc_periodic_spline_coefs (tmp_array, spline_densityCosTheta)

		! Update x_0 as x_0(t+dt) = x_0(t) + dt * v_x = x_0(t) + dt * V * l
		
		s = (bubble_points_array(N)%x + time_delta * bubble_points_array(N)%V * l) + &
			& int_spline_im1_i (tmp_array(1), tmp_array(N), &
			& spline_densityCosTheta(1), spline_densityCosTheta(N), N) * l
		bubble_points_array(1)%x = s

		do i=2,N
			s = s + int_spline_im1_i (tmp_array(i), tmp_array(i-1), &
	            & spline_densityCosTheta(i), spline_densityCosTheta(i-1), N) * l
			bubble_points_array(i)%x = s
		end do

		! Now we can calculate y as 
		! y = y_0 + l * int_0^dzeta (sin(theta) * density) d(dzeta)
		
		! Init spline for (cos(theta) * density)
		do i=1,N
			tmp_array(i) = sin(bubble_points_array(i)%theta) * density(i)
		end do	
		
		call calc_periodic_spline_coefs (tmp_array, spline_densitySinTheta)

		! y_0 = 0
		
		s = int_spline_im1_i (tmp_array(1), tmp_array(N), &
			& spline_densitySinTheta(1), spline_densitySinTheta(N), N) * l
		bubble_points_array(1)%y = s

		do i=2,N
			s = s + int_spline_im1_i (tmp_array(i), tmp_array(i-1), &
	            & spline_densitySinTheta(i), spline_densitySinTheta(i-1), N) * l
			bubble_points_array(i)%y = s
		end do

		! Update l as l(t+dt) = l(t) + dt * (l'/l)(t) * l(t)
		! And we have adams indexes
		
		l_rhs(adams_index) = l * dlol	

		! Adams method
		if (adams_counter == 1) then
			l = l + time_delta * l_rhs(adams_index)
		else if (adams_counter == 2) then
			l = l + time_delta * &
				& (1.5d0 * l_rhs(adams_index) - &
				& 0.5d0 * l_rhs(adams_index_m1))
		else if (adams_counter == 3) then
			l = l + time_delta * &
				& (23.0d0 / 12.0d0  * l_rhs(adams_index) - &
				& 4.0d0 / 3.0d0 * l_rhs(adams_index_m1) + &
				& 5.0d0 / 12.0d0 * l_rhs(adams_index_m2))
		else
			l = l + time_delta * &
				& (55.0d0 / 24.0d0 * l_rhs(adams_index) - &
				& 59.0d0 / 24.0d0 * l_rhs(adams_index_m1) + &
				& 37.0d0 / 24.0d0 * l_rhs(adams_index_m2) - &
				& 3.0d0 / 8.0d0 * l_rhs(adams_index_m3))
		end if

		! Update time
		time = time + time_delta

	end subroutine calc_next_state_test_1


	! Transform double point (x, y) into (dx, dy) for drawing area
	subroutine get_coordinates_for_drawing (x, y, dx, dy)
		real(kind=c_double), intent(in) :: x, y
		real(kind=c_double), intent(out) :: dx, dy

		dx = (x - x_min) / (x_max - x_min) * drawing_area_width
		dy = (1.0d0 - (y - y_min) / (y_max - y_min)) * drawing_area_height
			
	end subroutine get_coordinates_for_drawing

end module handlers


