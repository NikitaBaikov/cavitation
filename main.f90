! Main program

program bubble
	use handlers
	implicit none
	
	integer(kind=c_int) :: timeid, ret

	call gtk_init ()

	! Properties of the main window :
	main_window = hl_gtk_window_new (delete_event = c_funloc(quit))

	call gtk_window_set_default_size(main_window, width, height)
	call gtk_window_set_title(main_window, "Bubble GUI"//c_null_char)
	

	! Set border within window clear of any elements
	call gtk_container_set_border_width (main_window, 10_c_int)

	! Divide working domain within window into table with 15 rows and 5 columns
	! It works as a coordinate system for our GUI
	table = gtk_table_new (15_c_int, 20_c_int, TRUE)
	call gtk_table_set_row_spacings (table, 5_c_int)
	call gtk_table_set_col_spacings (table, 5_c_int)
	call gtk_container_add (main_window, table)

	! Drawing area
  	main_drawing_area = hl_gtk_drawing_area_new (scroll_event = c_funloc(rescale_domain_for_drawing))
	
	! Pause button
	button_pause = gtk_button_new_with_label ("Pause"//c_null_char)
	call g_signal_connect (button_pause, "clicked"//c_null_char, c_funloc(button_pause_handler))

	! Continue button
	button_continue = gtk_button_new_with_label ("Continue"//c_null_char)
	call g_signal_connect (button_continue, "clicked"//c_null_char, c_funloc(button_continue_handler))

	! Text area
	view = gtk_text_view_new ()
	text_buffer = gtk_text_view_get_buffer (view)	

	! Attach all widgets to table
	call gtk_table_attach_defaults(table, main_drawing_area, 0_c_int, 15_c_int, 0_c_int, 15_c_int)  
	call gtk_table_attach_defaults(table, button_pause, 15_c_int, 20_c_int, 0_c_int, 1_c_int)
	call gtk_table_attach_defaults(table, button_continue, 15_c_int, 20_c_int, 1_c_int, 2_c_int)
	call gtk_table_attach_defaults(table, view, 15_c_int, 20_c_int, 2_c_int, 15_c_int)

   	! Show 
	call gtk_widget_show_all (main_window)

	! Init all initial state
	call init_all ()
	ret = draw_current_state (main_drawing_area, c_null_ptr, c_null_ptr)

 	! Redraw every ... mseconds
	timeid = g_timeout_add(100_c_int, c_funloc(draw_current_state), main_drawing_area)

	! GTK main	
	call gtk_main ()
end program bubble

