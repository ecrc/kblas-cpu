# Defines the source files of the tests. Each file generates 1 test
file(GLOB tests_files "*.c")

include_directories( ${PROJECT_SOURCE_DIR}/include )


# Uses RUNPATH instead of RPATH
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--enable-new-dtags")

foreach( test_src ${tests_files} )
	get_filename_component( test_exe ${test_src}  NAME_WE )
	add_executable( ${test_exe} ${test_src} )
	# Files must be compiled with c++ compiler, regarding thier extension
	set_source_files_properties( ${test_src} PROPERTIES LANGUAGE CXX )
	target_link_libraries( ${test_exe} kblas ${CBLAS_LIBRARIES} )
	#target_link_libraries( ${test_exe} kblas -mkl=parallel  )
	install( TARGETS ${test_exe} DESTINATION lib/testing )

	# Add basic tests. To run them type: 'make test' in build directory
	add_test( NAME "${test_exe}:1" COMMAND ${test_exe} -N 2048:4096:1024 --niter 3 )

endforeach()



