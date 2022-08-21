set(TESTS
    test_delaunay_2d_1
    test_new_works_1
    test_find_works_1
    test_large_tris_1
    test_interp_works_1
)

foreach(var ${TESTS})
    add_executable(${var} "${var}.cpp")
    add_test(NAME ${var} COMMAND ${var})
    target_link_libraries(${var} PUBLIC ${TRILIB_NAME})
endforeach(var)
