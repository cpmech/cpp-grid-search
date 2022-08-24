set(TESTS
    test_delaunay_2d_1
    test_delaunay_3d_1
    test_grid_search_2d_1
    test_grid_search_2d_2
    test_tritet_coords
)

foreach(var ${TESTS})
    add_executable(${var} "${var}.cpp")
    add_test(NAME ${var} COMMAND ${var})
    target_link_libraries(${var} PUBLIC ${TRILIB_NAME} ${TETLIB_NAME})
endforeach(var)
