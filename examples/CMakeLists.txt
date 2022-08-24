set(EXAMPLES
    example_delaunay_2d
    example_delaunay_3d
    example_grid_search_2d
    example_grid_search_3d
)

foreach(var ${EXAMPLES})
    add_executable(${var} "${var}.cpp")
    target_link_libraries(${var} PUBLIC ${TRILIB_NAME} ${TETLIB_NAME})
endforeach(var)
