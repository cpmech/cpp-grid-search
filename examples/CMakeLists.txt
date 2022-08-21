set(EXAMPLES
    example_delaunay_2d
    example_grid_search
)

foreach(var ${EXAMPLES})
    add_executable(${var} "${var}.cpp")
    target_link_libraries(${var} PUBLIC ${TRILIB_NAME})
endforeach(var)
