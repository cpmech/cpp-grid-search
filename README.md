# Cpp implementation of GridSearch

## Concept

A grid is composed of containers. The container (aka "bin") is uniquely identified by its "left-most" corner (aka **key**). Any coordinate in 2D or 3D can be quickly compared with the "left-most" corner; therefore, we can find coordinates efficiently.

The main idea revolves around the following expression that calculates the index of a container where the coordinates `x` fall in (Eq. 8 of the Reference):

```text
ratio[i] = truncate((x[i] - xmin[i]) / Δx[i])
key = ratio[0] + ratio[1] * ndiv[0] 
```

Below is an illustration of a grid with two triangles:

![example_grid_two_triangles](https://github.com/cpmech/cpp-grid-search/raw/main/example_grid_two_triangles.svg)

Only the containers touched by the bounding box of the triangle are saved. Thus, we use a `map` from container key to a "set" of triangles belonging to this container. Therefore, all data is stored in the following data structure;

```c++
typedef std::set<size_t> Container_t;
typedef std::map<size_t, Container_t> Containers_t;
```

This data structure avoids repetition, saves space, and is somewhat efficient.

## Limitations

The `GridSearch` only works if the minimum container size (edge/side length) is greater than the maximum dimension of the largest triangle. Therefore, if one triangle is much larger that the other ones, the algorithm won't perform as well as it could. On the other hand, if the triangles have similar "sizes," then the search should be fast.

In summary, we need to make sure that:

* The container's `side_length` must be greater than the maximum dimension of the largest triangle

## Usage

You only need to copy-n-paste the file `grid_search.h` into your project and use it as [shown in the example](https://github.com/cpmech/cpp-grid-search/blob/main/example_triangles.cpp).

First, we allocate the grid:

```c++
auto grid = GridSearch::make_new(tris.size(), get_x);
```

where `get_x` is a lambda function that returns an access to the coordinates of point `m` of triangle `t`. This could be read from two lists of coordinates and connectivity or directly from the x-y values tabled in a nested vector data structure (see, e.g., `tris` in the example). The `get_x` lambda is such that:

```c++
auto get_x = [&tris](size_t t, size_t m) {
    return &tris[t][m];
};
```

Note that `get_x` returns a pointer to the innermost vector in `vector<vector<vector<double>>> tris`.

Second, we can search for cells containing points using the fragment below:

```c++
vector<double> x = {0.5, 0.5};
auto id = grid->find_cell(x, is_in_cell);
```

where `is_in_cell` is another lambda function that determines if the point is inside the geometric shape (triangle). For example:

```c++
auto is_in_cell = [&tris](size_t t, vector<double> const *x) {
    return in_triangle(tris[t][0], tris[t][1], tris[t][2], (*x));
};
```

where `in_triangle` implements the algorithm to detect if a point is inside the triangle, given its vertices.

## Full Example

Given the following mesh:

![example_triangles](https://github.com/cpmech/cpp-grid-search/raw/main/example_triangles.svg)

Find points using [example_triangles.cpp](https://github.com/cpmech/cpp-grid-search/blob/main/example_triangles.cpp):

Output:

```text
GridSearch
==========
number of non-empty containers = 4
container # 0: cells = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
container # 1: cells = [8, 10]
container # 2: cells = [2, 4, 8, 10]
container # 3: cells = [8, 10]

x = {0.5, 0.5}
found triangle with id = 6

x = {0.4, 0.2}
found triangle with id = 11
```

## Reference

* Durand, Farias, and Pedroso (2015) Computing intersections between
  non-compatible curves and finite elements, Computational Mechanics;
  [DOI=10.1007/s00466-015-1181-y](https://link.springer.com/article/10.1007/s00466-015-1181-y)
