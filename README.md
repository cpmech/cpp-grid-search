# Cpp implementation of GridSearch

## Concept

A grid is composed of containers. The container (aka "bin") is uniquely identified by its "left-most" corner (aka **key**). Any coordinate in 2D or 3D can be quickly compared with the "left-most" corner; therefore, we can find coordinates efficiently.

The main idea revolves around the following expression that calculates the index of a container where the coordinates `x` fall in (Eq. 8 of the Reference):

```text
ratio[i] = truncate((x[i] - xmin[i]) / Δx[i])
key = ratio[0] + ratio[1] * ndiv[0] 
```

Below is an illustration of a grid with two triangles:

![example_grid_two_triangles](https://github.com/cpmech/cpp-grid-search/raw/main/data/figures/example_grid_two_triangles.svg)

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

**Update.** Nonetheless, a performance improvement is available now where the large triangles are stored in a separate list. If 20% or fewer triangles are "too large," these large triangles are held in a separate list which will not affect the grid sizing. Otherwise, if more than 20% (GS_LARGE_CELLS_MAX_COUNT_PCT) triangles are "large," the "standard" algorithm takes place, i.e., no separate list is used, and all triangles affect the grid sizing. The definition of a "large" triangle is as follows: A large triangle has the largest dimension of its bounding box greater than or equal to 0.75 (GS_LARGE_CELLS_CUTOFF) times the size of the largest dimension of the largest triangle among all triangles.

The figures below show before and after the update. The darker yellow triangles are the "large" triangles in the mesh.

**Before**

![before](https://github.com/cpmech/cpp-grid-search/raw/main/data/figures/test_grid_search_tri_find_works_old.svg)

**After**

![after](https://github.com/cpmech/cpp-grid-search/raw/main/data/figures/test_grid_search_tri_find_works.svg)

## Usage

You only need to copy-n-paste the file `grid_search.h` into your project and use it as [shown in the example](https://github.com/cpmech/cpp-grid-search/blob/main/example_triangles.cpp).

First, we allocate the grid:

```c++
auto grid = GridSearch::make_new(coordinates, triangles);
```

where `coordinates` is a table with the x-y-temperature data; i.e., the third column contains the values to be interpolated. 

Second, we can interpolate values as follows:

```c++
vector<double> x = {0.5, 0.5};
auto temperature = grid->find_triangle_and_interpolate(x, coordinates, triangles);
```

where we must pass the same `coordinates` and `triangles` used in `new`.

## Full Example

Given the following mesh:

![example_grid_search](https://github.com/cpmech/cpp-grid-search/raw/main/data/figures/example_triangles.png)

Interpolate data using [example_grid_search.cpp](https://github.com/cpmech/cpp-grid-search/blob/main/examples/example_grid_search.cpp):

Output:

```text
x = {0.5, 0.5}
temperature = 0.785973
```

## Delaunay triangulation

In 2D, Delaunay triangulation is performed using the fantastic [Triangle](https://www.cs.cmu.edu/~quake/triangle.html). In 3D, the Delaunay tetrahedralization is performed with TetGen (and old version).

Given a "cloud" of points:

```c++
// generate Delaunay triangulation/tetrahedralization
auto triangles = delaunay_2d(cloud_2d, false);
auto tetrahedra = delaunay_3d(cloud_3d, false);
```

Full example [example_delaunay_2d.cpp](https://github.com/cpmech/cpp-grid-search/blob/main/examples/example_delaunay_2d.cpp)

Full example [example_delaunay_2d.cpp](https://github.com/cpmech/cpp-grid-search/blob/main/examples/example_delaunay_3d.cpp)

The code should generate a mesh like the one below:

![example_delaunay](https://github.com/cpmech/cpp-grid-search/raw/main/data/figures/doc_triangle_delaunay_1.svg)

## Reference

* Durand, Farias, and Pedroso (2015) Computing intersections between
  non-compatible curves and finite elements, Computational Mechanics;
  [DOI=10.1007/s00466-015-1181-y](https://link.springer.com/article/10.1007/s00466-015-1181-y)
