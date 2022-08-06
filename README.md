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

![example_triangles](https://github.com/cpmech/cpp-grid-search/raw/main/example_triangles.png)

Interpolate data using [example_triangles.cpp](https://github.com/cpmech/cpp-grid-search/blob/main/example_triangles.cpp):

Output:

```text
x = {0.5, 0.5}
temperature = 0.785973
```

## Reference

* Durand, Farias, and Pedroso (2015) Computing intersections between
  non-compatible curves and finite elements, Computational Mechanics;
  [DOI=10.1007/s00466-015-1181-y](https://link.springer.com/article/10.1007/s00466-015-1181-y)
