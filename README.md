# Cpp implementation of GridSearch

## Concept

A grid is composed of containers. The container (aka "bin") is uniquely identified by its "left-most" corner (aka **key**). Any coordinate in 2D or 3D can be quickly compared with the "left-most" corner; therefore, we can find coordinates efficiently.

To avoid repetition, save space, and be more efficient, we store all coordinates in a nested data structure. The outer structure is a map of containers that maps the container key to the inner data structure. The inner structure is a set that holds unique identifiers.

## Full Example

Output of [example_triangles.cpp](https://github.com/cpmech/cpp-grid-search/blob/main/example_triangles.cpp):

```text
Running example
number of non-empty containers = 15
container # 33: items = [1, 5, 11]
container # 34: items = [7]
container # 35: items = [5, 7, 11]
container # 43: items = [0, 1, 2, 3, 5, 11]
container # 44: items = [2, 3, 6, 7]
container # 45: items = [0, 5, 6, 7, 8, 9, 11]
container # 46: items = [8]
container # 53: items = [0, 3, 4, 10]
container # 54: items = [3, 6]
container # 55: items = [0, 4, 6, 9]
container # 56: items = [10]
container # 63: items = [2, 4, 10]
container # 64: items = [2]
container # 65: items = [4, 8]
container # 66: items = [8, 10]
found triangle with id = 6
Done
```
