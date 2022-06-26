# Cpp implementation of GridSearch

## Concept

A grid is composed of containers. The container (aka "bin") is uniquely identified by its
"left-most" corner. Any coordinated in 2D or 3D can be quickly compared with the "left-most"
corner; therefore we can find coordinates efficiently.

To avoid repetition (save space and be more efficient) we store all coordinates in a nested map.
The outer map is a set of containers whereas the inner map is a set of items. The outer map
maps `index-to-container`. The inner map maps `id-to-item`.

To use the grid, we fist insert items using the `insert` command. Then, we can find any
item by using the `find` command.

This idea is illustrated below:

```cpp
auto options = GridSearchOptions::make_new();
auto grid = GridSearch::make_new(options);
vector<double> x = {0.5, 0.5};
grid->insert(100, x);
auto id = grid->find(x);
// will return id = 100
```

## Full Example

Output of [example.cpp](https://github.com/cpmech/cpp-grid-search/blob/main/example.cpp):

```text
Running example
number of non-empty containers = 3
container # 0: items = [100:(0.5,0.5), 200:(0.8,0.5)]
container # 1: items = [300:(1,0.5)]
container # 8: items = [400:(2.9,2.9)]
found id=100 near (0.5,0.5)
found id=300 near (1,0.5)
found id=400 near (2.9,2.9)
found id=-1 near (3,3)
Done
```
