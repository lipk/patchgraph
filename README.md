# Prototype patch graph implementation in C++

[![License](https://img.shields.io/github/license/mashape/apistatus.svg)]()

A patch graph is a dynamically refinable rectangular grid for performing
scientific computations. Its main feature is the ability to maintain large
continuous patches.

This repository contains a pre-alpha, buggy, ugly and incomplete patch graph
implementation in C++. It's not usable for anything other than playing around.

## Dependencies

C++11 compatible compiler, CMake and SDL2. SDL2 is used for visualization only
and will be removed when early development is over.

## Usage

First you'll need an *upsample* and a *downsample* function. The downsampler
gets a bunch of cells from a higher resolution grid and combines them into one
value which is put in a lower resolution cell. Here is an example that simply
sums the input:

```cpp
auto downsample = [](const DataReader<double>& reader,
                     u32 fromX,
                     u32 fromY,
                     u32 toX,
                     u32 toY) -> double {
    double sum = 0;
    for (u32 y = fromY; y < toY; ++y) {
        for (u32 x = fromX; x < toX; ++x) {
            sum += reader.read(x, y);
        }
    }
    return sum;
};
```

`fromX`, `fromY`, `toX` and `toY` define the interval (left-inclusive) to be
sampled. In case of downsampling, this is always a square of some size. `reader`
can be used to obtain the values at a given x, y index. Technically you can read
cells outside the given interval, but it's not a good idea.

Upsampling is the reverse direction, converting a single low resolution cell
into several high resolution ones. This is a bit trickier, because the target
area is not always a square. There are three possibilities:

* When refining a cell, the target is an n-by-n square
* When synchronizing the boundaries of patches of different refinement levels,
  the target is a column (for vertical edges) or a row (for horizontal edges) of
  cells
* When synchronizing the corners between of patches, the target is a single cell

With that in mind, our example upsampler:

```cpp
auto upsample = [](DataWriter<double>& writer,
                   double value,
                   u32 fromX,
                   u32 fromY,
                   u32 toX,
                   u32 toY,
                   u8 focusDiff) {
    const double sampled = value / (focusDiff * focusDiff);
    for (u32 y = fromY; y < toY; ++y) {
        for (u32 x = fromX; x < toX; ++x) {
            writer.write(x, y, sampled);
        }
    }
};
```

`writer`, `fromX`, `fromY`, `toX` and `toY` have analogous meaning to their
counterparts in the downsampler. `value` is the input value that should be
refined. `focusDiff` is the ratio of the length of an edge of a cell in the
higher and the lower resolution grid.

Now we can create a graph:

```cpp
auto graph = createPatchGraph<double>(
    7U /*width*/, 7U /*height*/, std::move(downsample), std::move(upsample));
```

Initially every cell is set to 0. We can modify them with `write`:

```cpp
for (size_t y = 0; y < 7U; ++y) {
    for (size_t x = 0; x < 7U; ++x) {
        graph.write(x, y, 1);
    }
}
```

It's worth noting that no matter how the initial grid gets split up and refined,
`write` (and `read`) will always see it as a 7x7 array (or whatever size you
initially define).

Time to do some calculations! We can use `apply` for that:

```cpp
graph.apply(1 /*iterations*/,
            [](const DataReader<double>& reader,
               bool& focusHere,
               u32 x,
               u32 y) {
                focusHere = (x == 2 && y == 2) || (x == 6 && y == 6);
                return (reader.read(x, y) + reader.read(x-1, y)) / 2;
            });
```

The second argument is a function that calculates the new value of the `x`, `y`
cell. Currently PatchGraph assumes a 3x3 stencil. `focusHere` is an output
argument through which you can request the current cell to be refined.

That's all for now. You may wonder how can you "unfocus" cells. You can't. Yet.
