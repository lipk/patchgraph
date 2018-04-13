#pragma once

#include "Frac.hpp"
#include <functional>
#include <memory>
#include <vector>

template<typename T>
struct Patch;

template<typename T>
struct DataReader
{
    const Patch<T>& patch;
    inline T read(u32 x, u32 y) const;
    DataReader(const Patch<T>& patch)
        : patch(patch)
    {}
};

template<typename T>
struct DataWriter
{
    Patch<T>& patch;
    inline void write(u32 x, u32 y, T value);
    DataWriter(Patch<T>& patch)
        : patch(patch)
    {}
};

template<typename T, typename DownsampleFunc, typename UpsampleFunc>
class PatchGraph
{
    struct Grid
    {
        std::shared_ptr<Patch<T>> leftBound, rightBound, upBound, downBound;
        std::vector<std::shared_ptr<Patch<T>>> middle;

        Grid(u32 width, u32 height);
    };

    DownsampleFunc downsample;
    UpsampleFunc upsample;
    Grid grid1, grid2;
    Grid *source, *target;
    void synchronizeEdges();

    void focusAtPoints(
        size_t patchIndex,
        const std::vector<std::pair<u32, u32>>& points,
        std::vector<std::shared_ptr<Patch<T>>>& newSourcePatches,
        std::vector<std::shared_ptr<Patch<T>>>& newTargetPatches);

    void defocus(size_t patchIndex);

    std::tuple<std::shared_ptr<Patch<T>>, u32, u32> find(u32 x, u32 y) const;

  public:
    PatchGraph(u32 width,
               u32 height,
               DownsampleFunc&& downsample,
               UpsampleFunc&& upsample);

    void write(u32 x, u32 y, T value);
    T read(u32 x, u32 y) const;

    template<typename StencilFunc>
    void apply(size_t times, StencilFunc func);
    template<typename StencilFunc>
    void leftBound(StencilFunc func);
    template<typename StencilFunc>
    void rightBound(StencilFunc func);
    template<typename StencilFunc>
    void upBound(StencilFunc func);
    template<typename StencilFunc>
    void downBound(StencilFunc func);

    // TODO: remove this
    std::vector<std::shared_ptr<Patch<T>>>& getSource()
    {
        return this->source->middle;
    }

    void print() const;
};

template<typename T, typename DownsampleFunc, typename UpsampleFunc>
PatchGraph<T, DownsampleFunc, UpsampleFunc> createPatchGraph(
    u32 width,
    u32 height,
    DownsampleFunc downsample,
    UpsampleFunc upsample)
{
    return PatchGraph<T, DownsampleFunc, UpsampleFunc>(
        width, height, std::move(downsample), std::move(upsample));
}

#define ENABLE_PATCHGRAPH_IMPL_HPP
#include "PatchGraph.impl.hpp"
