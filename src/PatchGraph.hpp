#pragma once

#include "Frac.hpp"
#include <functional>
#include <memory>
#include <vector>

typedef uint8_t Side;
#define LEFT 0x00U
#define RIGHT 0xFFU
#define UP 0xF0U
#define DOWN 0x0FU

template<typename T>
struct Patch;

template<typename T>
struct Section
{
    std::weak_ptr<Patch<T>> leftOrUpPatch, rightOrDownPatch;
    frac1 length, leftOrUpPosition, rightOrDownPosition;
    std::weak_ptr<Patch<T>>& patch(Side side);
    frac1::lview position(Side side);
};

template<typename T>
struct Corner
{
    std::weak_ptr<Patch<T>> patch;
    frac1 position;
    Side side;
    Corner() = default;
    Corner(const std::weak_ptr<Patch<T>>& patch, frac1 position, Side side)
        : patch(patch)
        , position(position)
        , side(side)
    {
    }
};

template<typename T>
struct DataReader
{
    const Patch<T>& patch;
    inline T read(u32 x, u32 y) const;
    DataReader(const Patch<T>& patch)
        : patch(patch)
    {
    }
};

template<typename T>
struct DataWriter
{
    Patch<T>& patch;
    inline void write(u32 x, u32 y, T value);
    DataWriter(Patch<T>& patch)
        : patch(patch)
    {
    }
};
template<typename T>
struct Patch
{
    frac2 dimensions, position;
    std::vector<T> data;
    std::vector<std::shared_ptr<Section<T>>> leftEdge, upEdge, rightEdge,
        downEdge;
    std::shared_ptr<Corner<T>> upLeftCorner, upRightCorner, downLeftCorner,
        downRightCorner;
    std::tuple<std::shared_ptr<Corner<T>>&, std::shared_ptr<Corner<T>>&>
    corners(Side side);
    std::vector<std::shared_ptr<Corner<T>>> cornersOnLeftEdge, cornersOnUpEdge,
        cornersOnRightEdge, cornersOnDownEdge;
    std::vector<std::shared_ptr<Corner<T>>>& cornersOnEdge(Side side);
    std::vector<std::shared_ptr<Section<T>>>& edge(Side side);
    frac2::lview parallelDimension(Side side);
    frac2::lview orthogonalDimension(Side side);
    frac2::lview parallelPosition(Side side);
    frac2::lview orthogonalPosition(Side side);
    void prepareBuffer();
    T read(u32 x, u32 y) const;
    void write(u32 x, u32 y, T value);
    void writeEdge(u32 i, Side side, T value);
    std::tuple<u32, u32, u32, u32, u32>
    synchronizationParameters(FracRView pos, FracRView depth, Side side) const;

    template<typename DownsampleFunc, typename UpsampleFunc>
    void synchronizeSection(std::shared_ptr<Section<T>> section,
                            Side side,
                            DownsampleFunc downsample,
                            UpsampleFunc upsample);
    template<typename DownsampleFunc, typename UpsampleFunc>
    void synchronizeCorner(const std::shared_ptr<Corner<T>>& corner,
                           u32 thisX,
                           u32 thisY,
                           DownsampleFunc downsample,
                           UpsampleFunc upsample);
    template<typename DownsampleFunc, typename UpsampleFunc>
    void synchronizeEdges(DownsampleFunc downsample, UpsampleFunc upsample);

    u32 fracToLength(FracRView frac) const;
    void print() const;
};

template<typename T>
std::pair<std::shared_ptr<Patch<T>>, std::shared_ptr<Patch<T>>> splitAndFocus(
    std::shared_ptr<Patch<T>>&& patch,
    i32 where_i,
    bool vertical,
    u8 luFocus,
    u8 rdFocus);

template<typename T, typename DownsampleFunc, typename UpsampleFunc>
struct PatchGraph
{
    DownsampleFunc downsample;
    UpsampleFunc upsample;
    std::vector<std::shared_ptr<Patch<T>>> patches;
    PatchGraph(u32 width,
               u32 height,
               DownsampleFunc&& downsample,
               UpsampleFunc&& upsample);
    void synchronizeEdges();
    void splitAndFocus(size_t which,
                       size_t where,
                       bool vertical,
                       u8 luFocus,
                       u8 rdFocus);
    void print() const;
};

template<typename T, typename DownsampleFunc, typename UpsampleFunc>
PatchGraph<T, DownsampleFunc, UpsampleFunc> createPatchGraph(
    u32 width,
    u32 height,
    DownsampleFunc&& downsample,
    UpsampleFunc&& upsample)
{
    return PatchGraph<T, DownsampleFunc, UpsampleFunc>(
        width, height, std::move(downsample), std::move(upsample));
}

template<typename T>
void zoomIn(Patch<T>& patch, u8 level);

#define ENABLE_PATCHGRAPH_IMPL_HPP
#include "PatchGraph.impl.hpp"
