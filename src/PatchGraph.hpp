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

struct Patch;
struct Section
{
    std::weak_ptr<Patch> leftOrUpPatch, rightOrDownPatch;
    frac1 length, leftOrUpPosition, rightOrDownPosition;
    std::weak_ptr<Patch>& patch(Side side);
    frac1::lview position(Side side);
};

struct Corner
{
    std::weak_ptr<Patch> patch;
    frac1 position;
    Corner() = default;
    Corner(const std::weak_ptr<Patch>& patch, frac1 position)
        : patch(patch)
        , position(position)
    {
    }
};

struct Patch
{
    typedef double T;
    frac2 dimensions, position;
    std::vector<T> data;
    std::vector<std::shared_ptr<Section>> leftEdge, upEdge, rightEdge, downEdge;
    std::shared_ptr<Corner> upLeftCorner, upRightCorner, downLeftCorner,
        downRightCorner;
    std::tuple<std::shared_ptr<Corner>&, std::shared_ptr<Corner>&> corners(
        Side side);
    std::vector<std::shared_ptr<Corner>> cornersOnLeftEdge, cornersOnUpEdge,
        cornersOnRightEdge, cornersOnDownEdge;
    std::vector<std::shared_ptr<Corner>>& cornersOnEdge(Side side);
    std::vector<std::shared_ptr<Section>>& edge(Side side);
    frac2::lview parallelDimension(Side side);
    frac2::lview orthogonalDimension(Side side);
    frac2::lview parallelPosition(Side side);
    frac2::lview orthogonalPosition(Side side);
    void prepareBuffer();
    T read(u32 x, u32 y) const;
    void write(u32 x, u32 y, T value);
    void writeEdge(u32 i, Side side, T value);
    T sample(FracRView pos, FracRView size, Side side) const;
    T sum(frac2 pos, frac2 size) const;
    void synchronizeSection(std::shared_ptr<Section> section, Side side);
    u32 fracToLength(FracRView frac) const;
    void synchronizeEdges();
    void print() const;
};

std::pair<std::shared_ptr<Patch>, std::shared_ptr<Patch>> splitAndFocus(
    std::shared_ptr<Patch>&& patch,
    i32 where_i,
    bool vertical,
    u8 luFocus,
    u8 rdFocus);

struct DoublePatchGraph
{
    std::vector<std::shared_ptr<Patch>> patches;
    DoublePatchGraph(u32 width, u32 height);
    void synchronizeEdges();
    void print() const;
};

void zoomIn(Patch& patch, u8 level);
