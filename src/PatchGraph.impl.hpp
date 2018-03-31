#pragma once

#ifndef ENABLE_PATCHGRAPH_IMPL_HPP
#error "Don't include this file directly!"
#include "PatchGraph.hpp"
#else
#undef ENABLE_PATCHGRAPH_IMPL_HPP
#endif

#include <cassert>

template<typename T>
T DataReader<T>::read(u32 x, u32 y) const
{
    return this->patch.read(x, y);
}

template<typename T>
void DataWriter<T>::write(u32 x, u32 y, T value)
{
    this->patch.write(x, y, value);
}

template<typename T>
std::weak_ptr<Patch<T>>& Section<T>::patch(Side side)
{
    switch (side) {
        case LEFT:
        case UP:
            return leftOrUpPatch;
        case RIGHT:
        case DOWN:
            return rightOrDownPatch;
    }
    assert(false);
}

template<typename T>
frac1::lview Section<T>::position(Side side)
{
    switch (side) {
        case LEFT:
        case UP:
            return this->leftOrUpPosition[0];
        case RIGHT:
        case DOWN:
            return this->rightOrDownPosition[0];
    }
    assert(false);
}

template<typename T>
std::tuple<std::shared_ptr<Corner<T>>&, std::shared_ptr<Corner<T>>&>
Patch<T>::corners(Side side)
{
    switch (side) {
        case LEFT:
            return std::tie(this->upLeftCorner, this->downLeftCorner);
        case RIGHT:
            return std::tie(this->upRightCorner, this->downRightCorner);
        case UP:
            return std::tie(this->upLeftCorner, this->upRightCorner);
        case DOWN:
            return std::tie(this->downLeftCorner, this->downRightCorner);
        default:
            assert(false);
    }
}

template<typename T>
std::vector<std::shared_ptr<Corner<T>>>& Patch<T>::cornersOnEdge(Side side)
{
    switch (side) {
        case LEFT:
            return this->cornersOnLeftEdge;
        case RIGHT:
            return this->cornersOnRightEdge;
        case UP:
            return this->cornersOnUpEdge;
        case DOWN:
            return this->cornersOnDownEdge;
    }
    assert(false);
}

template<typename T>
std::vector<std::shared_ptr<Section<T>>>& Patch<T>::edge(Side side)
{
    switch (side) {
        case LEFT:
            return this->leftEdge;
        case RIGHT:
            return this->rightEdge;
        case UP:
            return this->upEdge;
        case DOWN:
            return this->downEdge;
    }
    assert(false);
}

template<typename T>
frac2::lview Patch<T>::parallelDimension(Side side)
{
    switch (side) {
        case UP:
        case DOWN:
            return this->dimensions[0];
        case LEFT:
        case RIGHT:
            return this->dimensions[1];
    }
    assert(false);
}

template<typename T>
frac2::lview Patch<T>::orthogonalDimension(Side side)
{
    switch (side) {
        case UP:
        case DOWN:
            return this->dimensions[1];
        case LEFT:
        case RIGHT:
            return this->dimensions[0];
    }
    assert(false);
}

template<typename T>
frac2::lview Patch<T>::parallelPosition(Side side)
{
    switch (side) {
        case UP:
        case DOWN:
            return this->position[0];
        case LEFT:
        case RIGHT:
            return this->position[1];
    }
    assert(false);
}

template<typename T>
frac2::lview Patch<T>::orthogonalPosition(Side side)
{
    switch (side) {
        case UP:
        case DOWN:
            return this->position[1];
        case LEFT:
        case RIGHT:
            return this->position[0];
    }
    assert(false);
}

template<typename T>
void Patch<T>::prepareBuffer()
{
    this->data.clear();
    this->data.resize(
        (this->dimensions[0].nom() + 2) * (this->dimensions[1].nom() + 2), 1.0);
}

template<typename T>
T Patch<T>::read(u32 x, u32 y) const
{
    u32 index = (this->dimensions[0].nom() + 2) * y + x;
    assert(index < this->data.size());
    return this->data[index];
}

template<typename T>
void Patch<T>::write(u32 x, u32 y, T value)
{
    u32 index = (this->dimensions[0].nom() + 2) * y + x;
    assert(index < this->data.size());
    this->data[index] = value;
}

template<typename T>
void Patch<T>::writeEdge(u32 i, Side side, T value)
{
    switch (side) {
        case LEFT:
            write(0, i, value);
            break;
        case RIGHT:
            write(this->dimensions[0].nom() + 1, i, value);
            break;
        case UP:
            write(i, 0, value);
            break;
        case DOWN:
            write(i, this->dimensions[1].nom() + 1, value);
            break;
        default:
            assert(false);
    }
}

template<typename T>
std::tuple<u32, u32, u32, u32, u32> Patch<T>::synchronizationParameters(
    FracRView pos,
    FracRView depth,
    Side side) const
{
    u32 cellSize = depth.denom() >= this->dimensions[0].denom()
                       ? 1U
                       : this->fracToLength(depth);
    u32 fromX, fromY, deltaX, deltaY;
    switch (side) {
        case LEFT:
            fromX = 1;
            fromY = 1 + this->fracToLength(pos);
            deltaX = 0;
            deltaY = cellSize;
            break;
        case RIGHT:
            fromX = 1 + this->dimensions[0].nom() - cellSize;
            fromY = 1 + this->fracToLength(pos);
            deltaX = 0;
            deltaY = cellSize;
            break;
        case UP:
            fromX = 1 + this->fracToLength(pos);
            fromY = 1;
            deltaX = cellSize;
            deltaY = 0;
            break;
        case DOWN:
            fromX = 1 + this->fracToLength(pos);
            fromY = 1 + this->dimensions[1].nom() - cellSize;
            deltaX = cellSize;
            deltaY = 0;
            break;
        default:
            assert(false);
    }

    return std::make_tuple(fromX, fromY, deltaX, deltaY, cellSize);
}

template<typename T>
u32 Patch<T>::fracToLength(FracRView frac) const
{
    u32 multiplier = 1;
    u32 divisor = 1;
    if (frac.denom() < this->dimensions.denom()) {
        multiplier *= this->dimensions.denom() / frac.denom();
    } else {
        divisor *= frac.denom() / this->dimensions.denom();
    }
    return frac.nom() * multiplier / divisor;
}

template<typename T>
template<typename DownsampleFunc, typename UpsampleFunc>
void Patch<T>::synchronizeSection(std::shared_ptr<Section<T>> section,
                                  Side side,
                                  DownsampleFunc downsample,
                                  UpsampleFunc upsample)
{
    auto other = section->patch(side).lock();
    u32 otherX, otherY, otherDeltaX, otherDeltaY, otherCellSize;
    std::tie(otherX, otherY, otherDeltaX, otherDeltaY, otherCellSize) =
        other->synchronizationParameters(
            section->position(side), this->dimensions.unit()[0], ~side);

    u32 thisCellSize = otherCellSize == 1U
                           ? this->fracToLength(other->dimensions.unit()[0])
                           : 1U;
    u32 thisX = 0, thisY = 0, thisDeltaX = 0, thisDeltaY = 0;
    switch (side) {
        case RIGHT:
            thisX = this->dimensions[0].nom() + 1;
        case LEFT:
            thisY = 1 + this->fracToLength(section->position(~side));
            thisDeltaY = thisCellSize;
            break;
        case DOWN:
            thisY = this->dimensions[1].nom() + 1;
        case UP:
            thisX = 1 + this->fracToLength(section->position(~side));
            thisDeltaX = thisCellSize;
            break;
        default:
            assert(false);
    }
    u32 thisCellWidth = std::max(1U, thisDeltaX);
    u32 thisCellHeight = std::max(1U, thisDeltaY);

    DataReader<T> reader(*other);
    DataWriter<T> writer(*this);
    for (u32 i = 0; i < this->fracToLength(section->length[0]);
         i += thisCellSize) {
        if (thisCellSize < otherCellSize) {
            this->write(thisX,
                        thisY,
                        downsample(reader,
                                   otherX,
                                   otherY,
                                   otherX + otherCellSize,
                                   otherY + otherCellSize));
        } else {
            upsample(writer,
                     other->read(otherX, otherY),
                     thisX,
                     thisY,
                     thisX + thisCellWidth,
                     thisY + thisCellHeight,
                     thisCellSize);
        }
        thisX += thisDeltaX;
        thisY += thisDeltaY;
        otherX += otherDeltaX;
        otherY += otherDeltaY;
    }
}

template<typename T>
template<typename DownsampleFunc, typename UpsampleFunc>
void Patch<T>::synchronizeCorner(const std::shared_ptr<Corner<T>>& corner,
                                 u32 thisX,
                                 u32 thisY,
                                 DownsampleFunc downsample,
                                 UpsampleFunc upsample)
{
    auto patch = corner->patch.lock();
    u32 otherX, otherY, otherCellSize;
    std::tie(otherX, otherY, std::ignore, std::ignore, otherCellSize) =
        patch->synchronizationParameters(
            corner->position[0], this->dimensions.unit()[0], corner->side);

    u32 thisCellSize = otherCellSize == 1U
                           ? this->fracToLength(patch->dimensions.unit()[0])
                           : 1U;

    DataReader<T> reader(*patch);
    DataWriter<T> writer(*this);
    if (thisCellSize < otherCellSize) {
        this->write(thisX,
                    thisY,
                    downsample(reader,
                               otherX,
                               otherY,
                               otherX + otherCellSize,
                               otherY + otherCellSize));
    } else {
        upsample(writer,
                 patch->read(otherX, otherY),
                 thisX,
                 thisY,
                 thisX + 1,
                 thisY + 1,
                 thisCellSize);
    }
}

template<typename T>
template<typename DownsampleFunc, typename UpsampleFunc>
void Patch<T>::synchronizeEdges(DownsampleFunc downsample,
                                UpsampleFunc upsample)
{
    for (const auto& section : this->leftEdge) {
        this->synchronizeSection(section, LEFT, downsample, upsample);
    }
    for (const auto& section : this->rightEdge) {
        this->synchronizeSection(section, RIGHT, downsample, upsample);
    }
    for (const auto& section : this->upEdge) {
        this->synchronizeSection(section, UP, downsample, upsample);
    }
    for (const auto& section : this->downEdge) {
        this->synchronizeSection(section, DOWN, downsample, upsample);
    }
    if (this->upLeftCorner != nullptr) {

        this->synchronizeCorner(this->upLeftCorner, 0, 0, downsample, upsample);
    }
    if (this->downLeftCorner != nullptr) {
        this->synchronizeCorner(this->downLeftCorner,
                                0,
                                this->dimensions[1].nom() + 1,
                                downsample,
                                upsample);
    }
    if (this->upRightCorner != nullptr) {
        this->synchronizeCorner(this->upRightCorner,
                                this->dimensions[0].nom() + 1,
                                0,
                                downsample,
                                upsample);
    }
    if (this->downRightCorner != nullptr) {
        this->synchronizeCorner(this->downRightCorner,
                                this->dimensions[0].nom() + 1,
                                this->dimensions[1].nom() + 1,
                                downsample,
                                upsample);
    }
}

template<typename T>
void Patch<T>::print() const
{
    for (u32 y = 0; y < this->dimensions[1].nom() + 2; ++y) {
        for (u32 x = 0; x < this->dimensions[0].nom() + 2; ++x) {
            std::cout << this->read(x, y) << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

template<typename T>
void zoomIn(Patch<T>& patch, u8 level)
{
    patch.dimensions.focus(level);
    patch.prepareBuffer();
}

template<typename T, typename DownsampleFunc, typename UpsampleFunc>
PatchGraph<T, DownsampleFunc, UpsampleFunc>::PatchGraph(
    u32 width,
    u32 height,
    DownsampleFunc&& downsample,
    UpsampleFunc&& upsample)
    : downsample(std::move(downsample))
    , upsample(std::move(upsample))
{
    auto patch = std::make_shared<Patch<T>>();
    patch->dimensions = frac2(width, height);
    patch->prepareBuffer();
    this->patches.push_back(std::move(patch));
}

template<typename T, typename DownsampleFunc, typename UpsampleFunc>
void PatchGraph<T, DownsampleFunc, UpsampleFunc>::synchronizeEdges()
{
    for (const auto& patch : this->patches) {
        patch->synchronizeEdges(this->downsample, this->upsample);
    }
}

template<typename T>
size_t splitEdge(std::vector<std::shared_ptr<Section<T>>>& edge,
                 frac1 where,
                 Side side)
{
    if (edge.size() == 0) {
        return 0;
    }
    frac1 length({ 0 });
    for (size_t i = 0; i < edge.size(); ++i) {
        frac1 sectionLength = edge[i]->length;
        if (length[0] + sectionLength[0] > where[0]) {
            frac1 overflow;
            overflow[0] = length[0] + sectionLength[0] - where[0];
            if (overflow[0] == sectionLength[0]) {
                return i;
            }
            edge[i]->length[0] = sectionLength[0] - overflow[0];
            auto newSection = std::make_shared<Section<T>>(*edge[i]);
            newSection->length[0] = frac1(0U)[0] + overflow[0];
            newSection->leftOrUpPosition[0] += sectionLength[0] - overflow[0];
            newSection->rightOrDownPosition[0] +=
                sectionLength[0] - overflow[0];
            edge.insert(edge.begin() + i + 1, newSection);
            auto otherPatch = edge[i]->patch(side);
            auto& otherEdge = otherPatch.lock()->edge(~side);
            for (size_t j = 0; j < otherEdge.size(); ++j) {
                if (otherEdge[j] == edge[i]) {
                    otherEdge.insert(otherEdge.begin() + j + 1,
                                     std::move(newSection));
                    return i + 1;
                }
            }
            assert(false);
        }
    }
    assert(false);
}

template<typename T, typename DownsampleFunc, typename UpsampleFunc>
void PatchGraph<T, DownsampleFunc, UpsampleFunc>::splitAndFocus(size_t which,
                                                                size_t where_,
                                                                bool vertical,
                                                                u8 luFocus,
                                                                u8 rdFocus)
{
    auto patch = this->patches[which];
    this->patches.erase(this->patches.begin() + which);

    Side side = vertical ? RIGHT : DOWN;
    Side pside = vertical ? UP : LEFT;
    frac1 where;
    where[0] = patch->dimensions.unit()[0] * where_;
    assert(where_ > 0);
    assert(where[0] < patch->orthogonalDimension(side));
    assert(where_ % patch->parallelDimension(side).denom() == 0);
    auto luPatch = std::make_shared<Patch<T>>();
    auto rdPatch = std::make_shared<Patch<T>>();
    luPatch->dimensions = patch->dimensions;
    luPatch->orthogonalDimension(side) = where[0];
    rdPatch->dimensions = patch->dimensions;
    rdPatch->orthogonalDimension(side) =
        patch->orthogonalDimension(side) - where[0];
    luPatch->position = patch->position;
    rdPatch->position = patch->position;
    rdPatch->orthogonalPosition(side) += where[0];
    luPatch->corners(~side) = patch->corners(~side);
    rdPatch->corners(side) = patch->corners(side);
    luPatch->cornersOnEdge(~side) = std::move(patch->cornersOnEdge(~side));
    for (auto& corner : luPatch->cornersOnEdge(~side)) {
        corner->patch = luPatch;
    }
    rdPatch->cornersOnEdge(side) = std::move(patch->cornersOnEdge(side));
    for (auto& corner : rdPatch->cornersOnEdge(side)) {
        corner->patch = rdPatch;
    }

    for (const auto& section : patch->edge(side)) {
        section->patch(side) = rdPatch;
        section->patch(~side) = luPatch;
    }

    size_t luSplitEdge = splitEdge(patch->edge(pside), where, pside);
    luPatch->edge(pside).assign(patch->edge(pside).begin(),
                                patch->edge(pside).begin() + luSplitEdge);
    rdPatch->edge(pside).assign(patch->edge(pside).begin() + luSplitEdge,
                                patch->edge(pside).end());
    for (auto& corner : patch->cornersOnEdge(pside)) {
        if (corner->position[0] < luPatch->parallelDimension(pside)) {
            corner->patch = luPatch;
            luPatch->cornersOnEdge(pside).push_back(std::move(corner));
        } else {
            corner->position[0] -= luPatch->parallelDimension(pside);
            corner->patch = rdPatch;
            rdPatch->cornersOnEdge(pside).push_back(std::move(corner));
        }
    }
    for (auto& corner : patch->cornersOnEdge(~pside)) {
        if (corner->position[0] < luPatch->parallelDimension(~pside)) {
            corner->patch = luPatch;
            luPatch->cornersOnEdge(~pside).push_back(std::move(corner));
        } else {
            corner->position[0] -= luPatch->parallelDimension(~pside);
            corner->patch = rdPatch;
            rdPatch->cornersOnEdge(~pside).push_back(std::move(corner));
        }
    }
    for (auto& section : luPatch->edge(pside)) {
        section->rightOrDownPatch = luPatch;
    }
    frac1 nearSideOffset;
    for (auto& section : rdPatch->edge(pside)) {
        section->rightOrDownPosition = nearSideOffset;
        section->rightOrDownPatch = rdPatch;
        nearSideOffset[0] += section->length[0];
    }
    size_t rdSplitEdge = splitEdge(patch->edge(~pside), where, ~pside);
    luPatch->edge(~pside).assign(patch->edge(~pside).begin(),
                                 patch->edge(~pside).begin() + rdSplitEdge);
    rdPatch->edge(~pside).assign(patch->edge(~pside).begin() + rdSplitEdge,
                                 patch->edge(~pside).end());
    for (auto& section : luPatch->edge(~pside)) {
        section->leftOrUpPatch = luPatch;
    }
    frac1 farSideOffset;
    for (auto& section : rdPatch->edge(~pside)) {
        section->leftOrUpPosition = farSideOffset;
        section->leftOrUpPatch = rdPatch;
        farSideOffset[0] += section->length[0];
    }

    auto luCorners = luPatch->corners(side);
    if (!rdPatch->edge(pside).empty()) {
        auto luSection = rdPatch->edge(pside)[0];
        std::get<0>(luCorners) = std::make_shared<Corner<T>>(
            luSection->leftOrUpPatch, luSection->leftOrUpPosition, ~pside);
    }
    if (!rdPatch->edge(~pside).empty()) {
        auto rdSection = rdPatch->edge(~pside)[0];
        std::get<1>(luCorners) = std::make_shared<Corner<T>>(
            rdSection->rightOrDownPatch, rdSection->rightOrDownPosition, pside);
    }
    auto rdCorners = rdPatch->corners(~side);
    if (!luPatch->edge(pside).empty()) {
        auto luSection = luPatch->edge(pside).back();
        std::get<0>(rdCorners) = std::make_shared<Corner<T>>(
            luSection->leftOrUpPatch, luSection->leftOrUpPosition, ~pside);
        std::get<0>(rdCorners)->position[0] +=
            luSection->leftOrUpPosition[0] + luPatch->dimensions.unit()[0];
    }
    if (!luPatch->edge(~pside).empty()) {
        auto rdSection = luPatch->edge(~pside)[0];
        std::get<1>(rdCorners) = std::make_shared<Corner<T>>(
            rdSection->rightOrDownPatch, rdSection->rightOrDownPosition, pside);
        std::get<1>(rdCorners)->position[0] +=
            rdSection->rightOrDownPosition[0] + rdPatch->dimensions.unit()[0];
    }

    auto sharedEdge = std::make_shared<Section<T>>();
    sharedEdge->patch(~side) = luPatch;
    sharedEdge->patch(side) = rdPatch;
    sharedEdge->length[0] = luPatch->parallelDimension(side);
    luPatch->edge(side).push_back(sharedEdge);
    rdPatch->edge(~side).push_back(sharedEdge);
    zoomIn(*luPatch, luFocus);
    zoomIn(*rdPatch, rdFocus);
    this->patches.push_back(std::move(luPatch));
    this->patches.push_back(std::move(rdPatch));
}

template<typename T, typename DownsampleFunc, typename UpsampleFunc>
void PatchGraph<T, DownsampleFunc, UpsampleFunc>::print() const
{
    for (const auto& patch : this->patches) {
        patch->print();
    }
}
