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
        case Side::Left:
        case Side::Up:
            return leftOrUpPatch;
        case Side::Right:
        case Side::Down:
            return rightOrDownPatch;
    }
    assert(false);
}

template<typename T>
frac1::lview Section<T>::position(Side side)
{
    switch (side) {
        case Side::Left:
        case Side::Up:
            return this->leftOrUpPosition[0];
        case Side::Right:
        case Side::Down:
            return this->rightOrDownPosition[0];
    }
    assert(false);
}

template<typename T>
std::tuple<std::shared_ptr<Corner<T>>&, std::shared_ptr<Corner<T>>&>
Patch<T>::corners(Side side)
{
    switch (side) {
        case Side::Left:
            return std::tie(this->upLeftCorner, this->downLeftCorner);
        case Side::Right:
            return std::tie(this->upRightCorner, this->downRightCorner);
        case Side::Up:
            return std::tie(this->upLeftCorner, this->upRightCorner);
        case Side::Down:
            return std::tie(this->downLeftCorner, this->downRightCorner);
        default:
            assert(false);
    }
}

template<typename T>
std::vector<std::shared_ptr<Corner<T>>>& Patch<T>::cornersOnEdge(Side side)
{
    switch (side) {
        case Side::Left:
            return this->cornersOnLeftEdge;
        case Side::Right:
            return this->cornersOnRightEdge;
        case Side::Up:
            return this->cornersOnUpEdge;
        case Side::Down:
            return this->cornersOnDownEdge;
    }
    assert(false);
}

template<typename T>
std::vector<std::shared_ptr<Section<T>>>& Patch<T>::edge(Side side)
{
    switch (side) {
        case Side::Left:
            return this->leftEdge;
        case Side::Right:
            return this->rightEdge;
        case Side::Up:
            return this->upEdge;
        case Side::Down:
            return this->downEdge;
    }
    assert(false);
}

template<typename T>
frac2::lview Patch<T>::parallelDimension(Side side)
{
    switch (side) {
        case Side::Up:
        case Side::Down:
            return this->dimensions[0];
        case Side::Left:
        case Side::Right:
            return this->dimensions[1];
    }
    assert(false);
}

template<typename T>
frac2::lview Patch<T>::orthogonalDimension(Side side)
{
    switch (side) {
        case Side::Up:
        case Side::Down:
            return this->dimensions[1];
        case Side::Left:
        case Side::Right:
            return this->dimensions[0];
    }
    assert(false);
}

template<typename T>
frac2::lview Patch<T>::parallelPosition(Side side)
{
    switch (side) {
        case Side::Up:
        case Side::Down:
            return this->position[0];
        case Side::Left:
        case Side::Right:
            return this->position[1];
    }
    assert(false);
}

template<typename T>
frac2::lview Patch<T>::orthogonalPosition(Side side)
{
    switch (side) {
        case Side::Up:
        case Side::Down:
            return this->position[1];
        case Side::Left:
        case Side::Right:
            return this->position[0];
    }
    assert(false);
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
        case Side::Left:
            write(0, i, value);
            break;
        case Side::Right:
            write(this->dimensions[0].nom() + 1, i, value);
            break;
        case Side::Up:
            write(i, 0, value);
            break;
        case Side::Down:
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
        case Side::Left:
            fromX = 1;
            fromY = 1 + this->fracToLength(pos);
            deltaX = 0;
            deltaY = cellSize;
            break;
        case Side::Right:
            fromX = 1 + this->dimensions[0].nom() - cellSize;
            fromY = 1 + this->fracToLength(pos);
            deltaX = 0;
            deltaY = cellSize;
            break;
        case Side::Up:
            fromX = 1 + this->fracToLength(pos);
            fromY = 1;
            deltaX = cellSize;
            deltaY = 0;
            break;
        case Side::Down:
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
        case Side::Right:
            thisX = this->dimensions[0].nom() + 1;
        case Side::Left:
            thisY = 1 + this->fracToLength(section->position(~side));
            thisDeltaY = thisCellSize;
            break;
        case Side::Down:
            thisY = this->dimensions[1].nom() + 1;
        case Side::Up:
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
        this->synchronizeSection(section, Side::Left, downsample, upsample);
    }
    for (const auto& section : this->rightEdge) {
        this->synchronizeSection(section, Side::Right, downsample, upsample);
    }
    for (const auto& section : this->upEdge) {
        this->synchronizeSection(section, Side::Up, downsample, upsample);
    }
    for (const auto& section : this->downEdge) {
        this->synchronizeSection(section, Side::Down, downsample, upsample);
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
template<typename UpsampleFunc>
void Patch<T>::focus(u8 level, UpsampleFunc upsample)
{
    if (level == 0) {
        return;
    }
    auto oldData = std::move(this->data);
    const u32 ratio = 1 << level;
    const u32 oldWidth = this->dimensions[0].nom();
    const u32 oldHeight = this->dimensions[1].nom();
    this->data.resize((ratio * oldWidth + 2) * (ratio * oldHeight + 2));
    DataWriter<T> writer(*this);
    this->dimensions.focus(level);
    for (u32 y = 0; y < oldHeight; ++y) {
        for (u32 x = 0; x < oldWidth; ++x) {
            upsample(writer,
                     oldData[(oldWidth + 2) * (y + 1) + x + 1],
                     1 + x * ratio,
                     1 + y * ratio,
                     1 + (x + 1) * ratio,
                     1 + (y + 1) * ratio,
                     ratio);
        }
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

template<typename T, typename DownsampleFunc, typename UpsampleFunc>
PatchGraph<T, DownsampleFunc, UpsampleFunc>::PatchGraph(
    u32 width,
    u32 height,
    DownsampleFunc&& downsample,
    UpsampleFunc&& upsample)
    : downsample(std::move(downsample))
    , upsample(std::move(upsample))
    , source(&this->patches1)
    , target(&this->patches2)
{
    auto patch1 = std::make_shared<Patch<T>>();
    patch1->dimensions = frac2(width, height);
    patch1->data.resize((width + 2) * (height + 2));
    this->source->push_back(std::move(patch1));

    auto patch2 = std::make_shared<Patch<T>>();
    patch2->dimensions = frac2(width, height);
    patch2->data.resize((width + 2) * (height + 2));
    this->target->push_back(std::move(patch2));
}

template<typename T, typename DownsampleFunc, typename UpsampleFunc>
void PatchGraph<T, DownsampleFunc, UpsampleFunc>::synchronizeEdges()
{
    for (const auto& patch : *this->source) {
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
void PatchGraph<T, DownsampleFunc, UpsampleFunc>::split(size_t which,
                                                        size_t where_,
                                                        bool vertical)
{
    auto patch1 = this->patches1[which];
    this->patches1.erase(this->patches1.begin() + which);
    auto newPatches1 = this->split(std::move(patch1), where_, vertical);
    this->patches1.push_back(std::move(newPatches1.first));
    this->patches1.push_back(std::move(newPatches1.second));

    auto patch2 = this->patches2[which];
    this->patches2.erase(this->patches2.begin() + which);
    auto newPatches2 = this->split(std::move(patch2), where_, vertical);
    this->patches2.push_back(std::move(newPatches2.first));
    this->patches2.push_back(std::move(newPatches2.second));
}

template<typename T, typename DownsampleFunc, typename UpsampleFunc>
std::pair<std::shared_ptr<Patch<T>>, std::shared_ptr<Patch<T>>>
PatchGraph<T, DownsampleFunc, UpsampleFunc>::split(
    std::shared_ptr<Patch<T>> patch,
    size_t where_,
    bool vertical)
{
    // the side parallel with the split
    Side parside = vertical ? Side::Right : Side::Down;
    // the side orthogonal to the split
    Side ortside = vertical ? Side::Up : Side::Left;

    frac1 where;
    where[0] = patch->dimensions.unit()[0] * where_;
    assert(where_ > 0);
    assert(where_ < patch->orthogonalDimension(parside).nom());
    assert(where_ % patch->parallelDimension(parside).denom() == 0);
    auto luPatch = std::make_shared<Patch<T>>();
    auto rdPatch = std::make_shared<Patch<T>>();

    // calculate dimensions
    luPatch->dimensions = patch->dimensions;
    luPatch->orthogonalDimension(parside) = where[0];
    rdPatch->dimensions = patch->dimensions;
    rdPatch->orthogonalDimension(parside) =
        patch->orthogonalDimension(parside) - where[0];
    luPatch->position = patch->position;
    rdPatch->position = patch->position;
    rdPatch->orthogonalPosition(parside) += where[0];

    // copy data
    luPatch->data.resize((luPatch->dimensions[0].nom() + 2) *
                         (luPatch->dimensions[1].nom() + 2));
    for (u32 y = 1; y < luPatch->dimensions[1].nom() + 1; ++y) {
        for (u32 x = 1; x < luPatch->dimensions[0].nom() + 1; ++x) {
            luPatch->write(x, y, patch->read(x, y));
        }
    }
    rdPatch->data.resize((rdPatch->dimensions[0].nom() + 2) *
                         (rdPatch->dimensions[1].nom() + 2));
    u32 xoff = vertical ? luPatch->dimensions[0].nom() : 0;
    u32 yoff = vertical ? 0 : luPatch->dimensions[1].nom();
    for (u32 y = 1; y < rdPatch->dimensions[1].nom() + 1; ++y) {
        for (u32 x = 1; x < rdPatch->dimensions[0].nom() + 1; ++x) {
            rdPatch->write(x, y, patch->read(x + xoff, y + yoff));
        }
    }

    // corners & sections on the parallel sides are moved to the new patches
    luPatch->corners(~parside) = patch->corners(~parside);
    rdPatch->corners(parside) = patch->corners(parside);
    luPatch->cornersOnEdge(~parside) =
        std::move(patch->cornersOnEdge(~parside));
    for (auto& corner : luPatch->cornersOnEdge(~parside)) {
        corner->patch = luPatch;
    }
    rdPatch->cornersOnEdge(parside) = std::move(patch->cornersOnEdge(parside));
    for (auto& corner : rdPatch->cornersOnEdge(parside)) {
        corner->patch = rdPatch;
    }

    for (const auto& section : patch->edge(parside)) {
        section->patch(parside) = rdPatch;
        section->patch(~parside) = luPatch;
    }

    // split one orthogonal side
    size_t luSplitEdge = splitEdge(patch->edge(ortside), where, ortside);
    luPatch->edge(ortside).assign(patch->edge(ortside).begin(),
                                  patch->edge(ortside).begin() + luSplitEdge);
    rdPatch->edge(ortside).assign(patch->edge(ortside).begin() + luSplitEdge,
                                  patch->edge(ortside).end());

    // modify existing corners to point into the new patches
    for (auto& corner : patch->cornersOnEdge(ortside)) {
        if (corner->position[0] < luPatch->parallelDimension(ortside)) {
            corner->patch = luPatch;
            luPatch->cornersOnEdge(ortside).push_back(std::move(corner));
        } else {
            corner->position[0] -= luPatch->parallelDimension(ortside);
            corner->patch = rdPatch;
            rdPatch->cornersOnEdge(ortside).push_back(std::move(corner));
        }
    }

    // move & adjust sections
    for (auto& section : luPatch->edge(ortside)) {
        section->rightOrDownPatch = luPatch;
    }
    frac1 nearSideOffset;
    for (auto& section : rdPatch->edge(ortside)) {
        section->rightOrDownPosition = nearSideOffset;
        section->rightOrDownPatch = rdPatch;
        nearSideOffset[0] += section->length[0];
    }

    // split the other orthogonal side
    size_t rdSplitEdge = splitEdge(patch->edge(~ortside), where, ~ortside);
    luPatch->edge(~ortside).assign(patch->edge(~ortside).begin(),
                                   patch->edge(~ortside).begin() + rdSplitEdge);
    rdPatch->edge(~ortside).assign(patch->edge(~ortside).begin() + rdSplitEdge,
                                   patch->edge(~ortside).end());

    // modify existing corners to point into the new patches
    for (auto& corner : patch->cornersOnEdge(~ortside)) {
        if (corner->position[0] < luPatch->parallelDimension(~ortside)) {
            corner->patch = luPatch;
            luPatch->cornersOnEdge(~ortside).push_back(std::move(corner));
        } else {
            corner->position[0] -= luPatch->parallelDimension(~ortside);
            corner->patch = rdPatch;
            rdPatch->cornersOnEdge(~ortside).push_back(std::move(corner));
        }
    }

    // move & adjust sections
    for (auto& section : luPatch->edge(~ortside)) {
        section->leftOrUpPatch = luPatch;
    }
    frac1 farSideOffset;
    for (auto& section : rdPatch->edge(~ortside)) {
        section->leftOrUpPosition = farSideOffset;
        section->leftOrUpPatch = rdPatch;
        farSideOffset[0] += section->length[0];
    }

    // create the new corners
    auto luCorners = luPatch->corners(parside);
    if (!rdPatch->edge(ortside).empty()) {
        auto luSection = rdPatch->edge(ortside)[0];
        std::get<0>(luCorners) = std::make_shared<Corner<T>>(
            luSection->leftOrUpPatch, luSection->leftOrUpPosition, ~ortside);
    }
    if (!rdPatch->edge(~ortside).empty()) {
        auto rdSection = rdPatch->edge(~ortside)[0];
        std::get<1>(luCorners) =
            std::make_shared<Corner<T>>(rdSection->rightOrDownPatch,
                                        rdSection->rightOrDownPosition,
                                        ortside);
    }
    auto rdCorners = rdPatch->corners(~parside);
    if (!luPatch->edge(ortside).empty()) {
        auto luSection = luPatch->edge(ortside).back();
        std::get<0>(rdCorners) = std::make_shared<Corner<T>>(
            luSection->leftOrUpPatch, luSection->leftOrUpPosition, ~ortside);
        std::get<0>(rdCorners)->position[0] +=
            luSection->leftOrUpPosition[0] + luPatch->dimensions.unit()[0];
    }
    if (!luPatch->edge(~ortside).empty()) {
        auto rdSection = luPatch->edge(~ortside)[0];
        std::get<1>(rdCorners) =
            std::make_shared<Corner<T>>(rdSection->rightOrDownPatch,
                                        rdSection->rightOrDownPosition,
                                        ortside);
        std::get<1>(rdCorners)->position[0] +=
            rdSection->rightOrDownPosition[0] + rdPatch->dimensions.unit()[0];
    }

    // create the new edge between the patches
    auto sharedEdge = std::make_shared<Section<T>>();
    sharedEdge->patch(~parside) = luPatch;
    sharedEdge->patch(parside) = rdPatch;
    sharedEdge->length[0] = luPatch->parallelDimension(parside);
    luPatch->edge(parside).push_back(sharedEdge);
    rdPatch->edge(~parside).push_back(sharedEdge);

    return { std::move(luPatch), std::move(rdPatch) };
}

template<typename T, typename DownsampleFunc, typename UpsampleFunc>
std::tuple<std::shared_ptr<Patch<T>>, u32, u32>
PatchGraph<T, DownsampleFunc, UpsampleFunc>::find(u32 x, u32 y) const
{
    frac1 fx(x);
    frac1 fy(y);
    for (const auto& patch : *this->source) {
        if (fx[0] >= patch->position[0] && fy[0] >= patch->position[1] &&
            fx[0] <= patch->position[0] + patch->dimensions[0] &&
            fy[0] <= patch->position[1] + patch->dimensions[1]) {
            fx[0] -= patch->position[0];
            fy[0] -= patch->position[1];

            u32 resX = 1 + fx.nom(0) / fx.denom() * patch->position.denom();
            u32 resY = 1 + fy.nom(0) / fy.denom() * patch->position.denom();

            return { patch, resX, resY };
        }
    }
    assert(false);
}

template<typename T, typename DownsampleFunc, typename UpsampleFunc>
void PatchGraph<T, DownsampleFunc, UpsampleFunc>::write(u32 x, u32 y, T value)
{
    std::shared_ptr<Patch<T>> patch;
    u32 patchX, patchY;
    std::tie(patch, patchX, patchY) = this->find(x, y);
    DataWriter<T> writer(*patch);
    u32 cellSize = patch->dimensions.denom();
    this->upsample(writer,
                   value,
                   patchX,
                   patchY,
                   patchX + cellSize,
                   patchY + cellSize,
                   cellSize);
}

template<typename T, typename DownsampleFunc, typename UpsampleFunc>
T PatchGraph<T, DownsampleFunc, UpsampleFunc>::read(u32 x, u32 y) const
{
    std::shared_ptr<Patch<T>> patch;
    u32 patchX, patchY;
    std::tie(patch, patchX, patchY) = this->find(x, y);
    DataReader<T> reader(*patch);
    u32 cellSize = 1 << patch->dimensions.denom();
    return this->downsample(
        reader, patchX, patchY, patchX + cellSize, patchY + cellSize);
}

template<typename T, typename DownsampleFunc, typename UpsampleFunc>
void PatchGraph<T, DownsampleFunc, UpsampleFunc>::focusAtPoints(
    size_t patchIndex,
    const std::vector<std::pair<u32, u32>>& points,
    std::vector<std::shared_ptr<Patch<T>>>& newSourcePatches,
    std::vector<std::shared_ptr<Patch<T>>>& newTargetPatches)
{
    assert(!points.empty());
    auto sourcePatch = (*this->source)[patchIndex];
    this->source->erase(this->source->begin() + patchIndex);
    auto targetPatch = (*this->target)[patchIndex];
    this->target->erase(this->target->begin() + patchIndex);
    const u32 width = sourcePatch->dimensions.nom(0);
    const u32 height = targetPatch->dimensions.nom(0);

    // TODO: better algorithm
    u32 xmin = width, xmax = 0, ymin = height, ymax = 0;
    for (auto point : points) {
        xmin = std::min(xmin, point.first);
        xmax = std::max(xmax, point.first + 1);
        ymin = std::min(ymin, point.second);
        ymax = std::max(ymax, point.second + 1);
    }

    u32 denom = sourcePatch->dimensions.denom();
    xmin = xmin / denom * denom;
    ymin = ymin / denom * denom;
    xmax = xmax % denom == 0 ? xmax : (xmax / denom + 1) * denom;
    xmax = ymax % denom == 0 ? ymax : (ymax / denom + 1) * denom;

    auto doSplit = [&](std::shared_ptr<Patch<T>>&& patch,
                       std::vector<std::shared_ptr<Patch<T>>>& patchList) {
        std::shared_ptr<Patch<T>> luPatch, rdPatch;
        if (xmax < width) {
            std::tie(luPatch, rdPatch) =
                this->split(std::move(patch), xmax, true);
            patchList.push_back(std::move(rdPatch));
        } else {
            luPatch = std::move(patch);
        }
        if (xmin > 0) {
            std::tie(luPatch, rdPatch) =
                this->split(std::move(luPatch), xmin, true);
            patchList.push_back(std::move(luPatch));
        } else {
            rdPatch = std::move(luPatch);
        }
        if (ymax < height) {
            std::tie(luPatch, rdPatch) =
                this->split(std::move(rdPatch), ymax, false);
            patchList.push_back(std::move(rdPatch));
        } else {
            luPatch = std::move(rdPatch);
        }
        if (ymin > 0) {
            std::tie(luPatch, rdPatch) =
                this->split(std::move(luPatch), ymin, false);
            rdPatch->focus(1, this->upsample);
            patchList.push_back(std::move(luPatch));
            patchList.push_back(std::move(rdPatch));
        } else {
            luPatch->focus(1, this->upsample);
            patchList.push_back(std::move(luPatch));
        }
    };

    doSplit(std::move(sourcePatch), newSourcePatches);
    doSplit(std::move(targetPatch), newTargetPatches);
}

template<typename T, typename DownsampleFunc, typename UpsampleFunc>
template<typename StencilFunc>
void PatchGraph<T, DownsampleFunc, UpsampleFunc>::apply(size_t times,
                                                        StencilFunc func)
{
    assert(this->source->size() == this->target->size());
    std::vector<std::pair<u32, u32>> focusPoints;
    std::vector<std::shared_ptr<Patch<T>>> newSourcePatches, newTargetPatches;
    for (size_t t = 0; t < times; ++t) {
        this->synchronizeEdges();
        for (size_t i = 0; i < this->source->size(); ++i) {
            auto sourcePatch = (*this->source)[i];
            auto targetPatch = (*this->target)[i];
            DataReader<T> reader(*sourcePatch);
            focusPoints.clear();
            for (u32 y = 1; y < sourcePatch->dimensions[1].nom() + 1; ++y) {
                for (u32 x = 1; x < targetPatch->dimensions[0].nom() + 1; ++x) {
                    bool focusHere = false;
                    targetPatch->write(x, y, func(reader, focusHere, x, y));
                    if (focusHere) {
                        focusPoints.emplace_back(x - 1, y - 1);
                    }
                }
            }
            if (!focusPoints.empty()) {
                this->focusAtPoints(
                    i, focusPoints, newSourcePatches, newTargetPatches);
                i--;
            }
        }
        this->source->insert(this->source->end(),
                             newSourcePatches.begin(),
                             newSourcePatches.end());
        newSourcePatches.clear();
        this->target->insert(this->target->end(),
                             newTargetPatches.begin(),
                             newTargetPatches.end());
        newTargetPatches.clear();
        std::swap(this->source, this->target);
    }
}

template<typename T, typename DownsampleFunc, typename UpsampleFunc>
void PatchGraph<T, DownsampleFunc, UpsampleFunc>::print() const
{
    for (const auto& patch : this->patches1) {
        patch->print();
    }

    for (const auto& patch : this->patches2) {
        patch->print();
    }
}
