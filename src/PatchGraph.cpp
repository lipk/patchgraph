#include "PatchGraph.hpp"

static size_t splitEdge(std::vector<std::shared_ptr<Section>>& edge,
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
            auto newSection = std::make_shared<Section>(*edge[i]);
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

std::pair<std::shared_ptr<Patch>, std::shared_ptr<Patch>> splitAndFocus(
    std::shared_ptr<Patch>&& patch,
    i32 where_i,
    bool vertical,
    u8 luFocus,
    u8 rdFocus)
{
    Side side = vertical ? RIGHT : DOWN;
    Side pside = vertical ? UP : LEFT;
    frac1 where;
    where[0] = patch->dimensions.unit()[0] * where_i;
    assert(where_i > 0);
    assert(where[0] < patch->orthogonalDimension(side));
    assert(where_i % patch->parallelDimension(side).denom() == 0);
    auto luPatch = std::make_shared<Patch>();
    auto rdPatch = std::make_shared<Patch>();
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
            luPatch->cornersOnEdge(pside).push_back(std::move(corner));
        } else {
            corner->position[0] -= luPatch->parallelDimension(pside);
            rdPatch->cornersOnEdge(pside).push_back(std::move(corner));
        }
    }
    for (auto& corner : patch->cornersOnEdge(~pside)) {
        if (corner->position[0] < luPatch->parallelDimension(~pside)) {
            luPatch->cornersOnEdge(~pside).push_back(std::move(corner));
        } else {
            corner->position[0] -= luPatch->parallelDimension(~pside);
            rdPatch->cornersOnEdge(~pside).push_back(std::move(corner));
        }
    }
    auto luCorners = luPatch->corners(side);
    if (!rdPatch->edge(pside).empty()) {
        auto luSection = rdPatch->edge(pside)[0];
        std::get<0>(luCorners) = std::make_shared<Corner>(
            luSection->leftOrUpPatch, luSection->leftOrUpPosition);
    }
    if (!rdPatch->edge(~pside).empty()) {
        auto rdSection = rdPatch->edge(~pside)[0];
        std::get<1>(luCorners) = std::make_shared<Corner>(
            rdSection->rightOrDownPatch, rdSection->rightOrDownPosition);
    }
    auto rdCorners = rdPatch->corners(~side);
    if (!luPatch->edge(pside).empty()) {
        auto luSection = luPatch->edge(pside).back();
        std::get<0>(rdCorners) = std::make_shared<Corner>(
            luSection->leftOrUpPatch, luSection->leftOrUpPosition);
        std::get<0>(rdCorners)->position[0] -= rdPatch->dimensions.unit()[0];
    }
    if (!luPatch->edge(~pside).empty()) {
        auto rdSection = luPatch->edge(~pside)[0];
        std::get<1>(rdCorners) = std::make_shared<Corner>(
            rdSection->rightOrDownPatch, rdSection->rightOrDownPosition);
        std::get<0>(rdCorners)->position[0] -= rdPatch->dimensions.unit()[0];
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

    auto sharedEdge = std::make_shared<Section>();
    sharedEdge->patch(~side) = luPatch;
    sharedEdge->patch(side) = rdPatch;
    sharedEdge->length[0] = luPatch->parallelDimension(side);
    luPatch->edge(side).push_back(sharedEdge);
    rdPatch->edge(~side).push_back(sharedEdge);
    zoomIn(*luPatch, luFocus);
    zoomIn(*rdPatch, rdFocus);
    return { luPatch, rdPatch };
}

std::weak_ptr<Patch>& Section::patch(Side side)
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

frac1::lview Section::position(Side side)
{
    switch (side) {
        case LEFT:
        case UP:
            return leftOrUpPosition[0];
        case RIGHT:
        case DOWN:
            return rightOrDownPosition[0];
    }
    assert(false);
}

std::tuple<std::shared_ptr<Corner>&, std::shared_ptr<Corner>&> Patch::corners(
    Side side)
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

std::vector<std::shared_ptr<Corner>>& Patch::cornersOnEdge(Side side)
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

std::vector<std::shared_ptr<Section>>& Patch::edge(Side side)
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

frac2::lview Patch::parallelDimension(Side side)
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

frac2::lview Patch::orthogonalDimension(Side side)
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

frac2::lview Patch::parallelPosition(Side side)
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

frac2::lview Patch::orthogonalPosition(Side side)
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

void Patch::prepareBuffer()
{
    this->data.clear();
    this->data.resize(
        (this->dimensions[0].nom() + 2) * (this->dimensions[1].nom() + 2), 1.0);
}

Patch::T Patch::read(u32 x, u32 y) const
{
    u32 index = (this->dimensions[0].nom() + 2) * y + x;
    assert(index < this->data.size());
    return this->data[index];
}

void Patch::write(u32 x, u32 y, Patch::T value)
{
    u32 index = (this->dimensions[0].nom() + 2) * y + x;
    assert(index < this->data.size());
    this->data[index] = value;
}

void Patch::writeEdge(u32 i, Side side, Patch::T value)
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

Patch::T Patch::sample(FracRView pos, FracRView size, Side side) const
{
    u32 sizeU32 = size.denom() >= this->dimensions[0].denom()
                      ? 1
                      : this->fracToLength(size);
    u32 fromX, fromY;
    switch (side) {
        case LEFT:
            fromX = 1;
            fromY = this->fracToLength(pos);
            break;
        case RIGHT:
            fromX = this->dimensions[0].nom() - sizeU32;
            fromY = this->fracToLength(pos);
            break;
        case UP:
            fromX = this->fracToLength(pos);
            fromY = 1;
            break;
        case DOWN:
            fromX = this->fracToLength(pos);
            fromY = this->dimensions[1].nom() - sizeU32;
            break;
        default:
            assert(false);
    }
    u32 sum = 0;
    for (u32 x = fromX; x < fromX + sizeU32; ++x) {
        for (u32 y = fromY; y < fromY + sizeU32; ++y) {
            sum += this->read(x, y);
        }
    }
    return sum;
}

u32 Patch::fracToLength(FracRView frac) const
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

void Patch::synchronizeSection(std::shared_ptr<Section> section, Side side)
{
    for (auto i = frac1(0U); i[0] < section->length[0];
         i[0] += this->dimensions.unit()[0]) {
        u32 coord = 1U + this->fracToLength(i[0] + section->position(~side));
        auto otherPatch = section->patch(side).lock();
        auto otherPos = section->position(side);
        this->writeEdge(coord,
                        side,
                        otherPatch->sample(otherPos + i[0],
                                           this->dimensions.unit()[0],
                                           ~side));
    }
}

void Patch::synchronizeEdges()
{
    for (const auto& section : this->leftEdge) {
        this->synchronizeSection(section, LEFT);
    }
    for (const auto& section : this->rightEdge) {
        this->synchronizeSection(section, RIGHT);
    }
    for (const auto& section : this->upEdge) {
        this->synchronizeSection(section, UP);
    }
    for (const auto& section : this->downEdge) {
        this->synchronizeSection(section, DOWN);
    }
}

void Patch::print() const
{
    for (u32 y = 0; y < this->dimensions[1].nom() + 2; ++y) {
        for (u32 x = 0; x < this->dimensions[0].nom() + 2; ++x) {
            std::cout << this->read(x, y) << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void zoomIn(Patch& patch, u8 level)
{
    patch.dimensions.focus(level);
    patch.prepareBuffer();
}

void DoublePatchGraph::synchronizeEdges()
{
    for (const auto& patch : this->patches) {
        patch->synchronizeEdges();
    }
}

void DoublePatchGraph::print() const
{
    for (const auto& patch : this->patches) {
        patch->print();
    }
}
