#include "Frac.hpp"
#include "PatchGraph.hpp"
#include "SDL.h"

int frac_to_pixel(FracRView frac)
{
    return (64 / frac.denom()) * frac.nom();
}

void draw_rect(SDL_Renderer* rnd, int x, int y, int w, int h)
{
    SDL_RenderDrawLine(rnd, x, y, x + w, y);
    SDL_RenderDrawLine(rnd, x + w, y, x + w, y + h);
    SDL_RenderDrawLine(rnd, x + w, y + h, x, y + h);
    SDL_RenderDrawLine(rnd, x, y + h, x, y);
}

void draw_block(SDL_Renderer* rnd, const Patch& block)
{
    size_t pos_coeff = 64 >> block.position[0].denominator;
    int x = pos_coeff * block.position[0].nominator;
    int y = pos_coeff * block.position[1].nominator;
    int rect_size = 64 >> block.dimensions[0].denominator;
    for (size_t i = 0; i < block.dimensions[0].nominator; ++i) {
        for (size_t j = 0; j < block.dimensions[1].nominator; ++j) {
            draw_rect(rnd,
                      x + i * rect_size,
                      y + j * rect_size,
                      rect_size,
                      rect_size);
        }
    }
}

void splitAndFocus(std::vector<std::shared_ptr<Patch>>& blocks,
                   size_t which,
                   size_t where,
                   bool vertical)
{
    auto bl = blocks[which];
    blocks.erase(blocks.begin() + which);
    auto nbl = splitAndFocus(std::move(bl), where, vertical, 0, 0);
    blocks.push_back(nbl.first);
    blocks.push_back(nbl.second);
}

void highlightCell(SDL_Renderer* rnd,
                   const std::vector<std::shared_ptr<Patch>>& patches,
                   int x,
                   int y)
{
    for (const auto& patch : patches) {
        int rx = frac_to_pixel(patch->position[0]);
        int ry = frac_to_pixel(patch->position[1]);
        int rw = frac_to_pixel(patch->dimensions[0]);
        int rh = frac_to_pixel(patch->dimensions[1]);
        if (rx <= x && x < rx + rw && ry <= y && y < ry + rh) {
            int ca = frac_to_pixel(patch->dimensions.unit()[0]);
            int cx = ((x - rx) / ca) * ca;
            int cy = ((y - ry) / ca) * ca;
            SDL_Rect rect = { rx + cx, ry + cy, ca, ca };
            SDL_SetRenderDrawColor(rnd, 255, 0, 0, SDL_ALPHA_OPAQUE);
            SDL_RenderFillRect(rnd, &rect);
        }
    }
}

int main()
{
    std::vector<std::shared_ptr<Patch>> blocks;
    auto bl = std::make_shared<Patch>();
    bl->dimensions = frac2(7, 7);

    blocks.push_back(bl);
    bl->prepareBuffer();
    blocks[0]->print();
    std::cout << "---" << std::endl;
    splitAndFocus(blocks, 0, 3, true);
    blocks[0]->print();
    blocks[1]->print();
    std::cout << "---" << std::endl;
    blocks[0]->synchronizeEdges();
    blocks[1]->synchronizeEdges();
    blocks[0]->print();
    blocks[1]->print();
    std::cout << "---" << std::endl;
    splitAndFocus(blocks, 0, 2, false);
    blocks[0]->print();
    blocks[1]->print();
    blocks[2]->print();
    std::cout << "---" << std::endl;
    blocks[0]->synchronizeEdges();
    blocks[1]->synchronizeEdges();
    blocks[2]->synchronizeEdges();
    blocks[0]->print();
    blocks[1]->print();
    blocks[2]->print();
    std::cout << "---" << std::endl;
    zoomIn(*blocks[0], 1);
    zoomIn(*blocks[1], 2);
    blocks[0]->print();
    blocks[1]->print();
    blocks[2]->print();
    std::cout << "---" << std::endl;
    blocks[0]->synchronizeEdges();
    blocks[1]->synchronizeEdges();
    blocks[2]->synchronizeEdges();
    blocks[0]->print();
    blocks[1]->print();
    blocks[2]->print();
    std::cout << "---" << std::endl;
    int hlx, hly;
    if (SDL_Init(SDL_INIT_VIDEO) == 0) {
        SDL_Window* window = NULL;
        SDL_Renderer* renderer = NULL;

        if (SDL_CreateWindowAndRenderer(640, 480, 0, &window, &renderer) == 0) {
            SDL_bool done = SDL_FALSE;

            while (!done) {
                SDL_Event event;

                SDL_SetRenderDrawColor(renderer, 0, 0, 0, SDL_ALPHA_OPAQUE);
                SDL_RenderClear(renderer);

                uint8_t color = 255;
                for (auto block : blocks) {
                    SDL_SetRenderDrawColor(renderer,
                                           color,
                                           255 - color,
                                           127 + color,
                                           SDL_ALPHA_OPAQUE);
                    draw_block(renderer, *block);
                    color -= 64;
                }
                highlightCell(renderer, blocks, hlx, hly);
                SDL_RenderPresent(renderer);

                while (SDL_PollEvent(&event)) {
                    if (event.type == SDL_QUIT) {
                        done = SDL_TRUE;
                    } else if (event.type == SDL_MOUSEMOTION) {
                        hlx = event.motion.x;
                        hly = event.motion.y;
                    }
                }
            }
        }

        if (renderer) {
            SDL_DestroyRenderer(renderer);
        }
        if (window) {
            SDL_DestroyWindow(window);
        }
    }
    SDL_Quit();
    return 0;
}
