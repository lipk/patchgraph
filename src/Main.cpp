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

void draw_block(SDL_Renderer* rnd, const Patch<double>& block)
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

void highlightCell(SDL_Renderer* rnd,
                   const std::vector<std::shared_ptr<Patch<double>>>& patches,
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

    auto graph = createPatchGraph<double>(7U, 7U, downsample, upsample);

    for (size_t y = 0; y < 7U; ++y) {
        for (size_t x = 0; x < 7U; ++x) {
            graph.write(x, y, 1);
        }
    }

    graph.apply(1,
                [](const DataReader<double>& reader,
                   bool& focusHere,
                   bool& defocusHere,
                   u32 x,
                   u32 y) {
                    focusHere = (x == 2 && y == 2) || (x == 6 && y == 6);
                    return reader.read(x, y);
                });
    graph.apply(1,
                [](const DataReader<double>& reader,
                   bool& focusHere,
                   bool& defocusHere,
                   u32 x,
                   u32 y) {
                    defocusHere = true;
                    return reader.read(x, y);
                });

    //    int i = 0;
    //    for (auto patch : graph.getSource()) {
    //        if (patch->dimensions.denom() > 1) {
    //            patch->defocus(1, downsample);
    //        }
    //        std::cout << i++ << ": " << patch << ' ' <<
    //        patch->dimensions.nom(0)
    //                  << ' ' << patch->dimensions.nom(1) << std::endl;
    //    }
    //    std::cout << graph.getSource()[2]->canMerge(Side::Up) << std::endl;

    //    auto other = graph.getSource()[2]->canMerge(Side::Up);
    //    auto merged = graph.getSource()[2]->merge(Side::Up);

    //    graph.getSource().erase(graph.getSource().begin() + 4);
    //    graph.getSource().erase(graph.getSource().begin() + 2);
    //    graph.getSource().push_back(merged);

    // graph.print();
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
                for (auto block : graph.getSource()) {
                    SDL_SetRenderDrawColor(renderer,
                                           color,
                                           255 - color,
                                           127 + color,
                                           SDL_ALPHA_OPAQUE);
                    draw_block(renderer, *block);
                    color -= 64;
                }
                highlightCell(renderer, graph.getSource(), hlx, hly);
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
