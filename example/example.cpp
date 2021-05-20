#include "SDL2/SDL.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <strstream>
#include <string>
#include <algorithm>
#include <list>
#include "engine.h"



int main()
{
    const int SCREEN_HEIGHT = 600; 
    const int SCREEN_WIDTH = 800;
    SDL_Window *window;
    SDL_Renderer *renderer;
    SDL_atomic_t running;
    window = SDL_CreateWindow("SFDL Example", 100, 100, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_METAL);
    renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_SOFTWARE);
    SDL_ShowWindow(window);

    float FOV = 90;
    float Zfar = 1000;
    float Znear = 0.1;
    mat4x4 Projection_Matrix = Matrix_Projection(SCREEN_WIDTH, SCREEN_HEIGHT, FOV, Zfar, Znear);
    float *pDepthBuffer = new float[SCREEN_WIDTH * SCREEN_HEIGHT];

    mesh Cube;
    Cube.LoadFile("./import/sphere.obj");

    player_t Camera = {0,0,0,0,0};

    vec3 light;
    light.x = -1;
    light.y = -0.75;
    light.z = -0.5;

    Cube.position.x = 0;
    Cube.position.y = 0;
    Cube.position.z = 5;

    Cube.size.x = 1;
    Cube.size.y = 1;
    Cube.size.z = 1;

    


    SDL_Event e;
    while (SDL_PollEvent(&e) || true)
    {
        if (e.type == SDL_QUIT)
        {
            SDL_AtomicSet(&running, 0);
            break;
        }

        Cube.rotation.x = (float)(SDL_GetTicks()/50);
        Cube.rotation.y = (float)(SDL_GetTicks()/50);
        Cube.rotation.z = (float)(SDL_GetTicks()/50);

        std::vector<mesh> mesh_collection;
        mesh_collection.push_back(Cube);
        Frame(renderer, mesh_collection, Projection_Matrix, Camera, pDepthBuffer, light, SCREEN_WIDTH, SCREEN_HEIGHT);
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);

    return 0;
}