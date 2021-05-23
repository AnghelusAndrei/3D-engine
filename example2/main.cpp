#include "SDL2/SDL.h"
#include "SDL2/SDL_image.h"
#include "SDL2/SDL_mutex.h"
#include "SDL2/SDL_ttf.h"
#include <cstdint>
#include <cstdlib>
#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <fstream>
#include <strstream>
#include <string>
#include <algorithm>
#include <list>
#include "engine.h"



//constants

    float FOV = 90;
    const std::uint64_t SCREEN_HEIGHT = 600; 
    const std::uint64_t SCREEN_WIDTH = 800;
    
    
    const float Zfar = 1000;
    const float Znear = 0.1;

    bool KEY_W = false;
    bool KEY_A = false;
    bool KEY_S = false;
    bool KEY_D = false;
    bool KEY_SPACE = false;
    bool KEY_SHIFT = false;
    bool FULLSCREEN = false;
    SDL_atomic_t running;

    vec3 light;

    bool FIRST_FRAME_END = false;
    int START_ANGLEX;
    int START_ANGLEY;
    int START_ANGLEY2;

    float PLAYER_SPEED = 0.00005;

    player_t player = {0,0,0,0,0};

//


std::vector<mesh> InitMesh(){
    std::vector<mesh> mesh_collection;

    mesh Object1, Object2, Object3;
    Object1.LoadFile("./import/world.obj", true);
    Object2.LoadFile("./import/sphere.obj", false);
    Object3.LoadFile("./import/cube.obj", false);
    
    Object1.texture = IMG_Load("./import/texture.jpg");

    mesh_collection.push_back(Object1);
    mesh_collection.push_back(Object2);
    mesh_collection.push_back(Object3);

    return mesh_collection;
}


std::vector<mesh> mesh_collection = InitMesh();


void Input(SDL_Event &e, SDL_Window* window){
    if(KEY_W){
        player.z += cos((FixAng(player.angx)) * M_PI / 180) * PLAYER_SPEED;
        player.x += sin((FixAng(player.angx)) * M_PI / 180) * PLAYER_SPEED;
    }
    if(KEY_A){
        player.z += cos((FixAng(player.angx-90)) * M_PI / 180) * PLAYER_SPEED;
        player.x += sin((FixAng(player.angx-90)) * M_PI / 180) * PLAYER_SPEED;
    }
    if(KEY_S){
        player.z += cos((FixAng(player.angx-180)) * M_PI / 180) * PLAYER_SPEED;
        player.x += sin((FixAng(player.angx-180)) * M_PI / 180) * PLAYER_SPEED;
    }
    if(KEY_D){
        player.z += cos((FixAng(player.angx+90)) * M_PI / 180) * PLAYER_SPEED;
        player.x += sin((FixAng(player.angx+90)) * M_PI / 180) * PLAYER_SPEED;
    }
    if(KEY_SPACE){
        player.y -= PLAYER_SPEED;
    }
    if(KEY_SHIFT){
        player.y += PLAYER_SPEED;
    }
    
        
    if (e.type == SDL_KEYDOWN)
    {

        if(e.key.keysym.sym == SDLK_w){
            KEY_W = true;
        }
        if(e.key.keysym.sym == SDLK_a){
            KEY_A = true;
        }
        if(e.key.keysym.sym == SDLK_s){
            KEY_S = true;
        }
        if(e.key.keysym.sym == SDLK_d){
            KEY_D = true;
        }
        if(e.key.keysym.sym == 32){
            KEY_SPACE = true;
        }
        if(e.key.keysym.sym == 1073742049){
            KEY_SHIFT = true;
        }
            if(e.key.keysym.sym == SDLK_ESCAPE){
                SDL_SetWindowFullscreen(window, 0);
                SDL_SetWindowSize(window, 800, 600);
                SDL_ShowCursor(1);
                FULLSCREEN = false;
            }
            if(e.key.keysym.sym == SDLK_f){
                SDL_SetWindowFullscreen(window, SDL_WINDOW_FULLSCREEN);
                SDL_ShowCursor(0);
                FULLSCREEN = true;
            }
    }
    else if(e.type == SDL_KEYUP){
        if(e.key.keysym.sym == SDLK_w){
            KEY_W = false;
        }
        if(e.key.keysym.sym == SDLK_a){
            KEY_A = false;
        }
        if(e.key.keysym.sym == SDLK_s){
            KEY_S = false;
        }
        if(e.key.keysym.sym == SDLK_d){
            KEY_D = false;
        }
        if(e.key.keysym.sym == SDLK_d){
            KEY_D = false;
        }
        if(e.key.keysym.sym == 32){
            KEY_SPACE = false;
        }
        if(e.key.keysym.sym == 1073742049){
            KEY_SHIFT = false;
            }
        }
    else if (e.type == SDL_MOUSEMOTION)
    {   
        if(FIRST_FRAME_END == false){
            FIRST_FRAME_END = true;
            START_ANGLEX = e.motion.x;
            START_ANGLEY = e.motion.y;

            player.angx = 0;
            player.angy = 0;
        }else{
            player.angx += FixAng(e.motion.x - START_ANGLEX);
            float yAng = e.motion.y - START_ANGLEY;
            player.angy += yAng;
            player.angy = FixY(player.angy);

            START_ANGLEX = e.motion.x;
            START_ANGLEY = e.motion.y;
        }

        if(FULLSCREEN){
            if(e.motion.x > SCREEN_WIDTH-10){
                START_ANGLEX = 10;
                SDL_WarpMouseInWindow(window, 10, e.motion.y);
            }
            if(e.motion.x < 10){
                START_ANGLEX = SCREEN_WIDTH-10;
                SDL_WarpMouseInWindow(window, SCREEN_WIDTH-10, e.motion.y);
            }

            if(e.motion.y > SCREEN_HEIGHT-10){
                START_ANGLEY = 10;
                SDL_WarpMouseInWindow(window, e.motion.x, 10);
            }
            if(e.motion.y < 10){
                START_ANGLEY = SCREEN_HEIGHT-10;
                SDL_WarpMouseInWindow(window, e.motion.x, SCREEN_HEIGHT-10);
            }
        }
    }  
}


static int Rendering(void *data){
    SDL_Renderer *renderer=(SDL_Renderer*)data;

    mat4x4 Projection_Matrix = Matrix_Projection(SCREEN_WIDTH, SCREEN_HEIGHT, FOV, Zfar, Znear);
    float *pDepthBuffer = new float[SCREEN_WIDTH * SCREEN_HEIGHT];


    while (SDL_AtomicGet(&running))
    {
        std::fill_n(pDepthBuffer, SCREEN_HEIGHT*SCREEN_WIDTH, 0.0f);
        float start = (float)SDL_GetTicks();
        Frame(renderer, mesh_collection, Projection_Matrix, player, pDepthBuffer, light, SCREEN_WIDTH, SCREEN_HEIGHT);
        float end = (float)SDL_GetTicks();
        std::cout<<(int)(1000/(end-start))<<" FPS"<<std::endl;
    }

    return 0;
}


int main()
{
    SDL_Window *window;
    SDL_Renderer *renderer;
    window = SDL_CreateWindow("SDL Render Window", 200, 200, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_METAL);
    renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_SOFTWARE);
    SDL_ShowWindow(window);


    SDL_AtomicSet(&running, 1);
    SDL_Thread *renderingThread;

    renderingThread = SDL_CreateThread(Rendering, "rendering", (void *)renderer);

    light.x = -1;
    light.y = -0.75;
    light.z = -0.5;

    mesh_collection[0].position.y = 12;
    mesh_collection[0].position.x = 0;
    mesh_collection[0].position.z = 0;
    mesh_collection[1].position.y = -5;
    mesh_collection[1].position.z = 15;
    mesh_collection[2].position.y = -7;
    mesh_collection[2].position.z = 15;

    mesh_collection[0].size.x = 4;
    mesh_collection[0].size.y = 4;
    mesh_collection[0].size.z = 4;
    mesh_collection[1].size.x = 0.5;
    mesh_collection[1].size.y = 0.5;
    mesh_collection[1].size.z = 0.5;
    mesh_collection[2].size.x = 0.5;
    mesh_collection[2].size.y = 0.5;
    mesh_collection[2].size.z = 0.5;

    SDL_Event e;
    while (SDL_PollEvent(&e) || true)
    {
        float Time = (float)SDL_GetTicks();

        if (e.type == SDL_QUIT)
        {
            SDL_AtomicSet(&running, 0);
            break;
        }

        mesh_collection[1].rotation.x = (float)(Time/50);
        mesh_collection[1].rotation.y = (float)(Time/50);
        mesh_collection[1].rotation.z = (float)(Time/50);
        mesh_collection[2].rotation.x = (float)(Time/50);
        mesh_collection[2].rotation.y = (float)(Time/50);
        mesh_collection[2].rotation.z = (float)(Time/50);

        Input(e, window);

        float end = (float)SDL_GetTicks();
        PLAYER_SPEED = abs(0.0001f+((end-Time)/10));

    }

    int status;
    SDL_WaitThread(renderingThread, &status);

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);

    return 0;
}