#include "SDL2/SDL.h"
#include "SDL2/SDL_image.h"
#include "SDL2/SDL_mutex.h"
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
#include "types.h"



//constants

    float FOV = 90;
    std::uint64_t SCREEN_HEIGHT = 600; 
    std::uint64_t SCREEN_WIDTH = 800;
        
    std::uint64_t DESKTOP_HEIGHT = 900; 
    std::uint64_t DESKTOP_WIDTH = 1440;
    
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

    vec3 light = { -1, -0.75, -0.5};

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
    Object1.LoadFile("./import/world.obj", false);
    Object2.LoadFile("./import/sphere.obj", false);
    Object3.LoadFile("./import/cube.obj", false);


    mesh_collection.push_back(Object1);
    mesh_collection.push_back(Object2);
    mesh_collection.push_back(Object3);

    return mesh_collection;
}



void Frame(SDL_Renderer *renderer, std::vector<mesh> mesh_collection, mat4x4 Projection_Matrix, player_t Camera, float *pDepthBuffer){


    vec3 yLookDir = (vec3){0, (float)cos((FixAng(-player.angy)) * M_PI / 180), (float)sin((FixAng(-player.angy)) * M_PI / 180)};

    float l_length = sqrt(light.x*light.x + light.y*light.y + light.z*light.z);
    light.x /= l_length; light.y /= l_length; light.z /= l_length;

    vec3 vLookDir = (vec3){0, 0, 1};
    vLookDir = Vec_RotateX(vLookDir, -player.angy);
    vLookDir = Vec_RotateY(vLookDir, player.angx);
    vec3 vUp = { 0,1,0 };
    vec3 vTarget = {Camera.x + vLookDir.x, Camera.y + vLookDir.y, Camera.z + vLookDir.z};

    mat4x4 matCamera = Matrix_CameraPos(Camera, vTarget, vUp);
    mat4x4 matView = Matrix_QuickInverse(matCamera);



    for(int obj = 0; obj < mesh_collection.size(); obj++){
        mesh Mesh = mesh_collection[obj];

        for(int i = 0; i < Mesh.tris.size(); i++){
            triangle transformation, projection, view;

            transformation = Mesh.tris[i];
            transformation = RotateX(transformation, 180);

            transformation = RotateY(transformation, Mesh.rotation.y);
            transformation = RotateZ(transformation, Mesh.rotation.z);
            transformation = RotateX(transformation, Mesh.rotation.x);

            transformation.p[0] = Scale(transformation.p[0], Mesh.size);
            transformation.p[1] = Scale(transformation.p[1], Mesh.size);
            transformation.p[2] = Scale(transformation.p[2], Mesh.size);

            transformation.p[0] = Translate(transformation.p[0], Mesh.position);
            transformation.p[1] = Translate(transformation.p[1], Mesh.position);
            transformation.p[2] = Translate(transformation.p[2], Mesh.position);

            vec3 normal, line1, line2;
            line1.x = transformation.p[1].x - transformation.p[0].x;
            line1.y = transformation.p[1].y - transformation.p[0].y;
            line1.z = transformation.p[1].z - transformation.p[0].z;

            line2.x = transformation.p[2].x - transformation.p[0].x;
            line2.y = transformation.p[2].y - transformation.p[0].y;
            line2.z = transformation.p[2].z - transformation.p[0].z;

            normal.x = line1.y * line2.z - line1.z * line2.y;
            normal.y = line1.z * line2.x - line1.x * line2.z;
            normal.z = line1.x * line2.y - line1.y * line2.x;

            float l = sqrt(normal.x*normal.x + normal.y*normal.y + normal.z*normal.z);
            normal.x /= l; normal.y /= l; normal.z /= l;

            if(normal.x * (transformation.p[0].x - Camera.x) + normal.y * (transformation.p[0].y - Camera.y) + normal.z * (transformation.p[0].z - Camera.z) < 0)
            {

                MultiplyMatrix(transformation.p[0], view.p[0], matView);
                MultiplyMatrix(transformation.p[1], view.p[1], matView);
                MultiplyMatrix(transformation.p[2], view.p[2], matView);

                float light_dp = normal.x*light.x+normal.y*light.y+normal.z*light.z;

                int shade = (int)(light_dp*245+10);
                if(shade < 10 || shade > 255){
                    shade = 10;
                }


                int nClippedTriangles = 0;
                triangle clipped[2];
                nClippedTriangles = Triangle_ClipAgainstPlane((vec3){ 0, 0, 0.1 }, (vec3){ 0, 0, 1 }, view, clipped[0], clipped[1]);

                for (int n = 0; n < nClippedTriangles; n++)
                {
                    MultiplyMatrix(clipped[n].p[0], projection.p[0], Projection_Matrix);
                    MultiplyMatrix(clipped[n].p[1], projection.p[1], Projection_Matrix);
                    MultiplyMatrix(clipped[n].p[2], projection.p[2], Projection_Matrix);

                    projection.t[0] = clipped[n].t[0];
                    projection.t[1] = clipped[n].t[1];
                    projection.t[2] = clipped[n].t[2];

                    projection.c = clipped[n].c;

                    projection.p[0].x *= 1;
                    projection.p[1].x *= 1;
                    projection.p[2].x *= 1;
                    projection.p[0].y *= 1;
                    projection.p[1].y *= 1;
                    projection.p[2].y *= 1;

                    projection.p[0].x += 1; projection.p[0].y += 1;
                    projection.p[1].x += 1; projection.p[1].y += 1;
                    projection.p[2].x += 1; projection.p[2].y += 1;

                    projection.p[0].x *= SCREEN_HEIGHT*0.5; projection.p[0].y *= SCREEN_HEIGHT*0.5;
                    projection.p[1].x *= SCREEN_HEIGHT*0.5; projection.p[1].y *= SCREEN_HEIGHT*0.5;
                    projection.p[2].x *= SCREEN_HEIGHT*0.5; projection.p[2].y *= SCREEN_HEIGHT*0.5;

                    projection.p[0].x += (SCREEN_WIDTH-SCREEN_HEIGHT)/2;
                    projection.p[1].x += (SCREEN_WIDTH-SCREEN_HEIGHT)/2;
                    projection.p[2].x += (SCREEN_WIDTH-SCREEN_HEIGHT)/2;

                    projection.c = (color){shade, shade, shade};


                    projection.t[0].x = projection.t[0].x / GetW(clipped[n].p[0], Projection_Matrix);
                    projection.t[1].x = projection.t[1].x / GetW(clipped[n].p[1], Projection_Matrix);
                    projection.t[2].x = projection.t[2].x / GetW(clipped[n].p[2], Projection_Matrix);

                    projection.t[0].y = projection.t[0].y / GetW(clipped[n].p[0], Projection_Matrix);
                    projection.t[1].y = projection.t[1].y / GetW(clipped[n].p[1], Projection_Matrix);
                    projection.t[2].y = projection.t[2].y / GetW(clipped[n].p[2], Projection_Matrix);

                    projection.t[0].z = 1/GetW(clipped[n].p[0], Projection_Matrix);
                    projection.t[1].z = 1/GetW(clipped[n].p[1], Projection_Matrix);
                    projection.t[2].z = 1/GetW(clipped[n].p[2], Projection_Matrix);

                    triangle clipped[2];
                    std::list<triangle> listTriangles;

                    listTriangles.push_back(projection);
                    int nNewTriangles = 1;

                    for (int p = 0; p < 4; p++)
                    {
                        int nTrisToAdd = 0;
                        while (nNewTriangles > 0)
                        {
                            triangle test = listTriangles.front();
                            listTriangles.pop_front();
                            nNewTriangles--;

                            switch (p)
                            {
                            case 0:	nTrisToAdd = Triangle_ClipAgainstPlane((vec3){ 0, 0, 0 }, (vec3){ 0, 1, 0 }, test, clipped[0], clipped[1]); break;
                            case 1:	nTrisToAdd = Triangle_ClipAgainstPlane((vec3){ 0, (float)SCREEN_HEIGHT - 1, 0 }, (vec3){ 0, -1, 0 }, test, clipped[0], clipped[1]); break;
                            case 2:	nTrisToAdd = Triangle_ClipAgainstPlane((vec3){ 0, 0, 0 }, (vec3){ 1, 0, 0 }, test, clipped[0], clipped[1]); break;
                            case 3:	nTrisToAdd = Triangle_ClipAgainstPlane((vec3){ (float)SCREEN_WIDTH - 1, 0, 0 }, (vec3){ -1, 0, 0 }, test, clipped[0], clipped[1]); break;
                            }
                            for (int w = 0; w < nTrisToAdd; w++)
                                listTriangles.push_back(clipped[w]);
                        }
                        nNewTriangles = listTriangles.size();
                    }


                    for (std::list<triangle>::iterator t = listTriangles.begin(); t != listTriangles.end(); t++)
                    {
                        SDL_SetRenderDrawColor(renderer, t->c.r, t->c.g, t->c.b, 255);
                        FillTriangle(t->p[0].x, t->p[0].y, t->p[1].x, t->p[1].y, t->p[2].x, t->p[2].y, t->t[0].z, t->t[1].z, t->t[2].z, renderer, pDepthBuffer, SCREEN_WIDTH, SCREEN_HEIGHT);
                    }
                }
            }
        }
    }
}


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
                SCREEN_WIDTH = 800;
                SCREEN_HEIGHT = 600;
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

    std::vector<mesh> mesh_collection = InitMesh();

    mat4x4 Projection_Matrix;

    Projection_Matrix.m[0][0] = (SCREEN_WIDTH/SCREEN_HEIGHT) * (1/tan(degToRad(FOV/2)));
    Projection_Matrix.m[1][1] = (1/tan(degToRad(FOV/2)));
    Projection_Matrix.m[2][2] = Zfar / (Zfar-Znear);
    Projection_Matrix.m[3][2] = (-Zfar * Znear) / (Zfar - Znear);
    Projection_Matrix.m[2][3] = 1;
    Projection_Matrix.m[3][3] = 0;

    float *pDepthBuffer = new float[SCREEN_WIDTH * SCREEN_HEIGHT];


    mesh_collection[0].position = (vec3){0.0f, 5.0f, 0.0f};
    mesh_collection[1].position = (vec3){0.0f, -5.0f, 15.0f};
    mesh_collection[2].position = (vec3){0.0f, -7.0f, 15.0f};

    mesh_collection[0].size = (vec3){2, 3, 2};
    mesh_collection[1].size = (vec3){0.5f, 0.5f, 0.5f};
    mesh_collection[2].size = (vec3){0.5f, 0.5f, 0.5f};


    while (SDL_AtomicGet(&running))
    {
        float Time = (float)SDL_GetTicks();
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);

        mesh_collection[1].rotation = (vec3){(float)(Time/50), (float)(Time/50), (float)(Time/50)};
        mesh_collection[2].rotation = (vec3){(float)(Time/10), (float)(Time/10), (float)(Time/10)};


        std::fill_n(pDepthBuffer, SCREEN_HEIGHT*SCREEN_WIDTH, 0.0f);
        Frame(renderer, mesh_collection, Projection_Matrix, player, pDepthBuffer);

        SDL_RenderPresent(renderer);
        
        float end = (float)SDL_GetTicks();
        PLAYER_SPEED = 0.0001f-((end-Time)/1000000);
    }

    return 0;
}


int main()
{
    SDL_Window *window;
    SDL_Renderer *renderer;
    window = SDL_CreateWindow("SDL Render Window", (int)((DESKTOP_WIDTH-SCREEN_WIDTH)/2), (int)((DESKTOP_HEIGHT-SCREEN_HEIGHT)/2), SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_METAL);
    renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_SOFTWARE);
    SDL_ShowWindow(window);


    SDL_AtomicSet(&running, 1);
    SDL_Thread *renderingThread;

    renderingThread = SDL_CreateThread(Rendering, "rendering", (void *)renderer);

    

    SDL_Event e;
    while (SDL_PollEvent(&e) || true)
    {
        if (e.type == SDL_QUIT)
        {
            SDL_AtomicSet(&running, 0);
            break;
        }

        Input(e, window);

    }

    int status;
    SDL_WaitThread(renderingThread, &status);

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);

    return 0;
}