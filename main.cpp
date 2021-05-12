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

    const double FOV = 90;
    const std::uint64_t SCREEN_HEIGHT = 600; 
    const std::uint64_t SCREEN_WIDTH = 800;
    const double Zfar = 1000;
    const double Znear = 0.1;

    bool KEY_W = false;
    bool KEY_A = false;
    bool KEY_S = false;
    bool KEY_D = false;
    bool KEY_SPACE = false;
    bool KEY_SHIFT = false;

    bool FIRST_FRAME_END = false;
    int START_ANGLEX;
    int START_ANGLEY;

    const double PLAYER_SPEED = 0.1;

    player_t player = {0,0,0,0,0};

//

std::vector<mesh> InitMesh(){
    std::vector<mesh> mesh_collection;

    mesh Object1, Object2, Object3;
    Object1.LoadFile("./import/world.obj");
    Object2.LoadFile("./import/sphere.obj");
    Object3.LoadFile("./import/cube.obj");


    mesh_collection.push_back(Object1);
    mesh_collection.push_back(Object2);
    mesh_collection.push_back(Object3);

    return mesh_collection;
}

void Frame(SDL_Renderer *renderer, std::vector<mesh> mesh_collection, mat4x4 Projection_Matrix, player_t Camera){
    std::vector<triangle> tris_collection;

    double TimeVar = ((double)SDL_GetTicks());


    for(int obj = 0; obj < mesh_collection.size(); obj++){
        mesh Mesh = mesh_collection[obj];

        for(int i = 0; i < Mesh.tris.size(); i++){

            vec3 tVec, sVec;
            triangle transformation = Mesh.tris[i];

            switch (obj)
            {
                case 0:
                //world

                    transformation = RotateX(transformation, 180);

                    tVec = (vec3){0, 0, 0};
                    sVec = (vec3){2, 3, 2};

                    transformation.p[0] = Scale(transformation.p[0], sVec);
                    transformation.p[1] = Scale(transformation.p[1], sVec);
                    transformation.p[2] = Scale(transformation.p[2], sVec);

                    transformation.p[0] = Translate(transformation.p[0], tVec);
                    transformation.p[1] = Translate(transformation.p[1], tVec);
                    transformation.p[2] = Translate(transformation.p[2], tVec);

                    break;
                case 1:
                //cube

                    transformation = RotateX(Mesh.tris[i], 180);

                    transformation = RotateY(transformation, TimeVar/10);
                    transformation = RotateZ(transformation, TimeVar/10);
                    transformation = RotateX(transformation, TimeVar/10);

                    tVec = (vec3){0, -5, 15};
                    sVec = (vec3){0.5, 0.5, 0.5};

                    transformation.p[0] = Scale(transformation.p[0], sVec);
                    transformation.p[1] = Scale(transformation.p[1], sVec);
                    transformation.p[2] = Scale(transformation.p[2], sVec);

                    transformation.p[0] = Translate(transformation.p[0], tVec);
                    transformation.p[1] = Translate(transformation.p[1], tVec);
                    transformation.p[2] = Translate(transformation.p[2], tVec);

                    break;
                case 2:
                //sphere

                    transformation = RotateX(Mesh.tris[i], 180);

                    transformation = RotateY(transformation, TimeVar/10);
                    transformation = RotateZ(transformation, TimeVar/10);
                    transformation = RotateX(transformation, TimeVar/50);

                    tVec = (vec3){0, -7, 15};
                    sVec = (vec3){0.5, 0.5, 0.5};

                    transformation.p[0] = Scale(transformation.p[0], sVec);
                    transformation.p[1] = Scale(transformation.p[1], sVec);
                    transformation.p[2] = Scale(transformation.p[2], sVec);

                    transformation.p[0] = Translate(transformation.p[0], tVec);
                    transformation.p[1] = Translate(transformation.p[1], tVec);
                    transformation.p[2] = Translate(transformation.p[2], tVec);

                    break;
                default:

                    break;
            }

            tris_collection.push_back(transformation);
        }
    }


    
    double yAngVec = player.angy*2;
    double yVec1 = sin((FixAng(yAngVec)) * M_PI / 180);
    double yVec2 = sin((FixAng(yAngVec+90)) * M_PI / 180);

    vec3 vLookDir = { sin((FixAng(player.angx)) * M_PI / 180) ,yVec1, cos((FixAng(player.angx)) * M_PI / 180) };
    vec3 vUp = { 0 ,1, 0 };
    vec3 vTarget = {Camera.x + vLookDir.x, Camera.y + vLookDir.y, Camera.z + vLookDir.z};

    mat4x4 matCamera = Matrix_CameraPos(Camera, vTarget, vUp);
    mat4x4 matView = Matrix_QuickInverse(matCamera);
        
    std::vector<triangle> triangleSorter;

    for(int i = 0; i < tris_collection.size(); i++){
        triangle projection, translation, view;


        translation = tris_collection[i];

        vec3 normal, line1, line2;
        line1.x = translation.p[1].x - translation.p[0].x;
        line1.y = translation.p[1].y - translation.p[0].y;
        line1.z = translation.p[1].z - translation.p[0].z;

        line2.x = translation.p[2].x - translation.p[0].x;
        line2.y = translation.p[2].y - translation.p[0].y;
        line2.z = translation.p[2].z - translation.p[0].z;

        normal.x = line1.y * line2.z - line1.z * line2.y;
        normal.y = line1.z * line2.x - line1.x * line2.z;
        normal.z = line1.x * line2.y - line1.y * line2.x;

		double l = sqrt(normal.x*normal.x + normal.y*normal.y + normal.z*normal.z);
		normal.x /= l; normal.y /= l; normal.z /= l;

        if(normal.x * (translation.p[0].x - Camera.x) + normal.y * (translation.p[0].y - Camera.y) + normal.z * (translation.p[0].z - Camera.z) < 0)
        {
            MultiplyMatrix(translation.p[0], view.p[0], matView);
            MultiplyMatrix(translation.p[1], view.p[1], matView);
            MultiplyMatrix(translation.p[2], view.p[2], matView);

            vec3 light = { -1, -0.75, -0.5};

			double l_length = sqrt(light.x*light.x + light.y*light.y + light.z*light.z);
			light.x /= l_length; light.y /= l_length; light.z /= l_length;
            double light_dp = normal.x*light.x+normal.y*light.y+normal.z*light.z;

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

				triangleSorter.push_back(projection);
			}
        }
    }


    std::sort(triangleSorter.begin(), triangleSorter.end(), Cr);

    for(int i = 0; i < triangleSorter.size(); i++){
		triangle clipped[2];
		std::list<triangle> listTriangles;

		listTriangles.push_back(triangleSorter[i]);
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
				case 1:	nTrisToAdd = Triangle_ClipAgainstPlane((vec3){ 0, SCREEN_HEIGHT - 1, 0 }, (vec3){ 0, -1, 0 }, test, clipped[0], clipped[1]); break;
				case 2:	nTrisToAdd = Triangle_ClipAgainstPlane((vec3){ 0, 0, 0 }, (vec3){ 1, 0, 0 }, test, clipped[0], clipped[1]); break;
				case 3:	nTrisToAdd = Triangle_ClipAgainstPlane((vec3){ SCREEN_WIDTH - 1, 0, 0 }, (vec3){ -1, 0, 0 }, test, clipped[0], clipped[1]); break;
				}
				for (int w = 0; w < nTrisToAdd; w++)
					listTriangles.push_back(clipped[w]);
			}
			nNewTriangles = listTriangles.size();
		}


		for (std::list<triangle>::iterator t = listTriangles.begin(); t != listTriangles.end(); t++)
		{
            SDL_SetRenderDrawColor(renderer, t->c.r, t->c.g, t->c.b, 255);
			FillTriangle(t->p[0].x, t->p[0].y, t->p[1].x, t->p[1].y, t->p[2].x, t->p[2].y, renderer);

            SDL_SetRenderDrawColor(renderer, t->c.r, t->c.g, t->c.b, 255);
            SDL_RenderDrawLine(renderer, static_cast<int>(floor(t->p[0].x)), static_cast<int>(floor(t->p[0].y)),  static_cast<int>(floor(t->p[1].x)), static_cast<int>(floor(t->p[1].y)));
            SDL_RenderDrawLine(renderer, static_cast<int>(floor(t->p[2].x)), static_cast<int>(floor(t->p[2].y)),  static_cast<int>(floor(t->p[1].x)), static_cast<int>(floor(t->p[1].y)));
            SDL_RenderDrawLine(renderer, static_cast<int>(floor(t->p[2].x)), static_cast<int>(floor(t->p[2].y)),  static_cast<int>(floor(t->p[0].x)), static_cast<int>(floor(t->p[0].y)));
            
		}
    }
    
}

void Input(SDL_Event &e){
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
            player.angx = FixAng(e.motion.x - START_ANGLEX);
            /*if(player.angx*2 > 90){player.angy = (e.motion.y - START_ANGLEY - (player.angy/2-90));}
            else if(player.angx*2 < -90){player.angy = (e.motion.y - START_ANGLEY - (player.angy/2+90));}
            else{player.angy = (e.motion.y - START_ANGLEY);}*/
        }
    }  
}

int main()
{
    SDL_Window *window;
    SDL_Renderer *renderer;
    window = SDL_CreateWindow("SDL Render Window", 100, 100, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_METAL);
    renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_SOFTWARE);
    SDL_ShowWindow(window);


    std::vector<mesh> mesh_collection = InitMesh();

    mat4x4 Projection_Matrix;

    Projection_Matrix.m[0][0] = (SCREEN_WIDTH/SCREEN_HEIGHT) * (1/tan(degToRad(FOV/2)));
    Projection_Matrix.m[1][1] = (1/tan(degToRad(FOV/2)));
    Projection_Matrix.m[2][2] = Zfar / (Zfar-Znear);
    Projection_Matrix.m[3][2] = (-Zfar * Znear) / (Zfar - Znear);
    Projection_Matrix.m[2][3] = 1;
    Projection_Matrix.m[3][3] = 0;

    

    SDL_Event e;
    while (SDL_PollEvent(&e) || true)
    {
        if (e.type == SDL_QUIT)
            break;

        Input(e);

        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);
        

        Frame(renderer, mesh_collection, Projection_Matrix, player);

        SDL_RenderPresent(renderer);
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);

    return 0;
}


//g++ engine.cpp main.cpp -lSDL2 -lSDL2_image --std=c++17
//./a.out