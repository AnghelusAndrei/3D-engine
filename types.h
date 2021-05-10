#pragma once

struct color
{
    int r,g,b;
};

struct vec3
{
    double x, y, z;
};

struct player_t
{
    double x, y, z, angx, angy;
};

struct triangle
{
    vec3 p[3];
    color c;
};

struct mat4x4
{
    double m[4][4] = { 0 };
};

struct mesh
{
    std::vector<triangle> tris;

    bool LoadFile(std::string filename){
        std::ifstream f(filename);
        if(!f.is_open()){
            return false;
        }

        std::vector<vec3> verts;

        while(!f.eof()){
            char line[128];
            f.getline(line, 128);

            std::strstream s;
            s<<line;

            char l;
            if(line[0] == 'v'){
                vec3 v;
                s>>l>>v.x>>v.y>>v.z;
                verts.push_back(v);
            }

            if(line[0] == 'f'){
                int f[3];
                s>>l>>f[0]>>f[1]>>f[2];
                triangle trf = {{verts[f[0]-1], verts[f[1]-1], verts[f[2]-1]}};
                tris.push_back(trf);
            }
        }

        return true;
    }
};



double Distance(vec3 vec1, vec3 vec2){
    return sqrt((vec1.x-vec2.x) * (vec1.x-vec2.x) + (vec1.y-vec2.y)*(vec1.y-vec2.y));
}


double degToRad(double ang)
{
    return ang * (M_PI/180);;
}

double radToDeg(double rad)
{
    return rad * (180/M_PI);
}

double FixAng(double ang)
{
    if (ang > 359)
    {
        double mult = floor(ang / 360);
        ang = ang - (360 * mult);
    }
    if (ang < 0)
    {
        double ap = abs(ang);
        double mult = ceil(ap / 360);
        ang = (360 * mult) - ap;
    }
    return ang;
}

uint8_t r(uint8_t n)
{
    return rand() % n;
}

bool Cr(triangle &t1, triangle &t2){
    double z1 = (t1.p[0].z + t1.p[1].z + t1.p[2].z)/3;
    double z2 = (t2.p[0].z + t2.p[1].z + t2.p[2].z)/3;
    return z1 > z2;
}


void MultiplyMatrix(vec3 &i, vec3 &o, mat4x4 &m){
    o.x = i.x * m.m[0][0] + i.y * m.m[1][0] + i.z * m.m[2][0] + m.m[3][0];
    o.y = i.x * m.m[0][1] + i.y * m.m[1][1] + i.z * m.m[2][1] + m.m[3][1];
    o.z = i.x * m.m[0][2] + i.y * m.m[1][2] + i.z * m.m[2][2] + m.m[3][2];
    double w = i.x * m.m[0][3] + i.y * m.m[1][3] + i.z * m.m[2][3] + m.m[3][3];

    if(w!=0){
        o.x /= w;
        o.y /= w;
        o.z /= w;
    }
}


triangle RotateX(triangle tris, double ang){
    mat4x4 Rotated_MatrixX;
    triangle NewTris;

    Rotated_MatrixX.m[0][0] = 1;
    Rotated_MatrixX.m[1][1] = cos(ang);
    Rotated_MatrixX.m[1][2] = sin(ang);
    Rotated_MatrixX.m[2][1] = -sin(ang);
    Rotated_MatrixX.m[2][2] = cos(ang);
    Rotated_MatrixX.m[3][3] = 1;

    MultiplyMatrix(tris.p[0], NewTris.p[0], Rotated_MatrixX);
    MultiplyMatrix(tris.p[1], NewTris.p[1], Rotated_MatrixX);
    MultiplyMatrix(tris.p[2], NewTris.p[2], Rotated_MatrixX);

    return NewTris;
}

triangle RotateY(triangle tris, double ang){
    mat4x4 Rotated_MatrixY;
    triangle NewTris;

    Rotated_MatrixY.m[0][0] = cos(ang);
    Rotated_MatrixY.m[1][1] = 1;
    Rotated_MatrixY.m[0][2] = -sin(ang);
    Rotated_MatrixY.m[2][0] = sin(ang);
    Rotated_MatrixY.m[2][2] = cos(ang);

    MultiplyMatrix(tris.p[0], NewTris.p[0], Rotated_MatrixY);
    MultiplyMatrix(tris.p[1], NewTris.p[1], Rotated_MatrixY);
    MultiplyMatrix(tris.p[2], NewTris.p[2], Rotated_MatrixY);

    return NewTris;
}

triangle RotateZ(triangle tris, double ang){
    mat4x4 Rotated_MatrixZ;
    triangle NewTris;

    Rotated_MatrixZ.m[0][0] = cos(ang);
    Rotated_MatrixZ.m[0][1] = sin(ang);
    Rotated_MatrixZ.m[1][0] = -sin(ang);
    Rotated_MatrixZ.m[1][1] = cos(ang);
    Rotated_MatrixZ.m[2][2] = 1;
    Rotated_MatrixZ.m[3][3] = 1;

    MultiplyMatrix(tris.p[0], NewTris.p[0], Rotated_MatrixZ);
    MultiplyMatrix(tris.p[1], NewTris.p[1], Rotated_MatrixZ);
    MultiplyMatrix(tris.p[2], NewTris.p[2], Rotated_MatrixZ);

    return NewTris;
}


mat4x4 Matrix_CameraPos(player_t &pos, vec3 &target, vec3 &up){

    vec3 newForward = {target.x-pos.x, target.y-pos.y, target.z-pos.z};
    double NFl = sqrt(newForward.x*newForward.x + newForward.y*newForward.y + newForward.z*newForward.z);
    newForward.x /= NFl; newForward.y /= NFl; newForward.z /= NFl;

    double NFUdp = up.x*newForward.x+up.y*newForward.y+up.z*newForward.z;
    vec3 a = {newForward.x * NFUdp, newForward.y * NFUdp, newForward.z * NFUdp};
    vec3 newUp = {up.x - a.x, up.y - a.y, up.z - a.z};
    double NUl = sqrt(newUp.x*newUp.x + newUp.y*newUp.y + newUp.z*newUp.z);
    newUp.x /= NFl; newUp.y /= NFl; newUp.z /= NFl;

    vec3 newRight = {newUp.y * newForward.z - newUp.z * newForward.y, newUp.z * newForward.x - newUp.x * newForward.z, newUp.x * newForward.y - newUp.y * newForward.x};

    mat4x4 matrix;
	matrix.m[0][0] = newRight.x;	matrix.m[0][1] = newRight.y;	matrix.m[0][2] = newRight.z;	matrix.m[0][3] = 0;
	matrix.m[1][0] = newUp.x;		matrix.m[1][1] = newUp.y;		matrix.m[1][2] = newUp.z;		matrix.m[1][3] = 0;
	matrix.m[2][0] = newForward.x;	matrix.m[2][1] = newForward.y;	matrix.m[2][2] = newForward.z;	matrix.m[2][3] = 0;
	matrix.m[3][0] = pos.x;			matrix.m[3][1] = pos.y;			matrix.m[3][2] = pos.z;			matrix.m[3][3] = 1;
	return matrix;

}

mat4x4 Matrix_QuickInverse(mat4x4 &m)
{
	mat4x4 matrix;
	matrix.m[0][0] = m.m[0][0]; matrix.m[0][1] = m.m[1][0]; matrix.m[0][2] = m.m[2][0]; matrix.m[0][3] = 0;
	matrix.m[1][0] = m.m[0][1]; matrix.m[1][1] = m.m[1][1]; matrix.m[1][2] = m.m[2][1]; matrix.m[1][3] = 0;
	matrix.m[2][0] = m.m[0][2]; matrix.m[2][1] = m.m[1][2]; matrix.m[2][2] = m.m[2][2]; matrix.m[2][3] = 0;
	matrix.m[3][0] = -(m.m[3][0] * matrix.m[0][0] + m.m[3][1] * matrix.m[1][0] + m.m[3][2] * matrix.m[2][0]);
	matrix.m[3][1] = -(m.m[3][0] * matrix.m[0][1] + m.m[3][1] * matrix.m[1][1] + m.m[3][2] * matrix.m[2][1]);
	matrix.m[3][2] = -(m.m[3][0] * matrix.m[0][2] + m.m[3][1] * matrix.m[1][2] + m.m[3][2] * matrix.m[2][2]);
	matrix.m[3][3] = 1;
	return matrix;
}



void FillTriangle(	int x1, int y1, int x2, int y2, int x3, int y3, SDL_Renderer *renderer, color c)
{

	if (y2 < y1)
	{
		std::swap(y1, y2);
		std::swap(x1, x2);
	}

	if (y3 < y1)
	{
		std::swap(y1, y3);
		std::swap(x1, x3);
	}

	if (y3 < y2)
	{
		std::swap(y2, y3);
		std::swap(x2, x3);
	}

	int dy1 = y2 - y1;
	int dx1 = x2 - x1;

	int dy2 = y3 - y1;
	int dx2 = x3 - x1;

	double tex_u, tex_v, tex_w;

	double dax_step = 0, dbx_step = 0,
		du1_step = 0, dv1_step = 0,
		du2_step = 0, dv2_step = 0,
		dw1_step=0, dw2_step=0;

	if (dy1) dax_step = dx1 / (double)abs(dy1);
	if (dy2) dbx_step = dx2 / (double)abs(dy2);

	if (dy1)
	{
		for (int i = y1; i <= y2; i++)
		{
			int ax = x1 + (double)(i - y1) * dax_step;
			int bx = x1 + (double)(i - y1) * dbx_step;

			if (ax > bx)
			{
				std::swap(ax, bx);
			}

			SDL_SetRenderDrawColor(renderer, c.r, c.g, c.b, 255);
            SDL_RenderDrawLine(renderer, ax, i, bx, i);
		}
	}

	dy1 = y3 - y2;
	dx1 = x3 - x2;

	if (dy1) dax_step = dx1 / (double)abs(dy1);
	if (dy2) dbx_step = dx2 / (double)abs(dy2);

	du1_step = 0, dv1_step = 0;

	if (dy1)
	{
		for (int i = y2; i <= y3; i++)
		{
			int ax = x2 + (double)(i - y2) * dax_step;
			int bx = x1 + (double)(i - y1) * dbx_step;


			if (ax > bx)
			{
				std::swap(ax, bx);
			}
			SDL_SetRenderDrawColor(renderer, c.r, c.g, c.b, 255);
            SDL_RenderDrawLine(renderer, ax, i, bx, i);
		}	
	}		
}