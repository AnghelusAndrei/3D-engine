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


vec3 Vector_Add(vec3 &v1, vec3 &v2)
{
	return { v1.x + v2.x, v1.y + v2.y, v1.z + v2.z };
}

vec3 Vector_Sub(vec3 &v1, vec3 &v2)
{
	return { v1.x - v2.x, v1.y - v2.y, v1.z - v2.z };
}

vec3 Vector_Mul(vec3 &v1, double k)
{
	return { v1.x * k, v1.y * k, v1.z * k };
}

vec3 Vector_Div(vec3 &v1, double k)
{
	return { v1.x / k, v1.y / k, v1.z / k };
}

double Vector_DotProduct(vec3 &v1, vec3 &v2)
{
	return v1.x*v2.x + v1.y*v2.y + v1.z * v2.z;
}

double Vector_Length(vec3 &v)
{
	return sqrt(Vector_DotProduct(v, v));
}

vec3 Vector_Normalise(vec3 &v)
{
	double l = Vector_Length(v);
	return { v.x / l, v.y / l, v.z / l };
}


vec3 Vector_IntersectPlane(vec3 &plane_p, vec3 &plane_n, vec3 &lineStart, vec3 &lineEnd)
{
	plane_n = Vector_Normalise(plane_n);
	double plane_d = -Vector_DotProduct(plane_n, plane_p);
	double ad = Vector_DotProduct(lineStart, plane_n);
	double bd = Vector_DotProduct(lineEnd, plane_n);
	double t = (-plane_d - ad) / (bd - ad);
	vec3 lineStartToEnd = Vector_Sub(lineEnd, lineStart);
	vec3 lineToIntersect = Vector_Mul(lineStartToEnd, t);
	return Vector_Add(lineStart, lineToIntersect);
}


double dist(vec3 &p, vec3 &plane_n, vec3 &plane_p){
	vec3 n = Vector_Normalise(p);
	return (plane_n.x * p.x + plane_n.y * p.y + plane_n.z * p.z - Vector_DotProduct(plane_n, plane_p));
}

int Triangle_ClipAgainstPlane(vec3 plane_p, vec3 plane_n, triangle &in_tri, triangle &out_tri1, triangle &out_tri2)
{
	plane_n = Vector_Normalise(plane_n);

	vec3* inside_points[3];  int nInsidePointCount = 0;
	vec3* outside_points[3]; int nOutsidePointCount = 0;
	double d0 = dist(in_tri.p[0], plane_n, plane_p);
	double d1 = dist(in_tri.p[1], plane_n, plane_p);
	double d2 = dist(in_tri.p[2], plane_n, plane_p);

	if (d0 >= 0) { inside_points[nInsidePointCount++] = &in_tri.p[0]; }
	else { outside_points[nOutsidePointCount++] = &in_tri.p[0]; }
	if (d1 >= 0) { inside_points[nInsidePointCount++] = &in_tri.p[1]; }
	else { outside_points[nOutsidePointCount++] = &in_tri.p[1]; }
	if (d2 >= 0) { inside_points[nInsidePointCount++] = &in_tri.p[2]; }
	else { outside_points[nOutsidePointCount++] = &in_tri.p[2]; }


	if (nInsidePointCount == 0)
	{

		return 0;
	}

	if (nInsidePointCount == 3)
	{
		out_tri1 = in_tri;

		return 1;
	}

	if (nInsidePointCount == 1 && nOutsidePointCount == 2)
	{
		out_tri1.c =  in_tri.c;

		out_tri1.p[0] = *inside_points[0];

		out_tri1.p[1] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[0]);
		out_tri1.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[1]);

		return 1;
	}

	if (nInsidePointCount == 2 && nOutsidePointCount == 1)
	{
		out_tri1.c =  in_tri.c;

		out_tri2.c =  in_tri.c;

		out_tri1.p[0] = *inside_points[0];
		out_tri1.p[1] = *inside_points[1];
		out_tri1.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[0]);

		out_tri2.p[0] = *inside_points[1];
		out_tri2.p[1] = out_tri1.p[2];
		out_tri2.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[1], *outside_points[0]);
		return 2;
	}

    return 0;
}



vec3 Translate(vec3 &v1, vec3 v2){
    vec3 vec;
    mat4x4 matrix;

	matrix.m[0][0] = 1;
	matrix.m[1][1] = 1;
	matrix.m[2][2] = 1;
	matrix.m[3][3] = 1;
	matrix.m[3][0] = v2.x;
	matrix.m[3][1] = v2.y;
	matrix.m[3][2] = v2.z;

    MultiplyMatrix(v1, vec, matrix);

    return vec;
}

vec3 Scale(vec3 &v1, vec3 v2){
    vec3 vec;
    mat4x4 matrix;

	matrix.m[0][0] = v2.x;
	matrix.m[1][1] = v2.y;
	matrix.m[2][2] = v2.z;
	matrix.m[3][3] = 1;

    MultiplyMatrix(v1, vec, matrix);

    return vec;
}


triangle RotateX(triangle &tris, double ang){
    mat4x4 Rotated_MatrixX;
    triangle NewTris;

    Rotated_MatrixX.m[0][0] = 1;
    Rotated_MatrixX.m[1][1] = cos(degToRad(ang));
    Rotated_MatrixX.m[1][2] = sin(degToRad(ang));
    Rotated_MatrixX.m[2][1] = -sin(degToRad(ang));
    Rotated_MatrixX.m[2][2] = cos(degToRad(ang));
    Rotated_MatrixX.m[3][3] = 1;

    MultiplyMatrix(tris.p[0], NewTris.p[0], Rotated_MatrixX);
    MultiplyMatrix(tris.p[1], NewTris.p[1], Rotated_MatrixX);
    MultiplyMatrix(tris.p[2], NewTris.p[2], Rotated_MatrixX);

    return NewTris;
}

triangle RotateY(triangle &tris, double ang){
    mat4x4 Rotated_MatrixY;
    triangle NewTris;

    Rotated_MatrixY.m[0][0] = cos(degToRad(ang));
    Rotated_MatrixY.m[1][1] = 1;
    Rotated_MatrixY.m[0][2] = -sin(degToRad(ang));
    Rotated_MatrixY.m[2][0] = sin(degToRad(ang));
    Rotated_MatrixY.m[2][2] = cos(degToRad(ang));

    MultiplyMatrix(tris.p[0], NewTris.p[0], Rotated_MatrixY);
    MultiplyMatrix(tris.p[1], NewTris.p[1], Rotated_MatrixY);
    MultiplyMatrix(tris.p[2], NewTris.p[2], Rotated_MatrixY);

    return NewTris;
}

triangle RotateZ(triangle &tris, double ang){
    mat4x4 Rotated_MatrixZ;
    triangle NewTris;

    Rotated_MatrixZ.m[0][0] = cos(degToRad(ang));
    Rotated_MatrixZ.m[0][1] = sin(degToRad(ang));
    Rotated_MatrixZ.m[1][0] = -sin(degToRad(ang));
    Rotated_MatrixZ.m[1][1] = cos(degToRad(ang));
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



void FillTriangle(	int x1, int y1, int x2, int y2, int x3, int y3, SDL_Renderer *renderer)
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
            SDL_RenderDrawLine(renderer, ax, i, bx, i);
		}	
	}		
}