#pragma once

struct vec3
{
	float x = 0.0f;
	float y = 0.0f;
	float z = 0.0f;
	float w = 1.0f;
};

struct color
{
    int r,g,b;
};


struct player_t
{
    float x, y, z, angx, angy;
};

struct triangle
{
    vec3 p[3];
	vec3 t[3];
    color c;
};

struct mat4x4
{
    float m[4][4] = { 0 };
};


struct mesh
{
    std::vector<triangle> tris;

	vec3 position;
	vec3 size;
	vec3 rotation;

	bool LoadFile(std::string sFilename, bool bHasTexture = false)
	{
		std::ifstream f(sFilename);
		if (!f.is_open())
			return false;

		std::vector<vec3> verts;
		std::vector<vec3> texs;

		while (!f.eof())
		{
			char line[128];
			f.getline(line, 128);

			std::strstream s;
			s << line;

			char junk;

			if (line[0] == 'v')
			{
				if (line[1] == 't')
				{
					vec3 v;
					s >> junk >> junk >> v.x >> v.y;
					texs.push_back(v);
				}
				else
				{
					vec3 v;
					s >> junk >> v.x >> v.y >> v.z;
					verts.push_back(v);
				}
			}

			if (line[0] == 'f')
			{
				int f[3];
				s >> junk >> f[0] >> f[1] >> f[2];
				tris.push_back((triangle){ verts[f[0] - 1], verts[f[1] - 1], verts[f[2] - 1] });
			}
			
		}
		return true;
	}
};



float Distance(vec3 vec1, vec3 vec2){
    return sqrt((vec1.x-vec2.x) * (vec1.x-vec2.x) + (vec1.y-vec2.y)*(vec1.y-vec2.y));
}


float degToRad(float ang)
{
    return ang * (M_PI/180);;
}

float radToDeg(float rad)
{
    return rad * (180/M_PI);
}

float FixAng(float ang)
{
    if (ang > 359)
    {
        float mult = floor(ang / 360);
        ang = ang - (360 * mult);
    }
    if (ang < 0)
    {
        float ap = abs(ang);
        float mult = ceil(ap / 360);
        ang = (360 * mult) - ap;
    }
    return ang;
}

float FixY(float ang){
    if (ang > 89)
    {
        ang = 89;
    }
    if (ang < -89)
    {
        ang = -89;
    }
    return ang;
}



vec3 Vector_Add(vec3 &v1, vec3 &v2)
{
	vec3 nv;
	nv.x = v1.x + v2.x;
	nv.y = v1.y + v2.y;
	nv.z = v1.z + v2.z;
	return nv;
}

vec3 Vector_Sub(vec3 &v1, vec3 &v2)
{
	vec3 nv;
	nv.x = v1.x - v2.x;
	nv.y = v1.y - v2.y;
	nv.z = v1.z - v2.z;
	return nv;
}

vec3 Vector_Mul(vec3 &v1, float k)
{
	vec3 nv;
	nv.x = v1.x * k;
	nv.y = v1.y * k;
	nv.z = v1.z * k;
	return nv;
}

vec3 Vector_Div(vec3 &v1, float k)
{
	vec3 nv;
	nv.x = v1.x / k;
	nv.y = v1.y / k;
	nv.z = v1.z / k;
	return nv;
}

float Vector_DotProduct(vec3 &v1, vec3 &v2)
{
	return v1.x*v2.x + v1.y*v2.y + v1.z * v2.z;
}

float Vector_Length(vec3 &v)
{
	return Vector_DotProduct(v, v);
}

vec3 Vector_Normalise(vec3 &v)
{ 
	float l = Vector_DotProduct(v,v);
	vec3 nv;
	nv.x = v.x / l;
	nv.y = v.y / l;
	nv.z = v.z / l;
	return nv;
}

vec3 Vector_IntersectPlane(vec3 &plane_p, vec3 &plane_n, vec3 &lineStart, vec3 &lineEnd, float &t)
{
		plane_n = Vector_Normalise(plane_n);
		float plane_d = -Vector_DotProduct(plane_n, plane_p);
		float ad = Vector_DotProduct(lineStart, plane_n);
		float bd = Vector_DotProduct(lineEnd, plane_n);
		t = (-plane_d - ad) / (bd - ad);
		vec3 lineStartToEnd = Vector_Sub(lineEnd, lineStart);
		vec3 lineToIntersect = Vector_Mul(lineStartToEnd, t);
		return Vector_Add(lineStart, lineToIntersect);
}

float dist(vec3 &p, vec3 &plane_n, vec3 &plane_p){
	vec3 n = Vector_Normalise(p);
	return (plane_n.x * p.x + plane_n.y * p.y + plane_n.z * p.z - Vector_DotProduct(plane_n, plane_p));
}

	int Triangle_ClipAgainstPlane(vec3 plane_p, vec3 plane_n, triangle &in_tri, triangle &out_tri1, triangle &out_tri2)
	{
		plane_n = Vector_Normalise(plane_n);


		vec3* inside_points[3];  int nInsidePointCount = 0;
		vec3* outside_points[3]; int nOutsidePointCount = 0;
		vec3* inside_tex[3]; int nInsideTexCount = 0;
		vec3* outside_tex[3]; int nOutsideTexCount = 0;


		float d0 = dist(in_tri.p[0], plane_n, plane_p);
		float d1 = dist(in_tri.p[1], plane_n, plane_p);
		float d2 = dist(in_tri.p[2], plane_n, plane_p);

		if (d0 >= 0) { inside_points[nInsidePointCount++] = &in_tri.p[0]; inside_tex[nInsideTexCount++] = &in_tri.t[0]; }
		else {
			outside_points[nOutsidePointCount++] = &in_tri.p[0]; outside_tex[nOutsideTexCount++] = &in_tri.t[0];
		}
		if (d1 >= 0) {
			inside_points[nInsidePointCount++] = &in_tri.p[1]; inside_tex[nInsideTexCount++] = &in_tri.t[1];
		}
		else {
			outside_points[nOutsidePointCount++] = &in_tri.p[1];  outside_tex[nOutsideTexCount++] = &in_tri.t[1];
		}
		if (d2 >= 0) {
			inside_points[nInsidePointCount++] = &in_tri.p[2]; inside_tex[nInsideTexCount++] = &in_tri.t[2];
		}
		else {
			outside_points[nOutsidePointCount++] = &in_tri.p[2];  outside_tex[nOutsideTexCount++] = &in_tri.t[2];
		}

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
			out_tri1.p[0] = *inside_points[0];
			out_tri1.t[0] = *inside_tex[0];
			out_tri1.c =  in_tri.c;


			float t;
			out_tri1.p[1] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[0], t);
			out_tri1.t[1].x = t * (outside_tex[0]->x - inside_tex[0]->x) + inside_tex[0]->x;
			out_tri1.t[1].y = t * (outside_tex[0]->y - inside_tex[0]->y) + inside_tex[0]->y;
			out_tri1.t[1].z = t * (outside_tex[0]->z - inside_tex[0]->z) + inside_tex[0]->z;

			out_tri1.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[1], t);
			out_tri1.t[2].x = t * (outside_tex[1]->x - inside_tex[0]->x) + inside_tex[0]->x;
			out_tri1.t[2].y = t * (outside_tex[1]->y - inside_tex[0]->y) + inside_tex[0]->y;
			out_tri1.t[2].z = t * (outside_tex[1]->z - inside_tex[0]->z) + inside_tex[0]->z;

			return 1;
		}

		if (nInsidePointCount == 2 && nOutsidePointCount == 1)
		{
			out_tri1.p[0] = *inside_points[0];
			out_tri1.p[1] = *inside_points[1];
			out_tri1.t[0] = *inside_tex[0];
			out_tri1.t[1] = *inside_tex[1];
			out_tri1.c =  in_tri.c;
			out_tri2.c =  in_tri.c;

			float t;
			out_tri1.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[0], t);
			out_tri1.t[2].x = t * (outside_tex[0]->x - inside_tex[0]->x) + inside_tex[0]->x;
			out_tri1.t[2].y = t * (outside_tex[0]->y - inside_tex[0]->y) + inside_tex[0]->y;
			out_tri1.t[2].z = t * (outside_tex[0]->z - inside_tex[0]->z) + inside_tex[0]->z;

			out_tri2.p[0] = *inside_points[1];
			out_tri2.t[0] = *inside_tex[1];
			out_tri2.p[1] = out_tri1.p[2];
			out_tri2.t[1] = out_tri1.t[2];
			out_tri2.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[1], *outside_points[0], t);
			out_tri2.t[2].x = t * (outside_tex[0]->x - inside_tex[1]->x) + inside_tex[1]->x;
			out_tri2.t[2].y = t * (outside_tex[0]->y - inside_tex[1]->y) + inside_tex[1]->y;
			out_tri2.t[2].z = t * (outside_tex[0]->z - inside_tex[1]->z) + inside_tex[1]->z;
			return 2;
		}
		return 0;
	}


vec3 Vector_MultiplyMatrix(vec3 &i, mat4x4 &m)
{
	vec3 v;
	v.x = i.x * m.m[0][0] + i.y * m.m[1][0] + i.z * m.m[2][0] + i.w * m.m[3][0];
	v.y = i.x * m.m[0][1] + i.y * m.m[1][1] + i.z * m.m[2][1] + i.w * m.m[3][1];
	v.z = i.x * m.m[0][2] + i.y * m.m[1][2] + i.z * m.m[2][2] + i.w * m.m[3][2];
	v.w = i.x * m.m[0][3] + i.y * m.m[1][3] + i.z * m.m[2][3] + i.w * m.m[3][3];
	return v;
}

mat4x4 Matrix_MultiplyMatrix(mat4x4 &m1, mat4x4 &m2)
{
	mat4x4 matrix;
	for (int c = 0; c < 4; c++)
		for (int r = 0; r < 4; r++)
			matrix.m[r][c] = m1.m[r][0] * m2.m[0][c] + m1.m[r][1] * m2.m[1][c] + m1.m[r][2] * m2.m[2][c] + m1.m[r][3] * m2.m[3][c];
	return matrix;
}


mat4x4 Matrix_MakeIdentity()
{
	mat4x4 matrix;
	matrix.m[0][0] = 1.0f;
	matrix.m[1][1] = 1.0f;
	matrix.m[2][2] = 1.0f;
	matrix.m[3][3] = 1.0f;
	return matrix;
}


mat4x4 Translate(vec3 v2){
    mat4x4 matrix;

	matrix.m[0][0] = 1;
	matrix.m[1][1] = 1;
	matrix.m[2][2] = 1;
	matrix.m[3][3] = 1;
	matrix.m[3][0] = v2.x;
	matrix.m[3][1] = v2.y;
	matrix.m[3][2] = v2.z;

    return matrix;
}

mat4x4 Scale(vec3 v2){
    mat4x4 matrix;

	matrix.m[0][0] = v2.x;
	matrix.m[1][1] = v2.y;
	matrix.m[2][2] = v2.z;
	matrix.m[3][3] = 1;

    return matrix;
}



mat4x4 Matrix_RotateX(float ang)
{
	float fAngleRad = degToRad(ang);
	mat4x4 matrix;
	matrix.m[0][0] = 1;
	matrix.m[1][1] = cos(fAngleRad);
	matrix.m[1][2] = sin(fAngleRad);
	matrix.m[2][1] = -sin(fAngleRad);
	matrix.m[2][2] = cos(fAngleRad);
	matrix.m[3][3] = 1;
	return matrix;
}

mat4x4 Matrix_RotateY(float ang)
{
	mat4x4 matrix;
    matrix.m[0][0] = cos(degToRad(ang));
    matrix.m[1][1] = 1;
    matrix.m[0][2] = -sin(degToRad(ang));
    matrix.m[2][0] = sin(degToRad(ang));
    matrix.m[2][2] = cos(degToRad(ang));
	matrix.m[3][3] = 1;
	return matrix;
}

mat4x4 Matrix_RotateZ(float ang)
{
	float fAngleRad = degToRad(ang);
	mat4x4 matrix;
	matrix.m[0][0] = cos(fAngleRad);
	matrix.m[0][1] = sin(fAngleRad);
	matrix.m[1][0] = -sin(fAngleRad);
	matrix.m[1][1] = cos(fAngleRad);
	matrix.m[2][2] = 1;
	matrix.m[3][3] = 1;
	return matrix;
}


mat4x4 Matrix_CameraPos(player_t &pos, vec3 &target, vec3 &up){

    vec3 newForward; 
	newForward.x = target.x-pos.x;
	newForward.y = target.y-pos.y;
	newForward.z = target.z-pos.z;
    float NFl = sqrt(newForward.x*newForward.x + newForward.y*newForward.y + newForward.z*newForward.z);
    newForward.x /= NFl; newForward.y /= NFl; newForward.z /= NFl;

    float NFUdp = up.x*newForward.x+up.y*newForward.y+up.z*newForward.z;
    vec3 a; 
	a.x = newForward.x*NFUdp;
	a.y = newForward.y*NFUdp;
	a.z = newForward.z*NFUdp;
    vec3 newUp; 
	newUp.x = up.x-a.x;
	newUp.y = up.y-a.y;
	newUp.z = up.z-a.z;
    float NUl = sqrt(newUp.x*newUp.x + newUp.y*newUp.y + newUp.z*newUp.z);
    newUp.x /= NUl; newUp.y /= NUl; newUp.z /= NUl;

    vec3 newRight; 
	newRight.x = newUp.y * newForward.z - newUp.z * newForward.y;
	newRight.y = newUp.z * newForward.x - newUp.x * newForward.z;
	newRight.z = newUp.x * newForward.y - newUp.y * newForward.x;

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





void FillTriangle(	int x1, int y1, int x2, int y2, int x3, int y3, float z1, float z2, float z3, SDL_Renderer *renderer, float *pDepthBuffer, int SCREEN_WIDTH, int SCREEN_HEIGHT)
{

	if (y2 < y1)
	{
		std::swap(y1, y2);
		std::swap(x1, x2);
		std::swap(z1, z2);
	}

	if (y3 < y1)
	{
		std::swap(y1, y3);
		std::swap(x1, x3);
		std::swap(z1, z3);
	}

	if (y3 < y2)
	{
		std::swap(y2, y3);
		std::swap(x2, x3);
		std::swap(z2, z3);
	}

	int dy1 = y2 - y1;
	int dx1 = x2 - x1;
	float dw1 = z2 - z1;

	int dy2 = y3 - y1;
	int dx2 = x3 - x1;
	float dw2 = z3 - z1;

	float tex_u, tex_v, tex_w, tex_c;

	float dax_step = 0, dbx_step = 0,
		du1_step = 0, dv1_step = 0,
		du2_step = 0, dv2_step = 0,
		dw1_step=0, dw2_step=0,
		dc1_step=0, dc2_step=0;

	if (dy1) dax_step = dx1 / (float)abs(dy1);
	if (dy2) dbx_step = dx2 / (float)abs(dy2);

	if (dy1) dw1_step = dw1 / (float)abs(dy1);
	if (dy2) dw2_step = dw2 / (float)abs(dy2);


	if (dy1)
	{
		for (int i = y1; i <= y2; i++)
		{
			int ax = x1 + (float)(i - y1) * dax_step;
			int bx = x1 + (float)(i - y1) * dbx_step;

			float tex_sw = z1 + (float)(i - y1) * dw1_step;
			float tex_ew = z1 + (float)(i - y1) * dw2_step;


			if (ax > bx)
			{
				std::swap(ax, bx);
				std::swap(tex_sw, tex_ew);
			}

			tex_w = tex_sw;

			float tstep = 1.0f / ((float)(bx - ax));
			float t = 0.0f;

			for (int j = ax; j < bx; j++)
			{
				tex_w = (1.0f - t) * tex_sw + t * tex_ew;
				if (tex_w > pDepthBuffer[i*SCREEN_WIDTH + j])
				{
					SDL_RenderDrawPoint(renderer, j, i);
					pDepthBuffer[i*SCREEN_WIDTH + j] = tex_w;
				}
				t += tstep;
			}
		}
	}

	dy1 = y3 - y2;
	dx1 = x3 - x2;
	dw1 = z3 - z2;

	if (dy1) dax_step = dx1 / (float)abs(dy1);
	if (dy2) dbx_step = dx2 / (float)abs(dy2);

	if (dy1) dw1_step = dw1 / (float)abs(dy1);


	if (dy1)
	{
		for (int i = y2; i <= y3; i++)
		{
			int ax = x2 + (float)(i - y2) * dax_step;
			int bx = x1 + (float)(i - y1) * dbx_step;

			float tex_sw = z2 + (float)(i - y2) * dw1_step;
			float tex_ew = z1 + (float)(i - y1) * dw2_step;


			if (ax > bx)
			{
				std::swap(ax, bx);
				std::swap(tex_sw, tex_ew);
			}

			tex_w = tex_sw;

			float tstep = 1.0f / ((float)(bx - ax));
			float t = 0.0f;

			for (int j = ax; j < bx; j++)
			{
				tex_w = (1.0f - t) * tex_sw + t * tex_ew;
				if (tex_w > pDepthBuffer[i*SCREEN_WIDTH + j])
				{
					SDL_RenderDrawPoint(renderer, j, i);
					pDepthBuffer[i*SCREEN_WIDTH + j] = tex_w;
				}
				t += tstep;
			}
		}	
	}		
}