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
	SDL_Surface* texture = NULL;
	color color = {255, 255, 255};

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

			if (!bHasTexture)
			{
				if (line[0] == 'f')
				{
					int f[3];
					s >> junk >> f[0] >> f[1] >> f[2];
					tris.push_back((triangle){ verts[f[0] - 1], verts[f[1] - 1], verts[f[2] - 1] });
				}
			}
			else
			{
				if (line[0] == 'f')
				{
					s >> junk;

					std::string tokens[6];
					int nTokenCount = -1;


					while (!s.eof())
					{
						char c = s.get();
						if (c == ' ' || c == '/')
							nTokenCount++;
						else
							tokens[nTokenCount].append(1, c);
					}

					tokens[nTokenCount].pop_back();


					tris.push_back((triangle){ verts[stoi(tokens[0]) - 1], verts[stoi(tokens[2]) - 1], verts[stoi(tokens[4]) - 1],
						texs[stoi(tokens[1]) - 1], texs[stoi(tokens[3]) - 1], texs[stoi(tokens[5]) - 1] });

				}

			}
		}
		return true;
	}
};


int LZ(int x){
	int z = x;
	if(z < 0){
		z = 0;
	}
	return z;
}


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


float dist(vec3 &p, vec3 &plane_n, vec3 &plane_p){
	vec3 n = Vector_Normalise(p);
	return (plane_n.x * p.x + plane_n.y * p.y + plane_n.z * p.z - Vector_DotProduct(plane_n, plane_p));
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
			out_tri1.c =  in_tri.c;

			out_tri1.p[0] = *inside_points[0];
			out_tri1.t[0] = *inside_tex[0];

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
			out_tri1.c =  in_tri.c;

			out_tri2.c =  in_tri.c;

			out_tri1.p[0] = *inside_points[0];
			out_tri1.p[1] = *inside_points[1];
			out_tri1.t[0] = *inside_tex[0];
			out_tri1.t[1] = *inside_tex[1];

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


mat4x4 Matrix_Projection(int SCREEN_WIDTH, int SCREEN_HEIGHT, float FOV, float Zfar, float Znear)
{
	mat4x4 matrix;
    matrix.m[0][0] = (SCREEN_HEIGHT/SCREEN_HEIGHT) * (1/tan(degToRad(FOV/2)));
    matrix.m[1][1] = (1/tan(degToRad(FOV/2)));
    matrix.m[2][2] = Zfar / (Zfar-Znear);
    matrix.m[3][2] = (-Zfar * Znear) / (Zfar - Znear);
    matrix.m[2][3] = 1;
    matrix.m[3][3] = 0;
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


Uint32 getpixel(SDL_Surface *surface, int x, int y)
{
    int bpp = surface->format->BytesPerPixel;
    Uint8 *p = (Uint8 *)surface->pixels + y * surface->pitch + x * bpp;

	switch (bpp)
	{
		case 1:
			return *p;
			break;

		case 2:
			return *(Uint16 *)p;
			break;

		case 3:
			if (SDL_BYTEORDER == SDL_BIG_ENDIAN)
				return p[0] << 16 | p[1] << 8 | p[2];
			else
				return p[0] | p[1] << 8 | p[2] << 16;
				break;

			case 4:
				return *(Uint32 *)p;
				break;

			default:
				return 0;
		}
}


void FillTriangle(	int x1, int y1, float u1, float v1, int x2, int y2, float u2, float v2, int x3, int y3, float u3, float v3, float z1, float z2, float z3, SDL_Renderer *renderer, SDL_Surface* texture, color c, float *pDepthBuffer, int SCREEN_WIDTH, int SCREEN_HEIGHT)
{

	if (y2 < y1)
	{
		std::swap(y1, y2);
		std::swap(x1, x2);
		std::swap(u1, u2);
		std::swap(v1, v2);
		std::swap(z1, z2);
	}

	if (y3 < y1)
	{
		std::swap(y1, y3);
		std::swap(x1, x3);
		std::swap(u1, u3);
		std::swap(v1, v3);
		std::swap(z1, z3);
	}

	if (y3 < y2)
	{
		std::swap(y2, y3);
		std::swap(x2, x3);
		std::swap(u2, u3);
		std::swap(v2, v3);
		std::swap(z2, z3);
	}

	int dy1 = y2 - y1;
	int dx1 = x2 - x1;
	float dv1 = v2 - v1;
	float du1 = u2 - u1;
	float dw1 = z2 - z1;

	int dy2 = y3 - y1;
	int dx2 = x3 - x1;
	float dv2 = v3 - v1;
	float du2 = u3 - u1;
	float dw2 = z3 - z1;

	float tex_u, tex_v, tex_w;

	float dax_step = 0, dbx_step = 0,
		du1_step = 0, dv1_step = 0,
		du2_step = 0, dv2_step = 0,
		dw1_step=0, dw2_step=0;

	if (dy1) dax_step = dx1 / (float)abs(dy1);
	if (dy2) dbx_step = dx2 / (float)abs(dy2);

	if (dy1) du1_step = du1 / (float)abs(dy1);
	if (dy1) dv1_step = dv1 / (float)abs(dy1);
	if (dy1) dw1_step = dw1 / (float)abs(dy1);

	if (dy2) du2_step = du2 / (float)abs(dy2);
	if (dy2) dv2_step = dv2 / (float)abs(dy2);
	if (dy2) dw2_step = dw2 / (float)abs(dy2);

	if (dy1)
	{
		for (int i = y1; i <= y2; i++)
		{
			int ax = x1 + (float)(i - y1) * dax_step;
			int bx = x1 + (float)(i - y1) * dbx_step;

			float tex_su = u1 + (float)(i - y1) * du1_step;
			float tex_sv = v1 + (float)(i - y1) * dv1_step;
			float tex_sw = z1 + (float)(i - y1) * dw1_step;

			float tex_eu = u1 + (float)(i - y1) * du2_step;
			float tex_ev = v1 + (float)(i - y1) * dv2_step;
			float tex_ew = z1 + (float)(i - y1) * dw2_step;

			if (ax > bx)
			{
				std::swap(ax, bx);
				std::swap(tex_su, tex_eu);
				std::swap(tex_sv, tex_ev);
				std::swap(tex_sw, tex_ew);
			}

			tex_u = tex_su;
			tex_v = tex_sv;
			tex_w = tex_sw;

			float tstep = 1.0f / ((float)(bx - ax));
			float t = 0.0f;

			for (int j = ax; j < bx; j++)
			{
				tex_u = (1.0f - t) * tex_su + t * tex_eu;
				tex_v = (1.0f - t) * tex_sv + t * tex_ev;
				tex_w = (1.0f - t) * tex_sw + t * tex_ew;
				if (tex_w > pDepthBuffer[i*SCREEN_WIDTH + j])
				{
					if(texture != NULL){
						Uint32 data = getpixel(texture, (int)((float)(tex_u/tex_w)*(float)(texture->w)), (int)((float)(tex_v/tex_w)*(float)(texture->h)));
						SDL_Color col;
						SDL_GetRGB(data, texture->format, &col.r, &col.g, &col.b);
						SDL_SetRenderDrawColor(renderer, LZ(c.r+col.r-255), LZ(c.g+col.g-255), LZ(c.b+col.b-255), 255);
					}
					SDL_RenderDrawPoint(renderer, j, i);
					pDepthBuffer[i*SCREEN_WIDTH + j] = tex_w;
				}
				t += tstep;
			}
		}
	}

	dy1 = y3 - y2;
	dx1 = x3 - x2;
	dv1 = v3 - v2;
	du1 = u3 - u2;
	dw1 = z3 - z2;

	if (dy1) dax_step = dx1 / (float)abs(dy1);
	if (dy2) dbx_step = dx2 / (float)abs(dy2);

	du1_step = 0, dv1_step = 0;
	if (dy1) du1_step = du1 / (float)abs(dy1);
	if (dy1) dv1_step = dv1 / (float)abs(dy1);
	if (dy1) dw1_step = dw1 / (float)abs(dy1);

	if (dy1)
	{
		for (int i = y2; i <= y3; i++)
		{
			int ax = x2 + (float)(i - y2) * dax_step;
			int bx = x1 + (float)(i - y1) * dbx_step;

			float tex_su = u2 + (float)(i - y2) * du1_step;
			float tex_sv = v2 + (float)(i - y2) * dv1_step;
			float tex_sw = z2 + (float)(i - y2) * dw1_step;

			float tex_eu = u1 + (float)(i - y1) * du2_step;
			float tex_ev = v1 + (float)(i - y1) * dv2_step;
			float tex_ew = z1 + (float)(i - y1) * dw2_step;

			if (ax > bx)
			{
				std::swap(ax, bx);
				std::swap(tex_su, tex_eu);
				std::swap(tex_sv, tex_ev);
				std::swap(tex_sw, tex_ew);
			}

			tex_u = tex_su;
			tex_v = tex_sv;
			tex_w = tex_sw;

			float tstep = 1.0f / ((float)(bx - ax));
			float t = 0.0f;

			for (int j = ax; j < bx; j++)
			{
				tex_u = (1.0f - t) * tex_su + t * tex_eu;
				tex_v = (1.0f - t) * tex_sv + t * tex_ev;
				tex_w = (1.0f - t) * tex_sw + t * tex_ew;
				if (tex_w > pDepthBuffer[i*SCREEN_WIDTH + j])
				{
					if(texture != NULL){
						Uint32 data = getpixel(texture, (int)((float)(tex_u/tex_w)*(float)(texture->w)), (int)((float)(tex_v/tex_w)*(float)(texture->h)));
						SDL_Color col;
						SDL_GetRGB(data, texture->format, &col.r, &col.g, &col.b);
						SDL_SetRenderDrawColor(renderer, LZ(c.r+col.r-255), LZ(c.g+col.g-255), LZ(c.b+col.b-255), 255);
					}
					SDL_RenderDrawPoint(renderer, j, i);
					pDepthBuffer[i*SCREEN_WIDTH + j] = tex_w;
				}
				t += tstep;
			}
		}	
	}		
}


void Frame(SDL_Renderer *renderer, std::vector<mesh> mesh_collection, mat4x4 Projection_Matrix, player_t Camera, float *pDepthBuffer, vec3 light, int SCREEN_WIDTH, int SCREEN_HEIGHT){

    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    SDL_RenderClear(renderer);

    std::fill_n(pDepthBuffer, SCREEN_HEIGHT*SCREEN_WIDTH, 0.0f);

    float l_length = sqrt(light.x*light.x + light.y*light.y + light.z*light.z);
    light.x /= l_length; light.y /= l_length; light.z /= l_length;

    vec3 vLookDir;
    vLookDir.z = 1;
    mat4x4 LookRx, LookRy;

    LookRx = Matrix_RotateX(-Camera.angy);
    LookRy = Matrix_RotateY(Camera.angx);
    LookRy = Matrix_MultiplyMatrix(LookRx, LookRy);
    vLookDir = Vector_MultiplyMatrix(vLookDir, LookRy);

    vec3 pc1;
    vec3 pc2;
    vec3 pc3;
    vec3 pc4;
    vec3 pc5;
    vec3 pc6;
    vec3 pc7;
    pc2.y = 1;
    pc3.y = (float)SCREEN_HEIGHT - 1;
    pc4.y = -1;
    pc5.x = 1;
    pc6.x = (float)SCREEN_WIDTH - 1;
    pc7.x = -1;

    

    vec3 vUp;
    vUp.y = 1;
    vec3 vTarget; 
    vTarget.x = Camera.x + vLookDir.x;
    vTarget.y = Camera.y + vLookDir.y;
    vTarget.z = Camera.z + vLookDir.z;

    mat4x4 matCamera = Matrix_CameraPos(Camera, vTarget, vUp);
    mat4x4 matView = Matrix_QuickInverse(matCamera);


    for(int obj = 0; obj < mesh_collection.size(); obj++){
        mesh Mesh = mesh_collection[obj];

        mat4x4 RotX, RotY, RotZ, Size, Move, Mesh_Matrix;

        RotX = Matrix_RotateX(Mesh.rotation.x+180);
        RotY = Matrix_RotateY(Mesh.rotation.y);
        RotZ = Matrix_RotateZ(Mesh.rotation.z);
        Size = Scale(Mesh.size);
        Move = Translate(Mesh.position);

        Mesh_Matrix = Matrix_MultiplyMatrix(RotZ, RotY);
        Mesh_Matrix = Matrix_MultiplyMatrix(Mesh_Matrix, RotX);
        Mesh_Matrix = Matrix_MultiplyMatrix(Mesh_Matrix, Size);
        Mesh_Matrix = Matrix_MultiplyMatrix(Mesh_Matrix, Move);

        for(int i = 0; i < Mesh.tris.size(); i++){
            triangle transformation, projection, view;
            vec3 normal, line1, line2;

            transformation = Mesh.tris[i];

            for(int j = 0; j < 3; j++){
                transformation.p[j] = Vector_MultiplyMatrix(transformation.p[j], Mesh_Matrix);
            }

            line1 = Vector_Sub(transformation.p[1], transformation.p[0]);
            line2 = Vector_Sub(transformation.p[2], transformation.p[0]);

            normal.x = line1.y * line2.z - line1.z * line2.y;
            normal.y = line1.z * line2.x - line1.x * line2.z;
            normal.z = line1.x * line2.y - line1.y * line2.x;

            float l = sqrt(normal.x*normal.x + normal.y*normal.y + normal.z*normal.z);
            normal.x /= l; normal.y /= l; normal.z /= l;

            if(normal.x * (transformation.p[0].x - Camera.x) + normal.y * (transformation.p[0].y - Camera.y) + normal.z * (transformation.p[0].z - Camera.z) < 0)
            {
                for(int j = 0; j < 3; j++){
                    view.p[j] = Vector_MultiplyMatrix(transformation.p[j], matView);
                }

                float light_dp = normal.x*light.x+normal.y*light.y+normal.z*light.z;

                int shade = (int)(light_dp*250 + 5);
                if(shade < 5 || shade > 255){
                    shade = 5;
                }


                int nClippedTriangles = 0;
                triangle clipped[2];
                vec3 p_1;
                vec3 p_2;
                p_2.z = 1;
                p_1.z = 0.1;
                nClippedTriangles = Triangle_ClipAgainstPlane(p_1, p_2, view, clipped[0], clipped[1]);

                for (int n = 0; n < nClippedTriangles; n++)
                {
                    projection.p[0] = Vector_MultiplyMatrix(clipped[n].p[0], Projection_Matrix);
                    projection.p[1] = Vector_MultiplyMatrix(clipped[n].p[1], Projection_Matrix);
                    projection.p[2] = Vector_MultiplyMatrix(clipped[n].p[2], Projection_Matrix);

                    projection.t[0] = transformation.t[0];
                    projection.t[1] = transformation.t[1];
                    projection.t[2] = transformation.t[2];

                    projection.t[0].x = projection.t[0].x / projection.p[0].w;
                    projection.t[1].x = projection.t[1].x / projection.p[1].w;
                    projection.t[2].x = projection.t[2].x / projection.p[2].w;

                    projection.t[0].y = projection.t[0].y / projection.p[0].w;
                    projection.t[1].y = projection.t[1].y / projection.p[1].w;
                    projection.t[2].y = projection.t[2].y / projection.p[2].w;

                    projection.t[0].z = 1/projection.p[0].w;
                    projection.t[1].z = 1/projection.p[1].w;
                    projection.t[2].z = 1/projection.p[2].w;

					projection.p[0] = Vector_Div(projection.p[0], projection.p[0].w);
					projection.p[1] = Vector_Div(projection.p[1], projection.p[1].w);
					projection.p[2] = Vector_Div(projection.p[2], projection.p[2].w);

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
                            case 0:	nTrisToAdd = Triangle_ClipAgainstPlane(pc1, pc2, test, clipped[0], clipped[1]); break;
                            case 1:	nTrisToAdd = Triangle_ClipAgainstPlane(pc3, pc4, test, clipped[0], clipped[1]); break;
                            case 2:	nTrisToAdd = Triangle_ClipAgainstPlane(pc1, pc5, test, clipped[0], clipped[1]); break;
                            case 3:	nTrisToAdd = Triangle_ClipAgainstPlane(pc6, pc7, test, clipped[0], clipped[1]); break;
                            }
                            for (int w = 0; w < nTrisToAdd; w++)
                                listTriangles.push_back(clipped[w]);
                        }
                        nNewTriangles = listTriangles.size();
                    }


                    for (std::list<triangle>::iterator t = listTriangles.begin(); t != listTriangles.end(); t++)
                    {
						if(!Mesh.texture){
							SDL_SetRenderDrawColor(renderer, LZ(t->c.r + Mesh.color.r - 255), LZ(t->c.g + Mesh.color.g - 255), LZ(t->c.b + Mesh.color.b - 255), 255);
						}
                        FillTriangle(t->p[0].x, t->p[0].y, t->t[0].x, t->t[0].y, t->p[1].x, t->p[1].y, t->t[1].x, t->t[1].y, t->p[2].x, t->p[2].y, t->t[2].x, t->t[2].y, t->t[0].z, t->t[1].z, t->t[2].z, renderer, Mesh.texture, (color){LZ(t->c.r + Mesh.color.r - 255), LZ(t->c.g + Mesh.color.g - 255), LZ(t->c.b + Mesh.color.b - 255)}, pDepthBuffer, SCREEN_WIDTH, SCREEN_HEIGHT);
                    }
                }
            }
        }
    }

    SDL_RenderPresent(renderer);
}