#include <stdio.h>
#include <stdlib.h>
#include <crtdbg.h>
#include <math.h>
#include <assert.h>
#include "utility.h"

// incorrect version
void IncorrectUniformHemiSphere_func(float *out, CSobol<4> &generator)
{
	float	Dtmp[4];
	while(TRUE)
	{
		generator.get(Dtmp);
		out[2] = cosf(C_PI * Dtmp[0]);
		float r = sqrtf(max(0.f, 1.f-out[2]*out[2]));
		float phi = 2.f * C_PI * Dtmp[1];
		out[0] = r * cosf(phi);
		out[1] = r * sinf(phi);

		if (out[2] > 0)	break;// on positive Z side
	}
	assert(dotvv(out, out) >0.9 && dotvv(out, out) < 1.1);
}
void CubeReject_func(float *out, CSobol<4> &generator)
{
	float	Dtmp[4];
	while(TRUE) {
		generator.get(Dtmp);

		// Sample a uniformly distributed point on a sphere
		Dtmp[0]	=	2*Dtmp[0]-1;
		Dtmp[1]	=	2*Dtmp[1]-1;
		Dtmp[2]	=	2*Dtmp[2]-1;

		if (dotvv(Dtmp, Dtmp) <= 1 && Dtmp[2] > 0)	break;
	}
	normalizev(Dtmp);
	movvv(out, Dtmp);
}
void CubeIncorrectReject_func(float *out, CSobol<4> &generator)
{
	float	Dtmp[4];
	while(TRUE) {
		generator.get(Dtmp);

		// Sample a uniformly distributed point on a sphere
		Dtmp[0]	=	2*Dtmp[0]-1;
		Dtmp[1]	=	2*Dtmp[1]-1;
		Dtmp[2]	=	2*Dtmp[2]-1;

		if (Dtmp[2] > 0)	break;
	}
	normalizev(Dtmp);
	movvv(out, Dtmp);
}
void PixieHemiSphere4D_func(float *R,const float *Z,CSobol<4> &generator) {
	float	P[4];
	float 	Po[3];
	float	cosa;
	float	sina;

	while(TRUE) {
		generator.get(P);

		// Sample a uniformly distributed point on a sphere
		P[0]	=	2*P[0]-1;
		P[1]	=	2*P[1]-1;
		P[2]	=	2*P[2]-1;

		// did we get something inside the unit sphere and non-zero
		const float l = dotvv(P,P);
		if (l < 1 && l > 1e-6f)	break;
	}

	cosa			=	1 - P[3]*(1 - (float) cos(C_PI/2.f));
	sina			=	sqrtf(1 - cosa*cosa);

	// Po is orthagonal to N
	crossvv(Po,P,Z);
	// Po is unit length
	normalizev(Po);
	// Construct the sample vector
	mulvf(R,Z,cosa);
	mulvf(Po,sina);
	addvv(R,Po);
}
void PixieHemiSphereRand_func(float *R,const float *Z) {
	float	P[4];
	float 	Po[3];
	float	cosa;
	float	sina;

	while(TRUE) {
		P[0] = rand() / (float)RAND_MAX;
		P[1] = rand() / (float)RAND_MAX;
		P[2] = rand() / (float)RAND_MAX;
		P[3] = rand() / (float)RAND_MAX;

		// Sample a uniformly distributed point on a sphere
		P[0]	=	2*P[0]-1;
		P[1]	=	2*P[1]-1;
		P[2]	=	2*P[2]-1;

		// did we get something inside the unit sphere and non-zero
		const float l = dotvv(P,P);
		if (l < 1 && l > 1e-6f)	break;
	}

	cosa			=	1 - P[3]*(1 - (float) cos(C_PI/2.f));
	sina			=	sqrtf(1 - cosa*cosa);

	// Po is orthagonal to N
	crossvv(Po,P,Z);
	// Po is unit length
	normalizev(Po);
	// Construct the sample vector
	mulvf(R,Z,cosa);
	mulvf(Po,sina);
	addvv(R,Po);
}
// Routine to convert a quaternion to a 4x4 matrix
// ( input: quat = float[4]  output: mat = float[4*4] )
void ConvertQuaternionToMatrix(const float *quat, float *mat)
{
	float yy2 = 2.0f * quat[1] * quat[1];
	float xy2 = 2.0f * quat[0] * quat[1];
	float xz2 = 2.0f * quat[0] * quat[2];
	float yz2 = 2.0f * quat[1] * quat[2];
	float zz2 = 2.0f * quat[2] * quat[2];
	float wz2 = 2.0f * quat[3] * quat[2];
	float wy2 = 2.0f * quat[3] * quat[1];
	float wx2 = 2.0f * quat[3] * quat[0];
	float xx2 = 2.0f * quat[0] * quat[0];
	mat[0*4+0] = - yy2 - zz2 + 1.0f;
	mat[0*4+1] = xy2 + wz2;
	mat[0*4+2] = xz2 - wy2;
	mat[0*4+3] = 0;
	mat[1*4+0] = xy2 - wz2;
	mat[1*4+1] = - xx2 - zz2 + 1.0f;
	mat[1*4+2] = yz2 + wx2;
	mat[1*4+3] = 0;
	mat[2*4+0] = xz2 + wy2;
	mat[2*4+1] = yz2 - wx2;
	mat[2*4+2] = - xx2 - yy2 + 1.0f;
	mat[2*4+3] = 0;
	mat[3*4+0] = mat[3*4+1] = mat[3*4+2] = 0;
	mat[3*4+3] = 1;
}

void ConcentricSampleDisk(float *out, CSobol<4> &generator)
{    
	float r, theta, tmp[4];

	generator.get(tmp);
    // Map uniform random numbers to $[-1,1]^2$
    float sx = 2 * tmp[0] - 1;
    float sy = 2 * tmp[1] - 1;

    // Map square to $(r,\theta)$

    // Handle degeneracy at the origin
    if (sx == 0.0 && sy == 0.0) {
        out[0] = 0.0;
        out[1] = 0.0;
        return;
    }
    if (sx >= -sy) {
        if (sx > sy) {
            // Handle first region of disk
            r = sx;
            if (sy > 0.0) theta = sy/r;
            else          theta = 8.0f + sy/r;
        }
        else {
            // Handle second region of disk
            r = sy;
            theta = 2.0f - sx/r;
        }
    }
    else {
        if (sx <= sy) {
            // Handle third region of disk
            r = -sx;
            theta = 4.0f - sy/r;
        }
        else {
            // Handle fourth region of disk
            r = -sy;
            theta = 6.0f + sx/r;
        }
    }
    theta *= C_PI / 4.f;
    out[0] = r * cosf(theta);
    out[1] = r * sinf(theta);
	out[2] = 0.f;
}
void ConsineSampleHemiSphere_func(float *out, CSobol<4> &generator)
{
	ConcentricSampleDisk(out, generator);
	out[2] = sqrtf(max(0.f, 1.f - out[0]*out[0] - out[1]*out[1]));
}
void UniformHemiSphereRand_func(float *out)
{
	float	Dtmp[2];
	while(TRUE)
	{
		Dtmp[0] = rand() / (float)RAND_MAX;
		Dtmp[1] = rand() / (float)RAND_MAX;
		out[2] = 1.f - 2.f * Dtmp[0];
		float r = sqrtf(max(0.f, 1.f-out[2]*out[2]));
		float phi = 2.f * C_PI * Dtmp[1];
		out[0] = r * cosf(phi);
		out[1] = r * sinf(phi);

		if (out[2] > 0)	break;// on positive Z side
	}
	assert(dotvv(out, out) >0.9 && dotvv(out, out) < 1.1);
}