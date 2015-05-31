#ifndef _UTILITY_H_
#define _UTILITY_H_
#include "random.h"

#define TRUE 1
#define FALSE 0
#define C_PI                                3.141592653589793238462643383279502884197169399375105820974944592308
#define			max(a,b)					((a) > (b) ? (a) : (b))

typedef enum enum_SphereSampleMethod{UniformHemiSphere4D, UniformHemiSphere2D, UniformHemiSphere3D, UniformHemiSphereRand, IncorrectUniformHemiSphere4D, CubeReject4D, CubeIncorrectReject4D, PixieHemiSphere4D, PixieHemiSphereRand, ConsineSample4D}SphereSampleMethod;
#define SPHERE_SAMPLE_METHOD_STR "UniformHemiSphere4D,UniformHemiSphere2D,UniformHemiSphere3D,UniformHemiSphereRand,IncorrectUniformHemiSphere4D,CubeReject4D,CubeIncorrectReject4D,PixieHemiSphere4D, PixieHemiSphereRand,ConsineSample4D"

typedef enum enum_TwoDGridSampleMethod{OneSobol1D, TwoSobol1D, Sobol2D, Sobol3D, Sobol4D, Rand, Concentric}TwoDGridSampleMethod;
#define TWO_D_GRID_SAMPLE_METHOD_STR "OneSobol1D, TwoSobol1D, Sobol2D, Sobol3D, Sobol4D, Rand, Concentric"

typedef enum enum_DisplayMethod{Sphere, TwoDGrid}DisplayMethod;
#define DISPLAY_METHOD_STR "Sphere, TwoDGrid"

inline	float	dotvv(const float *s1,const float *s2) {
	return (float) (s1[0]*s2[0] + s1[1]*s2[1] + s1[2]*s2[2]);
}

inline	void	normalizev(float *v) {
	const double	l	=	1 / sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);

	v[0]	=	(float) (v[0]*l);
	v[1]	=	(float) (v[1]*l);
	v[2]	=	(float) (v[2]*l);
}
inline	void	movvv(float *dest,const float *src) {
	dest[0]		=	src[0];
	dest[1]		=	src[1];
	dest[2]		=	src[2];
}
inline	void	crossvv(float *result,const float *s1,const float *s2) {
	result[0]	=	(s1[1]*s2[2] - s1[2]*s2[1]);
	result[1]	=	(s1[2]*s2[0] - s1[0]*s2[2]);
	result[2]	=	(s1[0]*s2[1] - s1[1]*s2[0]);
}
inline	void	mulvf(float *result,const float m) {
	result[0]	*=	m;
	result[1]	*=	m;
	result[2]	*=	m;
}
inline	void	addvv(float *result,const float *s1) {
	result[0]	+=	s1[0];
	result[1]	+=	s1[1];
	result[2]	+=	s1[2];
}
inline	void	mulvf(float *result,const float *v,const float m) {
	result[0]	=	v[0]*m;
	result[1]	=	v[1]*m;
	result[2]	=	v[2]*m;
}


// uniform sample hemisphere on normal N side using random sphere
template<int dimension>
void UniformHemiSphere_func(float *out, CSobol<dimension> &generator)
{
	float	Dtmp[dimension];
	while(TRUE)
	{
		generator.get(Dtmp);
		out[2] = 1.f - 2.f * Dtmp[0];
		float r = sqrtf(max(0.f, 1.f-out[2]*out[2]));
		float phi = 2.f * C_PI * Dtmp[1];
		out[0] = r * cosf(phi);
		out[1] = r * sinf(phi);

		if (out[2] > 0)	break;// on positive Z side
	}
	assert(dotvv(out, out) >0.9 && dotvv(out, out) < 1.1);
}
void UniformHemiSphereRand_func(float *out);
void IncorrectUniformHemiSphere_func(float *out, CSobol<4> &generator);
void CubeIncorrectReject_func(float *out, CSobol<4> &generator);
void CubeReject_func(float *out, CSobol<4> &generator);
void PixieHemiSphere4D_func(float *R,const float *Z,CSobol<4> &generator);
void PixieHemiSphereRand_func(float *R,const float *Z);
void ConsineSampleHemiSphere_func(float *out, CSobol<4> &generator);

void ConvertQuaternionToMatrix(const float *quat, float *mat);
void ConcentricSampleDisk(float *out, CSobol<4> &generator);

#endif
