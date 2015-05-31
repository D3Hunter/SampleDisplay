#include <Windows.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "AntTweakBar.h"
#include "glut.h"

#include "display.h"
#include "utility.h"

// #include <vld.h>

////////////////////////////////////////////////////////
// Variables used in AntTweakBar
////////////////////////////////////////////////////////
float g_zoom = 0.6, g_lineWidth = 0.05;
float g_solidColor[] = {0.3, 0.3, 0.3, 1.0};
float g_wiredColor[] = {0.6, 0.6, 0.6, 1.0};
float g_samplePointColor[] = {0.8, 0, 0, 1.0};
float g_rotation[] = {0, 0, 0, 1};
int g_level = MAX_LEVEL, // display level
	g_displayMethod = (int)Sphere, // display method
	g_sphereSampleMethod = (int)UniformHemiSphere4D, // sphere sample method
	g_2dGridSampleMethod = (int)OneSobol1D, // 2d grid sample method
	g_pointSize = 4, // size of sample point
	g_getCount = 0; // invoke Sobol's get method count times after get one sample(for debugging one bug)

float allSphereSampleData[MAX_POINT * 3], all2DGridSampleData[MAX_POINT * 3];
float *sphereDisplayPtr[MAX_LEVEL], *twoDGridDisplayPtr[MAX_LEVEL];

void GenerateSphereSampleData(SphereSampleMethod m)
{
	float N[3] = {0, 0, 1}, tmp[4];
	CSobol<2> random2d;
	CSobol<3> random3d;
	CSobol<4> random4d;
	float minRadius = 2.f;
	for(int i = 0; i < MAX_POINT; i++)
	{
		switch(m)
		{
		case UniformHemiSphere4D:
			UniformHemiSphere_func<4>(allSphereSampleData + 3 * i, random4d);
			for(int i = 0; i < g_getCount; i++) random4d.get(tmp);
			break;
		case UniformHemiSphere2D:
			UniformHemiSphere_func<2>(allSphereSampleData + 3 * i, random2d);
			for(int i = 0; i < g_getCount; i++) random2d.get(tmp);
			break;
		case UniformHemiSphere3D:
			UniformHemiSphere_func<3>(allSphereSampleData + 3 * i, random3d);
			for(int i = 0; i < g_getCount; i++) random3d.get(tmp);
			break;
		case UniformHemiSphereRand:
			UniformHemiSphereRand_func(allSphereSampleData + 3 * i);
			break;
		case IncorrectUniformHemiSphere4D:
			IncorrectUniformHemiSphere_func(allSphereSampleData + 3 * i, random4d);
			break;
		case CubeReject4D:
			CubeReject_func(allSphereSampleData + 3 * i, random4d);
			break;
		case CubeIncorrectReject4D:
			CubeIncorrectReject_func(allSphereSampleData + 3 * i, random4d);
			break;
		case PixieHemiSphere4D:
			PixieHemiSphere4D_func(allSphereSampleData + 3 * i, N, random4d);
			break;
		case PixieHemiSphereRand:
			PixieHemiSphereRand_func(allSphereSampleData + 3 * i, N);
			break;
		case ConsineSample4D:
			ConsineSampleHemiSphere_func(allSphereSampleData + 3 * i, random4d);
			break;
		}
		if(sqrtf(dotvv(allSphereSampleData + 3 * i, allSphereSampleData + 3 * i)) < minRadius) minRadius = sqrtf(dotvv(allSphereSampleData + 3 * i, allSphereSampleData + 3 * i));
	}
	printf("MIN RADIUS IS %f\n", minRadius);
}

void Generate2DGridSampleData(TwoDGridSampleMethod m)
{
	float tmp[4];
	CSobol<1> random1d_1, random1d_2;
	CSobol<2> random2d;
	CSobol<3> random3d;
	CSobol<4> random4d;

	for(int i = 0; i < MAX_POINT; i++)
	{
		switch(m)
		{
		case OneSobol1D:
			random1d_1.get(tmp + 0);
			random1d_1.get(tmp + 1);
			break;
		case TwoSobol1D:
			random1d_1.get(tmp + 0);
			random1d_2.get(tmp + 1);
			break;
		case Sobol2D:
			random2d.get(tmp);
			break;
		case Sobol3D:
			random3d.get(tmp);
			break;
		case Sobol4D:
			random4d.get(tmp);
			break;
		case Rand:
			tmp[0] = rand() / (float)RAND_MAX;
			tmp[1] = rand() / (float)RAND_MAX;
			break;
		case Concentric:
			ConcentricSampleDisk(&all2DGridSampleData[3 * i], random4d);
			all2DGridSampleData[3 * i + 0] *= 0.5f;
			all2DGridSampleData[3 * i + 1] *= 0.5f;
			continue;
		}
		all2DGridSampleData[3 * i + 0] = tmp[0]-0.5f;
		all2DGridSampleData[3 * i + 1] = tmp[1]-0.5f;
		all2DGridSampleData[3 * i + 2] = 0.f;
	}
}
void InitDisplayData()
{
	int base = 0;
	for(int i = 0; i < MAX_LEVEL; i++)
	{
		sphereDisplayPtr[i] = allSphereSampleData + base;
		twoDGridDisplayPtr[i] = all2DGridSampleData + base;
		base += 3*(1 << i);
	}

	// default Uniform sample
	GenerateSphereSampleData((SphereSampleMethod)g_sphereSampleMethod);
	Generate2DGridSampleData((TwoDGridSampleMethod)g_2dGridSampleMethod);
}

void init(void)
{
	// below three is necessary for light0 to work
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);

	glDepthFunc(GL_LESS);
	glShadeModel(GL_SMOOTH);
	glClearColor(0, 0, 0, .5);
	glShadeModel(GL_SMOOTH);

	InitDisplayData();
}

void displaySphereSampleData()
{
	float colorAmb[4]    = { 0.1, 0.1, 0.1, 1.0 };
	glMaterialfv(GL_FRONT, GL_DIFFUSE, g_solidColor);
	glMaterialfv(GL_FRONT, GL_AMBIENT, colorAmb);
	glutSolidSphere(BACK_SPHERE_RADIUS - 0.002, 32, 32);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, g_wiredColor);
	glLineWidth(g_lineWidth);
	glutWireSphere(BACK_SPHERE_RADIUS, 32, 32);

	glMaterialfv(GL_FRONT, GL_DIFFUSE, g_samplePointColor);
	//glMaterialfv(GL_FRONT, GL_AMBIENT, colorNone_a);
	glPointSize(g_pointSize);
	glBegin(GL_POINTS);
	for(int i = 0; i < (1 << (g_level-1)); i++)
	{
		glNormal3fv(sphereDisplayPtr[(g_level-1)] + 3 * i);
		glVertex3fv(sphereDisplayPtr[(g_level-1)] + 3 * i);
	}
	glEnd();
}
void display2DGridSampleData()
{
	float colorAmb[4]    = { 0.1, 0.1, 0.1, 1.0 };
	glMaterialfv(GL_FRONT, GL_DIFFUSE, g_solidColor);
	glMaterialfv(GL_FRONT, GL_AMBIENT, colorAmb);
	
	glBegin(GL_POLYGON);
		glVertex3f(0.5f, 0.5f, -0.02f);
		glVertex3f(-0.5f, 0.5f, -0.02f);
		glVertex3f(-0.5f, -0.5f, -0.02f);
		glVertex3f(0.5f, -0.5f, -0.02f);
	glEnd();

	glMaterialfv(GL_FRONT, GL_DIFFUSE, g_samplePointColor);
	//glMaterialfv(GL_FRONT, GL_AMBIENT, colorNone_a);
	glPointSize(g_pointSize);
	glNormal3f(0, 0, 1);
	glBegin(GL_POINTS);
	for(int i = 0; i < (1 << (g_level-1)); i++)
	{
		glVertex3fv(twoDGridDisplayPtr[(g_level-1)] + 3 * i);
	}
	glEnd();
}
void display(void)
{
	long t1,t2;
	float mat[4*4]; // rotation matrix
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	ConvertQuaternionToMatrix(g_rotation, mat);

	glPushMatrix();
	glMultMatrixf(mat);
	glScalef(g_zoom, g_zoom, g_zoom);

	switch(g_displayMethod)
	{
	case Sphere:
		displaySphereSampleData();
		break;
	case TwoDGrid:
		display2DGridSampleData();
		break;
	}

	glPopMatrix();

	// Draw tweak bars
	TwDraw();

	// Present frame buffer
	glutSwapBuffers();

	// Recall Display at next frame
	glutPostRedisplay();
}

void reShape(int w, int h)
{
	float ratio = (0 == h) ? w : (1.0 * w) / h;
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	glViewport(0, 0, w, h);
	gluPerspective(15, ratio, MIN_PLAN_OF_EYE, 1000);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(0, 0, 5,  0, 0, -1,  0, 1, 0);

	GLfloat lightpos[] = { 0, 0, 30, 1 };
	GLfloat amb[] = { 1.f, 1.f, 1.f, 1 };
	glLightfv(GL_LIGHT0, GL_POSITION, lightpos);
	glLightfv(GL_LIGHT0, GL_AMBIENT, amb);

	// Send the new window size to AntTweakBar
	TwWindowSize(w, h);
}
// Function called at exit
void Terminate(void)
{ 
	//	glDeleteLists(SHAPE_TEAPOT, NUM_SHAPES);
	TwTerminate();
}
static void TW_CALL cb_SetGetCount(const void *value, void *clientData)
{
	*(float *)clientData = *(float *)value;

	// reinit sample data
	GenerateSphereSampleData((SphereSampleMethod)g_sphereSampleMethod);
}
static void TW_CALL cb_GetGetCount(void *value, void *clientData)
{
	*(float *)value = *(float *)clientData;
}
static void TW_CALL cb_SetSphereSampleMethod(const void *value, void *clientData)
{
	*(float *)clientData = *(float *)value;

	// reinit sample data
	GenerateSphereSampleData((SphereSampleMethod)g_sphereSampleMethod);
}
static void TW_CALL cb_GetSphereSampleMethod(void *value, void *clientData)
{
	*(float *)value = *(float *)clientData;
}
static void TW_CALL cb_Set2DGridSampleMethod(const void *value, void *clientData)
{
	*(float *)clientData = *(float *)value;

	// reinit sample data
	Generate2DGridSampleData((TwoDGridSampleMethod)g_2dGridSampleMethod);
}
static void TW_CALL cb_Get2DGridSampleMethod(void *value, void *clientData)
{
	*(float *)value = *(float *)clientData;
}
int main(int argc, char **argv)
{	
	TwBar *bar; // Pointer to the tweak bar

	setlocale(LC_ALL,"Chinese-simplified");//设置中文环境

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
	glutInitWindowPosition(0, 0);
	glutCreateWindow(argv[0]);
	init();
	glutDisplayFunc(display);
	glutReshapeFunc(reShape);

	atexit(Terminate);  // Called after glutMainLoop ends
	// Initialize AntTweakBar
	TwInit(TW_OPENGL, NULL);
	// Set GLUT event callbacks
	// - Directly redirect GLUT mouse button events to AntTweakBar
	glutMouseFunc((GLUTmousebuttonfun)TwEventMouseButtonGLUT);
	// - Directly redirect GLUT mouse motion events to AntTweakBar
	glutMotionFunc((GLUTmousemotionfun)TwEventMouseMotionGLUT);
	// - Directly redirect GLUT mouse "passive" motion events to AntTweakBar (same as MouseMotion)
	glutPassiveMotionFunc((GLUTmousemotionfun)TwEventMouseMotionGLUT);
	// - Directly redirect GLUT key events to AntTweakBar
	glutKeyboardFunc((GLUTkeyboardfun)TwEventKeyboardGLUT);
	// - Directly redirect GLUT special key events to AntTweakBar
	glutSpecialFunc((GLUTspecialfun)TwEventSpecialGLUT);
	// - Send 'glutGetModifers' function pointer to AntTweakBar;
	//   required because the GLUT key event functions do not report key modifiers states.
	TwGLUTModifiersFunc(glutGetModifiers);
	// Create a tweak bar
	bar = TwNewBar("ControlBar");
	TwDefine(" GLOBAL help='Loop And Modified Butterfly Subdivision Methods Demo' "); // Message added to the help bar.
	TwDefine(" ControlBar size='200 400' color='96 216 224' "); // change default tweak bar size and color

	// rotation, zoom, sample level, PointSize
	TwAddVarRW(bar, "ObjRotation", TW_TYPE_QUAT4F, &g_rotation, 
		" label='Object rotation' opened=true help='Change the object orientation.' ");
	TwAddVarRW(bar, "Zoom", TW_TYPE_FLOAT, &g_zoom, "min=0.01 max=10 step=0.01");
	TwAddVarRW(bar, "Level", TW_TYPE_INT32, &g_level, "min=1 max="MAX_LEVEL_S" step=1");
	TwAddSeparator(bar, "Sep1", "");
	TwAddVarRW(bar, "PointSize", TW_TYPE_INT32, &g_pointSize, "min=1 max=20 step=1");
	TwAddVarCB(bar, "GetAfterRandom", TW_TYPE_INT32, cb_SetGetCount, cb_GetGetCount, &g_getCount, "min=0 max=20 step=1");
	TwAddVarRW(bar, "LineWidth", TW_TYPE_FLOAT, &g_lineWidth, "min=0.01 max=1.5 step=0.02");
	TwAddVarRW(bar, "SolidColor", TW_TYPE_COLOR4F, &g_solidColor,  "");
	TwAddVarRW(bar, "WiredColor", TW_TYPE_COLOR4F, &g_wiredColor,  "");
	TwAddVarRW(bar, "SampleColor", TW_TYPE_COLOR4F, &g_samplePointColor, "");

	TwType clusterType1;
	clusterType1 = TwDefineEnumFromString("DisplayMethod", DISPLAY_METHOD_STR);
	TwAddVarRW(bar, "DisplayMethod", clusterType1, &g_displayMethod, "");
	clusterType1 = TwDefineEnumFromString("SphereSampleMethod", SPHERE_SAMPLE_METHOD_STR);
	TwAddVarCB(bar, "SphereSampleMethod", clusterType1, cb_SetSphereSampleMethod, cb_GetSphereSampleMethod, &g_sphereSampleMethod, "");
	clusterType1 = TwDefineEnumFromString("2DGridSampleMethod", TWO_D_GRID_SAMPLE_METHOD_STR);
	TwAddVarCB(bar, "2DGridSampleMethod", clusterType1, cb_Set2DGridSampleMethod, cb_Get2DGridSampleMethod, &g_2dGridSampleMethod, "");


	glutMainLoop();
	float *a = (float *)malloc(100);
	float *b = new float;

	return 0;
}
