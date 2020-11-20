

#include "Render.h"

#include <Windows.h>
#include <GL\GL.h>
#include <GL\GLU.h>
#include <corecrt_math.h>
#include <corecrt_math_defines.h>
#include < math.h >
#include <vector>
double direction = 1;
#include "Function.h"



using std::vector;


double t_max = 0;

void Render(double delta_time)
{
	t_max += delta_time / 5*direction; //t_max становится = 1 за 5 секунд
	if (t_max > 1) direction = -1;
	if (t_max < 0) direction = 1;

	//double P1[] = { 0,0,0 };
	//double P4[] = { 5,6,3 };
	//double R1[] = { -1,-1,-1 };
	//double R4[] = { 1,0,1 };
	//HermiteCurve(P1, P4, R1, R4);
	//PlaneCurse(P1,P4,R1,R4,t_max);


	double cP1[] = { 0,0,0 };
	double cP4[] = { 5,4,9 };
	double cR1[] = { 2,3,-2};
	double cR4[] = { 1,1,0 };

	HermiteCurve(cP1, cP4, cR1, cR4);
	PlaneCurse(cP1, cP4, cR1, cR4, t_max);
	//double P0[] = { 0,0,0 };
	//double P1с[] = { 6,2,3 };
	//double P2[] = { -1,-1,-1 };
	//double P3[] = { 1,0,1 };

	//BezierCurve(P0,P1с,P2,P3);

	//double cP0[] = { 0,0,0 };
	//double cP1с[] = { 10,3,3 };
	//double cP2[] = { 5,10,-1 };
	//double cP3[] = { 1,3,1 };

	//BezierCurve(cP0, cP1с, cP2, cP3);
	

	double points[4][4][3] = {
	{{0, 9, 0}, {3, 9, -3}, {6, 9, -3}, {9, 9, 0}},
	{{0, 6, -3}, {3, 6, -3}, {6, 6, -3}, {9, 6, -3}},
	{{0, 3, -3}, {3, 3, -3}, {6, 3, 5}, {9, 3, -3}},
	{{0, 0, -2}, {3, 0, -3}, {6, 0, -3}, {9, 0, -2}} };
	BezierSurface(points);
	
}   

