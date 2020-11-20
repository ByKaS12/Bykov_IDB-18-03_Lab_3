#pragma once
#include <vector>

#define v3f glVertex3f
#define v2f glVertex2f

using std::vector;
void PlaneCurse();
void taskONe(double A[], double B[], double C[]) {
	glBegin(GL_LINES);

	glVertex3dv(A);
	glVertex3dv(B);

	glVertex3dv(B);
	glVertex3dv(C);

	glVertex3dv(C);
	glVertex3dv(A);


	glEnd();

}
double* FindNormal(double T1[], double T2[], double T3[]) {
	double A[] = { T2[0] - T1[0],T2[1] - T1[1] ,T2[2] - T1[2] };
	double B[] = { T3[0] - T1[0],T3[1] - T1[1] ,T3[2] - T1[2] };
	double* array = new double[3];
	array[0] = (A[1] * B[2]) - (B[1] * A[2]);
	array[1] = (-(A[0] * B[2])) + (B[0] * A[2]);
	array[2] = (A[0] * B[1]) - (B[0] * A[1]);
	return array;

}
void createQuad(double A[], double B[], double C[], double D[], double color[]) {
	glBegin(GL_QUADS);
	glColor3dv(color);
	glVertex3dv(A);
	glVertex3dv(B);
	glVertex3dv(C);
	glVertex3dv(D);
	glEnd();

}
void taskFour(double A[], double B[], double C[]) {
	glBegin(GL_TRIANGLES);
	glVertex3dv(A);
	glColor3d(0.7, 0.2, 0.4);
	glVertex3dv(B);
	glColor3d(0.7, 0.2, 0);
	glVertex3dv(C);
	glColor3d(0.3, 0.1, 0.8);
	glEnd();

}
void taskFive(int i) {
	if (i == 0) {
		double A[] = { 1, 1, 0 };
		double B[] = { -1, 1, 0 };
		double C[] = { -1,-1, 0 };
		double D[] = { 1,-1, 0 };
		glBegin(GL_QUADS);
		glColor3d(0.2, 0.7, 0.7);
		glVertex3dv(A);
		glVertex3dv(B);
		glVertex3dv(C);
		glVertex3dv(D);
		glEnd();
	}
	else {
		double A[] = { 1, 1, 1 };
		double B[] = { -1, 1, 1 };
		double C[] = { -1,-1, 1 };
		double D[] = { 1,-1, 1 };
		glBegin(GL_QUADS);
		glColor3d(0.2, 0.7, 0.7);
		glVertex3dv(A);
		glVertex3dv(B);
		glVertex3dv(C);
		glVertex3dv(D);
		glEnd();
	}


}

void createTriangle(double A[], double B[], double C[], double color[]) {
	glBegin(GL_TRIANGLES);
	glColor3dv(color);
	glVertex3dv(A);
	glVertex3dv(B);
	glVertex3dv(C);

	glEnd();


}
void DrawCircle(double cx, double cy, double C[], double r, int num_segments, float z) {
	float x = 0.0f;
	float y = 0.0f;
	float prevX = x;
	float prevY = y - float(r);
	float floatR = float(r);
	glBegin(GL_TRIANGLES);
	glColor3f(sin(z), sin(z), cos(z));
	for (int ii = 0; ii <= num_segments; ii++) {
		float theta = 3.1415926f * float(ii) / float(num_segments);

		float newX = floatR * sin(theta);
		float newY = -floatR * cos(theta);
		glBegin(GL_TRIANGLES);
		glVertex3d(0, 0, z);
		glVertex3d(prevX, prevY, z);
		glVertex3d(newX, newY, z);
		glEnd();
		prevX = newX;
		prevY = newY;
	}

}

void point40(double A[], double B[], float z) {
	double C[] = { (B[0] + A[0]) / 2,(B[1] + A[1]) / 2 };


	double r = (B[0] - A[0]) / 2;
	int num_segments = 30;

	DrawCircle(A[0], A[1], C, r, num_segments, z);

}

void DrawCircle(float cx, float cy, float r, int num_segments, double z)
{
	glBegin(GL_LINE_LOOP);
	for (int ii = 0; ii < num_segments; ii++)
	{
		float theta = 3.1415926f * float(ii) / float(num_segments);//get the current angle

		float x = r * cosf(theta);//calculate the x component
		float y = r * sinf(theta);//calculate the y component
		if ((x + cx) != 6 && (y + cy) != 0)
			glVertex3d(x + cx, y + cy, z);//output vertex

	}
	glEnd();
}
void taskSix() {

	//
	double A[] = { -1, 1, 0 };
	double B[] = { -1, 1, 1 };
	double C[] = { 1, 1, 1 };
	double D[] = { 1,1, 0 };
	double color1[] = { 0.1,0.2,0,3 };
	//

	//
	double revA[] = { 1, -1, 0 };
	double revB[] = { 1, -1, 1 };
	double revC[] = { -1, -1, 1 };
	double revD[] = { -1,-1, 0 };
	double color2[] = { 0.2,0.4,0.7 };
	//
	double color3[] = { 0.1,0.1,0.2 };
	double color4[] = { 0.5,0.3,0.9 };
	taskFive(0);
	createQuad(A, B, C, D, color1);
	createQuad(revA, revB, revC, revD, color2);
	createQuad(revA, revB, C, D, color3);
	createQuad(A, B, revC, revD, color4);
	taskFive(1);
	double Zero[] = { 0,0,3 };

	taskFour(B, C, Zero);
	taskFour(revB, C, Zero);
	taskFour(B, revC, Zero);
	taskFour(revB, revC, Zero);

}

void taskFinnaly() {
	double A[] = { 1, -2, 3 };
	double B[] = { -4, 5, 8 };
	double C[] = { 1,2,3 };


	//taskONe(A, B, C);
	//taskFour(A, B, C);
	double angle = 7.2;
	for (int i = 0; i < 50; i++)
	{


		glPushMatrix();
		glRotated(angle, 0, 0, 1);
		glTranslated(0, 100, 0);
		glRotated(angle, 1, 0, 0);
		glScaled(5, 5, 5);
		taskSix();

		glPopMatrix();
		angle += 7.2;


	}
}
void labTwo() {
	double z = 0;
	double A1[] = { 1,3,z };//
	double B1[] = { 1,10,z };
	double C12[] = { 2,7,z };
	double color[] = { 0.2,0.3,0.5 };
	createTriangle(A1, B1, C12, color);
	double B2[] = { 6,5,z };
	double A2[] = { 5,12,z };
	createTriangle(C12, B2, A2, color);
	double D1[] = { 6,0,z };//
	createQuad(D1, A1, C12, B2, color);
	double C3[] = { 8,10,z };
	double D3[] = { 12,10,z };
	createQuad(D1, B2, C3, D3, color);
	z = 3;
	double A13[] = { 1,3,z };//
	double B13[] = { 1,10,z };
	double C123[] = { 2,7,z };
	createTriangle(A13, B13, C123, color);
	double B23[] = { 6,5,z };
	double A23[] = { 5,12,z };
	createTriangle(C123, B23, A23, color);
	double D13[] = { 6,0,z };//
	createQuad(D13, A13, C123, B23, color);
	double C33[] = { 8,10,z };
	double D33[] = { 12,10,z };
	createQuad(D13, B23, C33, D33, color);
	double color1[] = { 0.7,0.5,0.3 };
	createQuad(A1, B1, B13, A13, color1);
	createQuad(C12, B1, B13, C123, color1);
	createQuad(C12, C123, A23, A2, color1);
	createQuad(A2, A23, B23, B2, color1);
	createQuad(B2, B23, C33, C3, color1);
	//CREATE SHAR
	glPushMatrix();
	glTranslated(10, 10, 0);
	glRotated(90, 0, 0, 10);
	float i = 0.0f;
	point40(C3, D3, i);
	glPopMatrix();
	//glPushMatrix();

	////glTranslated(2.21429, -0.64286, 0);
	////glRotated(90, 0, 0, 10);
	//float ii = 0.0f;
	//while (ii <= 3.0f)
	//{
	//	
	//	DrawCircle(2.21429, -0.64286,  3.84, 50, ii);
	//	ii += 0.005f;
	//}
	//glPopMatrix();

	//END SHAR
	createQuad(D3, D33, D13, D1, color1);
	//createQuad(D1, D13, A13, A1, color1);


}

void point50(double z, vector <GLfloat> &X, vector <GLfloat> &Y) {
	glPushMatrix();
	GLfloat theta;
	GLfloat pi = acos(-1.0);
	GLfloat radius = 3.84f; // радиус
	GLfloat step = 1.0f; // чем больше шаг тем хуже диск
	int i = 0;
		glBegin(GL_LINE_STRIP);
		for (GLfloat a = 5.0f; a < 55.0f; a += step) {
			theta = 2.0f * pi * a / 180.0f;
			glColor4f(a / 360.0f, 1.0f, 1.0f - a / 360.0f, 1.0f);
			glVertex3f(radius * cos(theta), radius * sin(theta), z);
			if (i == X.size() || i == Y.size())
			{
				X.resize(X.size() + 5);
				Y.resize(Y.size() + 5);
			}
			X[i] = radius * cos(theta);
			Y[i] = radius * sin(theta);

			i++;

		}
		glEnd();

	glPopMatrix();

}
void Side3() {	
	vector <GLfloat> X(50);
	vector <GLfloat> Y(50);
	vector <GLfloat> X2(5);
	vector <GLfloat> Y2(5);
	double A[3];
	double B[3];
	double C[3];
	double D[3];
	double color[3] = {0,1,1};
	glPushMatrix();
	glTranslatef(2.21429f, -0.64286f, 0.0f);
	point50(0.0,X,Y);
	point50(3.0, X2, Y2);
	for (int i = 0; i < 49; i++)
	{
			A[0] = X[i];
			A[1] = Y[i];
			A[2] = 0;
			B[0] = X[i];
			B[1] = Y[i];
			B[2] = 3;
			C[0] = X[i+1];
			C[1] = Y[i+1];
			C[2] = 3;
			D[0] = X[i+1];
			D[1] = Y[i+1];
			D[2] = 0;
		createQuad(A,B,C,D,color);
	}
	glPopMatrix();

}
inline double HermiteCurvePoint(double P1, double P4, double R1, double R4, double t) {
	return  (P1 * ((2 * t * t * t) - (3 * t * t) + 1)) + (P4 * ((-2 * t * t * t) + 3 * t * t)) + (R1 * ((t * t * t) - (2 * t * t) + t)) + (R4 * (t * t * t) - (t * t));
}

void HermiteCurve(double P1[], double P4[], double R1[], double R4[]) {


	glLineWidth(3);
	glBegin(GL_LINE_STRIP);
	double Plast[3];
	for (double t = 0; t < 1.0001; t += 0.01)
	{
		double P[3];
		P[0] = HermiteCurvePoint(P1[0], P4[0], R1[0], R4[0], t);
		P[1] = HermiteCurvePoint(P1[1], P4[1], R1[1], R4[1], t);
		P[2] = HermiteCurvePoint(P1[2], P4[2], R1[2], R4[2], t);
		//glColor3d(sin(t),cos(t),tan(t));
		glVertex3dv(P);
		//PlaneCurse();
		if (t <= 1) {
			Plast[0] = P[0];
			Plast[1] = P[1];
			Plast[2] = P[2];
		}
			
	}
	glEnd();
	glLineWidth(1);
	glBegin(GL_LINE_STRIP);
	glVertex3dv(P1);
	glVertex3dv(R1);
	glEnd();
	glBegin(GL_LINE_STRIP);
	glVertex3dv(Plast);
	glVertex3dv(R4);
	glEnd();


}
inline double BezierCurvePoint(double P0, double P1, double P2, double P3, double t) {
	return (pow((1 - t), 3) * P0) + (3 * t * pow((1 - t), 2) * P1) + (3 * t * t * (1 - t) * P2)+(t * t * t * P3);
}
void BezierCurve(double P0[], double P1[], double P2[], double P3[]) {
	glLineWidth(3);
	glBegin(GL_LINE_STRIP);

	for (double t = 0; t < 1.001; t += 0.01)
	{
		double P[3];
		P[0] = BezierCurvePoint(P0[0], P1[0], P2[0], P3[0], t);
		P[1] = BezierCurvePoint(P0[1], P1[1], P2[1], P3[1], t);
		P[2] = BezierCurvePoint(P0[2], P1[2], P2[2], P3[2], t);
		//glColor3d(sin(t),cos(t),tan(t));
		glVertex3dv(P);
		
	}
	glEnd();
	glPointSize(5);
	glBegin(GL_POINTS);
	glVertex3dv(P0);
	glVertex3dv(P1);
	glVertex3dv(P2);
	glVertex3dv(P3);
	glEnd();
	glLineWidth(1);
	glBegin(GL_LINE_STRIP);
	glVertex3dv(P0);
	glVertex3dv(P1);
	glVertex3dv(P2);
	glVertex3dv(P3);
	glEnd();


}

double  PointsDistance(double P1[3], double P2[3]) {
	return sqrt(pow(P2[0] - P1[0], 2) + pow(P2[1] - P1[0], 2) + pow(P2[2] - P1[2], 2));
}
double* NormalizedVector(double V1[3]) {
	double nulcord[3] = { 0, 0, 0 };
	double length = PointsDistance(nulcord, V1);

	double relative_length = (1 / length);
	double x = V1[0] * relative_length;
	double y = V1[1] * relative_length;
	double z = V1[2] * relative_length;

	double* normalized = new double[3]{ x, y, z };
	return normalized;
}
double* CreateVector(double P1[3], double P2[3]) {
	double* V = new double[3];
	V[0] = P2[0] - P1[0];
	V[1] = P2[1] - P1[1];
	V[2] = P2[2] - P1[2];

	return V;
}
double ScalarProduct(double V1[3], double V2[3]) {
	double c = V1[0] * V2[0] + V1[1] * V2[1] + V1[2] * V2[2];
	double d = sqrt(pow(V1[0], 2) + pow(V1[1], 2) + pow(V1[2], 2)) * sqrt(pow(V2[0], 2) + pow(V2[1], 2) + pow(V2[2], 2));
	return c / d;
}

double* _crossProduct(double V1[3], double V2[3]) {
	double* V = new double[3];
	V[0] = V1[1] * V2[2] - V1[2] * V2[1];
	V[1] = V1[2] * V2[0] - V1[0] * V2[2];
	V[2] = V1[0] * V2[1] - V1[1] * V2[0];
	return V;
}
void PlaneCurse(double P1[], double P4[], double R1[], double R4[], double t_max) {
	boolean flag = true;

	double Vec[3] = { HermiteCurvePoint(P1[0], P4[0], R1[0], R4[0], t_max) ,HermiteCurvePoint(P1[1], P4[1], R1[1], R4[1], t_max),HermiteCurvePoint(P1[2], P4[2], R1[2], R4[2], t_max) };
	double	VecD[3] = { HermiteCurvePoint(P1[0], P4[0], R1[0], R4[0], t_max+0.01*direction) ,HermiteCurvePoint(P1[1], P4[1], R1[1], R4[1], t_max + 0.01* direction),HermiteCurvePoint(P1[2], P4[2], R1[2], R4[2], t_max + 0.01*direction) };
	double* dirV = CreateVector(Vec, VecD);
	double* dir = NormalizedVector(dirV);
	double orig[3] = { 1,0,0 };
	double rotXV[3] = { dir[0],dir[1],0 };
	double* rotX = NormalizedVector(rotXV);
	double cosU = ScalarProduct(orig, rotX);
	double* vecpr = _crossProduct(orig, rotX);
	double sinSign = vecpr[2] / abs(vecpr[2]);
	double U = acos(cosU) * 180 / acos(-1) * sinSign;
	double origZ[3] = { 0, 0, 1 };
	double cosZU = ScalarProduct(origZ, dir);
	double ZU = acos(dir[2]) * 180.0 / acos(-1) - 90;
	glPushMatrix();
	glTranslated(Vec[0], Vec[1], Vec[2]);
	glRotated(U, 0, 0, 1);
	glRotated(ZU, 0, 1, 0);
	glScalef(1.0 / 16.0, 1.0 / 16.0, 1.0 / 16.0);


// низ
	//glBegin(GL_POLYGON);
	//glColor3d(0, 0, 0);
	// glVertex3f(0.5, -0.5, -0.5);      
	// glVertex3f(0.5, 0.5, -0.5);      
 //    glVertex3f(-0.5, 0.5, -0.5);     
 //   glVertex3f(-0.5, -0.5, -0.5);      

	//glEnd();
	//glBegin(GL_LINES);
	//glColor3d(0, 0, 1);
	//glVertex3f(0, 0, 0);
	//glVertex3f(0, 0, -1);
	//glEnd();
	////верх
	//glBegin(GL_POLYGON);
	//glColor3d(0, 0, 0);
	//glVertex3f(0.5, -0.5, 0.5);
	//glVertex3f(0.5, 0.5, 0.5);
	//glVertex3f(-0.5, 0.5, 0.5);
	//glVertex3f(-0.5, -0.5, 0.5);
	//glEnd();
	//glBegin(GL_LINES);
	//glColor3d(0, 0, 1);
	//glVertex3f(0, 0, 0);
	//glVertex3f(0, 0, 1);
	//glEnd();

	//glBegin(GL_POLYGON);
	//glColor3d(0, 0, 0);
	//glVertex3f(0.5, -0.5, -0.5);
	//glVertex3f(0.5, 0.5, -0.5);
	//glVertex3f(0.5, 0.5, 0.5);
	//glVertex3f(0.5, -0.5, 0.5);
	//glEnd();
	//glBegin(GL_LINES);
	//glColor3d(1, 0, 0);
	//glVertex3f(0, 0, 0);
	//glVertex3f(1, 0, 0);
	//glEnd();
	//glBegin(GL_POLYGON);
	//glColor3d(0, 0, 0);
	//glVertex3f(-0.5, -0.5, 0.5);
	//glVertex3f(-0.5, 0.5, 0.5);
	//glVertex3f(-0.5, 0.5, -0.5);
	//glVertex3f(-0.5, -0.5, -0.5);
	//glEnd();
	//glBegin(GL_LINES);
	//glColor3d(1, 0, 0);
	//glVertex3f(0, 0, 0);
	//glVertex3f(-1, 0, 0);
	//glEnd();

	//glBegin(GL_POLYGON);
	//glColor3d(0, 0, 0);
	//glVertex3f(0.5, 0.5, 0.5);
	//glVertex3f(0.5, 0.5, -0.5);
	//glVertex3f(-0.5, 0.5, -0.5);
	//glVertex3f(-0.5, 0.5, 0.5);
	//glEnd();
	//glBegin(GL_LINES);
	//glColor3d(0, 1, 0);
	//glVertex3f(0, 0, 0);
	//glVertex3f(0, 1, 0);
	//glEnd();

	//glBegin(GL_POLYGON);
	//glColor3d(0, 0, 0);
	//glVertex3f(0.5, -0.5, -0.5);
	//glVertex3f(0.5, -0.5, 0.5);
	//glVertex3f(-0.5, -0.5, 0.5);
	//glVertex3f(-0.5, -0.5, -0.5);
	//glEnd();
	//glBegin(GL_LINES);
	//glColor3d(0, 1, 0);
	//glVertex3f(0, 0, 0);
	//glVertex3f(0, -1, 0);
	//glEnd();
	glRotated(-90,0,0,1);
	glBegin(GL_TRIANGLE_STRIP);
	v3f(-7.0, 0.0, 2.0);
	v3f(-1.0, 0.0, 3.0);
	v3f(-1.0, 7.0, 3.0);
	v3f(0.0, 0.0, 0.0);
	v3f(0.0, 8.0, 0.0);
	v3f(1.0, 0.0, 3.0);
	v3f(1.0, 7.0, 3.0);
	v3f(7.0, 0.0, 2.0);
	glEnd();
	glPopMatrix();
}

double FactorialNum(double num) {
	double res = 1;
	for (int i = 1; i <= num; i++) {
		res *= i;
	}
	return res;
}

double BernsteinPolynomial(double n, double i, double u) {
	return (FactorialNum(n) / (FactorialNum(i) * FactorialNum(n - i))) * pow(u, i) * pow(1 - u, n - i);
}

// для поверхности Безье третьего порядка
double* BezierSurfaceCalc(double points[4][4][3], double u, double v) {
	double* P = new double[3];
	//double P[3];
	P[0] = 0;
	P[1] = 0;
	P[2] = 0;

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			for (int c = 0; c < 3; c++) {
				P[c] += BernsteinPolynomial(3, i, u) * BernsteinPolynomial(3, j, v) * points[i][j][c];
			}
		}
	}
	return P;
}
void BezierSurface(double points[4][4][3]) {
	const double step = 0.1;
	const int size = 1 / 0.1 + 1; // если меняется step, то поменять здесь соотв. значение

	double* P;
	double PointCalc[size][size][3];
	int i = 0, j = 0;

	for (double u = 0; u <= 1; u += step) {
		for (double v = 0; v <= 1; v += step) {
			P = BezierSurfaceCalc(points, u, v);
			PointCalc[i][j][0] = P[0];
			PointCalc[i][j][1] = P[1];
			PointCalc[i][j][2] = P[2];

			delete[] P;
			j++;
		}
		j = 0;
		i++;
	}

	glPointSize(4);
	glColor3d(0, 0, 0);
	glBegin(GL_POINTS);
	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {
			glVertex3dv(PointCalc[i][j]);
		}
	}
	glEnd();

	glColor3d(1, 0, 1);
	for (i = 0; i < size; i++) {
		glBegin(GL_LINE_STRIP);
		for (j = 0; j < size; j++) {
			glVertex3dv(PointCalc[i][j]);
		}
		glEnd();
	}

	for (j = 0; j < size; j++) {
		glBegin(GL_LINE_STRIP);
		for (i = 0; i < size; i++) {
			glVertex3dv(PointCalc[i][j]);
		}
		glEnd();
	}
}