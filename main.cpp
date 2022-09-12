#include "BasicGraphics.h"
#include <iostream>
#include "Triangle.h"
#include "Square.h"
#include "RegularPolygon.h"
#include "Sphere.h"
#include "Helix.h"
#include "Vector3D.h"
#include "BasicPhysics.h"
#include "SpringConstraint.h"
#include <cstdlib> 

using namespace std;

#ifndef FU_H
#define FU_H
#include <functional>
#endif
const int X = 150;
const int Y = 300;
#define M_PI acos(-1.0)

const double width = 0.1;
const double height = 0.1;
const double sml = 1E-2;
double deltaX = width / X;
double deltaTime = 1.0/50;
double c = deltaX / deltaTime;
double nu = 0.5E-5;
double tau = 3.0 * nu * deltaTime / (deltaX * deltaX) + 0.5;
int layout[X][Y];
double rho[X][Y];
double u[X][Y][2];
double f[2][X][Y][9];
double feq[X][Y][9];
int old = 0;
int current = 1;
double rfreq[20][2];
int RCOS = 20;
double heatmaprange = 0.01;
double hmcentre = 3.235;
const double suction = 0.6;
const double inlet = 2.0;
const bool displayvel =false;
const double clamp_epsilon = 0.05;
const double clamp_a = clamp_epsilon / exp(1);
const double clamp_b = 1 / clamp_epsilon;
const int nTestBodies = 500;
double testBodies [nTestBodies][3];
const double testSpeed = 200.0;
const double spawnHeight = 200;
const double testLifetime = 10;
const double testLifeTimeRange = 5;

double w_i(int i) {
	if (i == 0) return 4.0/9;
	if (i <= 4) return 1.0 / 9;
	return 1.0 / 36;
}


double ei_dot_u(int i, double u_x, double u_y) {
	double e_x = (i == 0 || i == 2 || i == 4) ? 0.0 :
		(i == 1 || i == 5 || i == 8) ? 1.0 : -1.0;
	double e_y = (i == 0 || i == 1 || i == 3) ? 0.0 :
		(i == 2 || i == 5 || i == 6) ? 1.0 : -1.0;
	return e_x * u_x + e_y * u_y;
}

double si_u(int i, double u_x, double u_y) {
	double eiDotU = ei_dot_u(i, u_x, u_y);
	return w_i(i) * (6 * c * eiDotU + 9 * eiDotU * eiDotU - 3 * (u_x * u_x + u_y * u_y)) / (2 * c * c);
}
double smoothclamp(double x) {
	if (x > clamp_epsilon) return x;
	return clamp_a * exp(clamp_b * x);
}

float* heatMap(double val, double centre, double range,double alpha) {
	return new float[] {(float) pow(2, -pow(min(0.0,val - (centre + range)),2)*alpha / (range* range)), 

		(float) pow(2, -pow(val - centre, 2) *2.0* alpha / (range * range)), 

		(float)pow(2, -pow(max(0.0,val - (centre - range)), 2) * alpha / (range * range))};
}
void updateVelocities(int x, int y, int i, int t_x, int t_y, int r_i) {
	int t_tile = layout[t_x][t_y];
	if (t_tile == 0) f[current][t_x][t_y][i] = f[old][x][y][i];
	else if (t_tile == 1) f[current][x][y][r_i] = f[old][x][y][i];
	else if (t_tile == 2) f[current][x][y][r_i] = suction;
	else if (t_tile == 3) f[current][x][y][r_i] = inlet;

}

void inittestBody(int index) {
		testBodies[index][0] = 1.0 * X * rand() / RAND_MAX;
		testBodies[index][1] = 1.0 * spawnHeight* rand() / RAND_MAX;
		testBodies[index][2] = testLifetime + testLifeTimeRange*rand()/RAND_MAX;
}


void updateTestBodies(float deltaTime) {
	for (int i = 0; i < nTestBodies; i++) {
		double x = testBodies[i][0];
		double y = testBodies[i][1];
		/*
		double xr = x - floor(x);
		double yr = y - floor(y);
		int yc = (int)x;
		int xc = (int)y;
		double x1 = yr + xr;
		double x2 = yr - xr + 1.0;
		double u_x = (x1 * u[xc + 1][yc + 1][0] + (2.0 - x1) * u[xc][yc][0] + x2 * u[xc][yc + 1][0] + (2.0 - x2) * u[xc + 1][yc][0]) / 4.0;
		double u_y = (x1 * u[xc + 1][yc + 1][1] + (2.0 - x1) * u[xc][yc][1] + x2 * u[xc][yc + 1][1] + (2.0 - x2) * u[xc + 1][yc][1]) / 4.0;*/

		double newx = x + u[(int)round(x)][(int)round(y)][0] * testSpeed;
		double newy = y + u[(int)round(x)][(int)round(y)][1]*testSpeed;
		testBodies[i][2] -= deltaTime;
		while (testBodies[i][2] <0.0 || newx < 0.0 || newy < 0.0 || newx > X-1 || newy > Y-1 || layout[(int)round(newx)][(int)round(newy)] != 0) {
			inittestBody(i);
			newx = testBodies[i][0];
			newy = testBodies[i][1];
		}
		
		testBodies[i][0] = newx;
		testBodies[i][1] = newy;
		
	}
}
void update(float deltaTime) {
	double dt = deltaTime;

	for (int x = 0; x < X; x++) 
		for (int y = 0; y < Y; y++) {
			if (layout[x][y] != 0) continue;

			updateVelocities(x, y, 1, x + 1, y, 3);
			updateVelocities(x, y, 2, x, y + 1, 4);
			updateVelocities(x, y, 3, x - 1, y, 1);
			updateVelocities(x, y, 4, x, y - 1, 2);
			updateVelocities(x, y, 5, x + 1, y + 1, 7);
			updateVelocities(x, y, 6, x - 1, y + 1, 8);
			updateVelocities(x, y, 7, x - 1, y - 1, 5);
			updateVelocities(x, y, 8, x + 1, y - 1, 6);

		}


	double avrho = 0.0;
	double rhovar = 0.0;
	int N = 0;
	for (int x = 0; x < X; x++)
		for (int y = 0; y < Y; y++)
		{
			if (layout[x][y] != 0) continue;
			rho[x][y] = 0.0;
			for (int i = 0; i < 9; i++)
				rho[x][y] += f[current][x][y][i];
			N++;
			avrho += rho[x][y];
			rhovar += (rho[x][y] - avrho/N) * (rho[x][y] - avrho/N);
			u[x][y][0] = c*(f[current][x][y][1] + f[current][x][y][5] + f[current][x][y][8]
				- f[current][x][y][3] - f[current][x][y][7] - f[current][x][y][6])/rho[x][y];

			u[x][y][1] = c*(f[current][x][y][2] + f[current][x][y][5] + f[current][x][y][6]
				- f[current][x][y][4] - f[current][x][y][7] - f[current][x][y][8])/rho[x][y];

			for (int i = 0; i < 9; i++) {
				double feq= w_i(i) * rho[x][y] + rho[x][y] * si_u(i, u[x][y][0], u[x][y][1]);
				f[current][x][y][i] = smoothclamp((1 - 1.0 / tau) * f[current][x][y][i] + feq / tau);
				if (f[current][x][y][i] < 0) cout << "ono";
			}


		}

	updateTestBodies(deltaTime);
	for (int x = 0; x < X; x++)
		for (int y = 0; y < Y; y++) {
			if (layout[x][y] != 0)continue;
			
			delete BasicGraphics::grid[x][y];
			if (displayvel)
				BasicGraphics::grid[x][y] = (float*)heatMap(atan2(u[x][y][0], u[x][y][1]), 0.0, M_PI,1.8);
			else BasicGraphics::grid[x][y] =(float*)  heatMap(rho[x][y],hmcentre,heatmaprange,0.3);
		}
	for (int i = 0; i < nTestBodies; i++) {
		delete BasicGraphics::grid[(int)round(testBodies[i][0])][(int)round(testBodies[i][1])];
		BasicGraphics::grid[(int)round(testBodies[i][0])][(int)round(testBodies[i][1])] = new float[] {0.0, 0.0, 0.0};
	}
	old = (old + 1) % 2;
	current = (old + 1) % 2;
	cout << avrho / N << endl;
	hmcentre = avrho / N;
	heatmaprange = 0.5*sqrt(rhovar / N);
}

double initdist(int x, int y,int i) {
	double xd = max(0.0,1.0 - 4.0 * (1.0*x - X / 2) * (1.0*x - X / 2) / (1.0*X * X));
	double yd = max(0.0,1.0 - 4.0 * (1.0 * y - 1.0 * Y / 2) * (1.0 * y - 1.0 * Y / 2) / (1.0 * Y * Y));
	return 0.01*xd*yd + 2.0*rand()/RAND_MAX;
	//double rval = 0.0;
	//for (int i = 0; i < RCOS; i++) {
	//	rval += rfreq[i][1] * cos(rfreq[i][0] * x / X);
	//}
	//return max(0.0,rval + 1.0 -1.0* x * x / (1.0* X * X));
	//return (x < X / 3) ? 1.0*rand()/RAND_MAX : 0.0;
	double r = 1.0 * 2 * sqrt((x - X / 2) * (x - X / 2) + (y - Y / 2) * (y - Y / 2)) / X;
	double theta = x == X/2 ? 0.0 : atan2(y - Y / 2, x - X / 2);
	double u_x = sin(theta)*r;
	double u_y = cos(theta) * r;
	return 0.01*max(0.001,ei_dot_u(i, u_x, u_y));
}


int main(int argc, char** argv) {
	srand (time(NULL));
	cout << (int) 2.9 << endl;

	BasicGraphics::initialise(argc, argv);

	BasicGraphics::setGridSize(X, Y);
	for (int i = 0; i < nTestBodies; i++)
		inittestBody(i);
	for (int i = 0; i < RCOS; i++) {
		rfreq[i][0] = 1.0* rand() / RAND_MAX;
		rfreq[i][1] = 0.01 * rand() / RAND_MAX;
	}

	for (int x = 0; x < X; x++)
		for (int y = 0; y < Y; y++)
			if ((x == 0 || x == X - 1 || y == 0 || y == Y-1) && (x < 2*X/10 || x > 8*X/10) || (x-X/2)*(x-X/2)+ (y-Y/2)*(y-Y/2) < 500) {
				layout[x][y] = 1;
				BasicGraphics::grid[x][y] = new float[] {1.0, 1.0, 1.0};
			}
			else if (y == Y - 1 && x == min(max(2*X/10,x),8*X/10)) {
				layout[x][y] = 2;
			}else if (y == 0 && x == min(max(2*X/10,x),8*X/10)) {
				layout[x][y] = 3;
			}
			else {
				layout[x][y] = 0;
				for (int i = 0; i < 9; i++)	
					f[old][x][y][i] = 2.0 * initdist(x, y,i);
				}

	BasicGraphics::startVisuals(new function<void(float)>(update));
	return 0;
}
