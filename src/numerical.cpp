//// deshi includes ////
#include "defines.h"
#include "deshi.h"
#include "core/memory.h"
#include "utils/string.h"
#include "utils/array.h"
#include "utils/map.h"
#include "math/Math.h"
#include "core/logging.h"



//
// Global Variables
//

u32 steps = 50; 

double L = 10;

double m   = 1;
double mm  = m*m;
double mmm = m*m*m;

double g   = 0.2;
double gg  = g*g;
double ggg = g*g*g;

double k = 20;
double delta = 1;
double t0 = 2;

//
// Math Helper Functions
//

typedef double (*MathFunc)(double p);
double riemannSum(double start, double end, double dx, MathFunc function) {
	double sum = 0;
	for (double i = start; i < end; i += dx)
		sum += function(i) * dx;
	return sum;
}

double sech(double p) {
	return 1 / cosh(p);
}

//
//Initial Condition Functions
//

double P(double t) {
	if (t < t0 - delta / 2) return -k / 2;
	if (t > t0 + delta / 2) return  k / 2;
	double inside = (t - t0) / (sqrt((t - t0 + delta / 2) * (t0 - t + delta / 2)));
	return k / 2 * tanh(inside);
}


array<double> PSample(double dt) {
	array<double> values;
	for (int i = 0; i < steps; i++) values.add(P(i * dt));
	return values;
}


double initialPos(double rho) {
	double M0 = mmm / (3 * M_SQRT_TWO * gg);
	double gamma = sqrt(1 + (k / 2) * (k / 2) / (M0 * M0));
	return m / (2 * g) * tanh(gamma * m * rho / M_SQRT_TWO);
}

array<double> initialPosSample() {
	double dx = L / (2 * steps + 1);
	array<double> values;
	for (int i = 0; i < 2 * steps + 1; i++) values.add(initialPos(-L / 2 + i * dx));
	return values;
}

double IIRiemannSum(double hk, double M0, double gamma) {
	double sum = 0;
	double dx = 1 / (100 * m);
	for (double i = -L * 100; i < L * 100; i += dx) {
		double s0 = sech(i);
		double s1 = sech(gamma * i);
		sum += s0 * s0 * s1 * s1 * dx;
	}
	return sum;
}

double initialMomentum(double rho) {
	double hk = k / 2;
	double M0 = mmm / (3 * M_SQRT_TWO * gg);
	double gamma = sqrt(1 + hk * hk / (M0 * M0));
	double II = IIRiemannSum(hk, M0, gamma);
	double sech0 = sech(sqrt(1 + hk * hk / (M0 * M0) * m * rho / M_SQRT_TWO));
	double sech1 = sech(m * rho / M_SQRT_TWO);
	return 3 * g * k / (4 * m) * ( sech0 * sech0 - (3 * II * sech1 * sech1) / 4);
}

array<double> initialMomentumSample() {
	double dx = L / (2 * steps + 1);
	array<double> values;
	for (int i = 0; i < 2 * steps + 1; i++) values.add(initialMomentum(-L / 2 + i * dx));
	return values;
}

array<double> psinotSample() {
	double dx = L / (2 * steps + 1);
	double coeff = sqrt(3 * m / (4 * (M_SQRT_TWO)));
	array<double> values;
	for (int i = 0; i < 2 * steps; i++) {
		double sh = sech(m / M_SQRT_TWO * (-L / 2 + i * dx));
		values.add(coeff * sh * sh);
	}
	return values;
}

array<double> dpsinotSample() {
	double dx = L / (2 * steps + 1);
	double coeff = -sqrt(3 * mmm / (2 * (M_SQRT_TWO)));
	array<double> values;
	for (int i = 0; i < 2 * steps; i++) {
		double inner = m / M_SQRT_TWO * (-L / 2 + i * dx);
		double sh = sech(inner);
		double th = tanh(inner);
		values.add(coeff * sh * sh * th);
	}
	return values;
}


void do_math() {
	


}


int main() {
	
	//generate and gather all of our initial data
	array<double> pvals       = PSample(.05);
	array<double> initPosVals = initialPosSample();
	array<double> initMomVals = initialMomentumSample();
	array<double> psinotVals  = psinotSample();
	array<double> dpsinotVals = dpsinotSample();
	
	//init deshi
	deshi::init();
	Render::UseDefaultViewProjMatrix();

	//start main loop
	TIMER_START(t_d); TIMER_START(t_f);
	while (!deshi::shouldClose()) {
		DeshiImGui::NewFrame();                    //place imgui calls after this
		Memory::Update();
		DeshTime->Update();
		DeshWindow->Update();
		DeshInput->Update();
		DeshConsole->Update(); Console2::Update();

		{//debug area
			

		}

		UI::Update();
		Render::Update();                          //place imgui calls before this
		DeshTime->frameTime = TIMER_END(t_f); TIMER_RESET(t_f);
	}



	//cleanup deshi
	deshi::cleanup();
}