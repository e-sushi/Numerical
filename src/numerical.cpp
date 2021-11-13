//// deshi includes ////
#include "defines.h"
#include "deshi.h"
#include "core/memory.h"
#include "utils/string.h"
#include "utils/array.h"
#include "utils/map.h"
#include "math/Math.h"
#include "core/logging.h"

#include <fstream>

struct Point {
	//for the new equations 
	// x is rho
	// v is cursive pi
	// P is cursive phi

	//in the paper K is momentum g is position

	double pos = 0, mom = 0;
	//iteration vars, so i dont have to keep calculating them for each point when i already have
	double posi[4] = { 0 }, momi[4] = { 0 };

	//copy operator cuz c-arrays werent being copied 
	void operator = (const Point& rhs) {
		pos = rhs.pos; mom = rhs.mom;
		memcpy(posi, rhs.posi, 4 * sizeof(double));
		memcpy(momi, rhs.momi, 4 * sizeof(double));
	}
};

//
// Global Variables
//

double dt ;//= 0.01; // time step value
u32    Nt ;//= 50;   // how many time steps you're going to do 
double Bi ;//= 0;    // initial velocity as a fraction of speed of light
double Bf ;//= 0.9;  // final velocity
double L  ;//= 10;   // length of world
u32    Np ;//= 5;    // 2 * Np + 1 is the number of grid points

double m   ;//= 1;
double mm  ;//= m*m;
double mmm ;//= m*m*m;
		
double g   ;//= 0.2;
double gg  ;//= g*g;
double ggg ;//= g*g*g;

u32    NN = 0;    // num of time steps taken

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

//non linear term helper func
inline double NLT(double curr) {
	return -mm * curr - 4 * gg * (curr * curr * curr);
}

//
//Initial Condition Functions
//

//double P(double t) {
//	if (t < t0 - delta / 2) return -k / 2;
//	if (t > t0 + delta / 2) return  k / 2;
//	double inside = (t - t0) / (sqrt((t - t0 + delta / 2) * (t0 - t + delta / 2)));
//	return k / 2 * tanh(inside);
//}
//
//
//array<double> PSample(double dt) {
//	array<double> values;
//	for (int i = 0; i < Nt; i++) values.add(P(i * dt));
//	return values;
//}
//
//
//double initialPos(double rho) {
//	double M0 = mmm / (3 * M_SQRT_TWO * gg);
//	double gamma = sqrt(1 + (k / 2) * (k / 2) / (M0 * M0));
//	return m / (2 * g) * tanh(gamma * m * rho / M_SQRT_TWO);
//}
//
//array<double> initialPosSample() {
//	double dx = L / (2 * Nt + 1);
//	array<double> values;
//	for (int i = 0; i < 2 * Nt + 1; i++) values.add(initialPos(-L / 2 + i * dx));
//	return values;
//}
//
//double IIRiemannSum(double hk, double M0, double gamma) {
//	double sum = 0;
//	double dx = 1 / (100 * m);
//	for (double i = -L * 100; i < L * 100; i += dx) {
//		double s0 = sech(i);
//		double s1 = sech(gamma * i);
//		sum += s0 * s0 * s1 * s1 * dx;
//	}
//	return sum;
//}
//
//double initialMomentum(double rho) {
//	double hk = k / 2;
//	double M0 = mmm / (3 * M_SQRT_TWO * gg);
//	double gamma = sqrt(1 + hk * hk / (M0 * M0));
//	double II = IIRiemannSum(hk, M0, gamma);
//	double sech0 = sech(sqrt(1 + hk * hk / (M0 * M0) * m * rho / M_SQRT_TWO));
//	double sech1 = sech(m * rho / M_SQRT_TWO);
//	return 3 * g * k / (4 * m) * ( sech0 * sech0 - (3 * II * sech1 * sech1) / 4);
//}
//
//array<double> initialMomentumSample() {
//	double dx = L / (2 * Nt + 1);
//	array<double> values;
//	for (int i = 0; i < 2 * Nt + 1; i++) values.add(initialMomentum(-L / 2 + i * dx));
//	return values;
//}
//
//array<double> psinotSample() {
//	double dx = L / (2 * Nt + 1);
//	double coeff = sqrt(3 * m / (4 * (M_SQRT_TWO)));
//	array<double> values;
//	for (int i = 0; i < 2 * Nt; i++) {
//		double sh = sech(m / M_SQRT_TWO * (-L / 2 + i * dx));
//		values.add(coeff * sh * sh);
//	}
//	return values;
//}
//
//array<double> dpsinotSample() {
//	double dx = L / (2 * Nt + 1);
//	double coeff = -sqrt(3 * mmm / (2 * (M_SQRT_TWO)));
//	array<double> values;
//	for (int i = 0; i < 2 * Nt; i++) {
//		double inner = m / M_SQRT_TWO * (-L / 2 + i * dx);
//		double sh = sech(inner);
//		double th = tanh(inner);
//		values.add(coeff * sh * sh * th);
//	}
//	return values;
//}

//initial data arrays
array<double> PList;
array<double> PsiZeroList;
array<double> dPsiZeroList;
array<double> ccList;
array<double> d2ccList;
array<Point>  initialPointList;

//iterated arrays
array<Point>  pointIter; 
array<Point>  nextPointIter; 

//system info arrays
array<double> Hbar; //total energy of system 
array<double> Hsc;  //energy of just the soliton

//visualization arrays
array<float> posshow;

//4th order Runge Kutta scheme
void do_math() {
	double dx = L / (2 * Np + 1);

	//Runge Kutta iterations 
	for (u32 ITERidx = 0; ITERidx < pointIter.count; ITERidx++) {
		Point& curp = pointIter[ITERidx];

		for (u32 RKidx = 0; RKidx < 4; RKidx++) {
			
			array<Point>  dummyPoints;
			array<Point>  dummyDPoints;
			array<Point>  dummyD2PosPoints;
			array<double> dummyNLTList;
			double RS1; //Riemann Sum
			double RS2; //Riemann Sum
			double RS3; //Riemann Sum
			double RS4; //Riemann Sum
			RS1 = RS2 = RS3 = RS4 = 0;

			{//set initial lists according to what RK iteration we are on
				if (!RKidx) {
					dummyPoints = pointIter;
				}
				else if (RKidx < 4) {
					//if we are in between first and last iterations then we must add half of the last
					for (Point& p : pointIter) {
						Point nu;
						nu.pos + 0.5 * p.posi[RKidx];
						nu.mom + 0.5 * p.momi[RKidx];

						dummyPoints.add(nu);
					}
				}
				else {
					for (Point& p : pointIter) {
						Point nu;
						nu.pos + p.posi[RKidx];
						nu.mom + p.momi[RKidx];

						dummyPoints.add(nu);
					}
				}
			}


			{//fill in our dummy derivative using SLAC derivative 
				for (u32 j = 0; j < 2 * Np + 1; j++) {
					Point sum1, sum2;
					for (u32 kk = 0; kk < j; kk++) {
						sum1.pos += ccList[j - 1 - kk] * dummyPoints[kk].pos;
						sum1.mom += ccList[j - 1 - kk] * dummyPoints[kk].mom;
					}
					for (u32 kk = j; kk < 2 * Np; kk++) {
						sum2.pos += ccList[2 * Np - 1 - kk + j] * dummyPoints[kk].pos;
						sum2.mom += ccList[2 * Np - 1 - kk + j] * dummyPoints[kk].mom;
					}

					sum1.pos -= sum2.pos;
					sum1.mom -= sum2.mom;

					dummyDPoints.add(sum1);
				}
			}

			{//fill in our dummy second derivative
				for (u32 j = 0; j < 2 * Np + 1; j++) {
					Point sum1, sum2;
					for (u32 kk = 0; kk < j; kk++) {
						sum1.pos += d2ccList[j - 1 - kk] * dummyPoints[kk].pos;
					}
					for (u32 kk = j; kk < 2 * Np; kk++) {
						sum2.pos += d2ccList[2 * Np - 1 - kk + j] * dummyPoints[kk].pos;
					}

					sum1.pos -= sum2.pos;
					sum1.pos += (-(M_PI * M_PI) / (3 * dx * dx) * (1 - 1 / (2 * Np + 1)) * (1 - 2 / (2 * Np + 1))) * dummyPoints[j].pos;

					dummyD2PosPoints.add(sum1);
				}
			}

			{//Non Linear Term List
				for (u32 j = 0; j < 2 * Np + 1; j++) {
					dummyNLTList.add(NLT(dummyPoints[j].pos));
				}
			}

			{//Riemann Sums
				for (u32 j = 0; j < 2 * Np; j++) {
					RS1 += dx * ((dummyPoints[j].mom + dummyPoints[j + 1].mom) * (dummyDPoints[j].pos + dummyDPoints[j + 1].pos) / 2);
					RS2 += dx * ((PsiZeroList[j] + PsiZeroList[j + 1]) * (dummyDPoints[j].pos + dummyDPoints[j + 1].pos) / 2);
					RS3 += dx * ((PsiZeroList[j] + PsiZeroList[j + 1]) * ((dummyD2PosPoints[j].pos - dummyNLTList[j]) + (dummyD2PosPoints[j + 1].pos - dummyNLTList[j + 1])) / 2);
					RS4 += dx * ((PsiZeroList[j] + PsiZeroList[j + 1]) * (dummyDPoints[j].mom + dummyDPoints[j + 1].mom) / 2);
				}

				RS1 += dx * (dummyPoints.last->mom - dummyPoints[0].mom) * (dummyDPoints.last->pos - dummyDPoints[0].pos) / 2;
				RS2 += dx * (*PsiZeroList.last - PsiZeroList[0]) * (dummyDPoints.last->pos - dummyDPoints[0].pos) / 2;
				RS3 += dx * (*PsiZeroList.last - PsiZeroList[0]) * ((dummyD2PosPoints.last->pos - *dummyNLTList.last) + (dummyD2PosPoints[0].pos - dummyNLTList[0])) / 2;
				RS4 += dx * (*PsiZeroList.last - PsiZeroList[0]) * (dummyDPoints.last->mom - dummyDPoints[0].mom) / 2;
			}


			{//Hbar
				double hbar;

				double t0n = PList[NN + 1] + RS1;
				double term0 = t0n * t0n / (2 * RS2 * RS2);

				double term1 = 0;
				for (u32 j = 0; j < 2 * Np; j++) {
					double lt1 = (dummyPoints[j].mom + dummyPoints[j + 1].mom)   / 2;
					double lt2 = (dummyDPoints[j].pos + dummyDPoints[j + 1].pos) / 2;
					double lt3 = (dummyPoints[j].pos + dummyPoints[j + 1].pos)   / 2;
					double lt4 = (gg * lt3 * lt3 - 0.25 * mm);
					double lt5 = lt4 * lt4 / gg;

					term1 += dx * lt1 / 2 + lt2 / 2 + lt5;
				}

				double lt1 = (dummyPoints.last->mom  - dummyPoints[0].mom)  / 2;
				double lt2 = (dummyDPoints.last->pos - dummyDPoints[0].pos) / 2;
				double lt3 = (dummyPoints.last->pos  - dummyPoints[0].pos)  / 2;
				double lt4 = (gg * lt3 * lt3 - 0.25 * mm);
				double lt5 = lt4 * lt4 / gg;

				term1 += dx * lt1 / 2 + lt2 / 2 + lt5;

				hbar = term0 + term1;

				Hbar.add(hbar);
			}

			{//Hsc
				double hbar;

				double t0n = PList[NN + 1] * PList[NN + 1] + RS1 * RS1;
				double term0 = t0n * t0n / (2 * RS2 * RS2);

				double term1 = 0;
				for (u32 j = 0; j < 2 * Np; j++) {
					double lt1 = (dummyPoints[j].mom + dummyPoints[j + 1].mom)   / 2;
					double lt2 = (dummyDPoints[j].pos + dummyDPoints[j + 1].pos) / 2;
					double lt3 = (dummyPoints[j].pos + dummyPoints[j + 1].pos)   / 2;
					double lt4 = (gg * lt3 * lt3 - 0.25 * mm);
					double lt5 = lt4 * lt4 / gg;

					term1 += dx * lt1 / 2 + lt2 / 2 + lt5;
				}

				double lt1 = (dummyPoints.last->mom - dummyPoints[0].mom)   / 2;
				double lt2 = (dummyDPoints.last->pos - dummyDPoints[0].pos) / 2;
				double lt3 = (dummyPoints.last->pos - dummyPoints[0].pos)   / 2;
				double lt4 = (gg * lt3 * lt3 - 0.25 * mm);
				double lt5 = lt4 * lt4 / gg;

				term1 += dx * lt1 / 2 + lt2 / 2 + lt5;

				hbar = term0 + term1;

				Hsc.add(hbar);
			}


			for (u32 j = 0; j < 2 * Np + 1; j++) {
				pointIter[RKidx].pos = dt * ((PList[NN + 1] + RS1) / (RS2 * RS2) * (dummyDPoints[j].pos - RS2 * PsiZeroList[j]) + dummyPoints[j].mom);
				pointIter[RKidx].mom = 
					dt * 
					((PList[NN + 1] + RS1) / (RS2 * RS2) * (dummyDPoints[j].mom - RS4 * PsiZeroList[j]) - 
					(((PList[NN + 1] + RS1) * ((PList[NN + 1] + RS1)) / (RS2 * RS2 * RS2) * dPsiZeroList[j] +
										(dummyD2PosPoints[j].pos - dummyNLTList[j]) - RS3 * PsiZeroList[j])));
			}
		}
	}

	posshow.clear();
	nextPointIter.clear();

	//final Runge Kutta iteration
	for(Point& curp : pointIter) {
		Point nu;
		nu.pos = curp.pos + (1 / 6) * (curp.posi[0] + curp.posi[1] + curp.posi[2] + curp.posi[3]);
		nu.mom = curp.mom + (1 / 6) * (curp.momi[0] + curp.momi[1] + curp.momi[2] + curp.momi[3]);
		nextPointIter.add(nu);
		posshow.add(nu.pos);
	}

}


int main() {
	
	//generate and gather all of our initial data
	//pvals       = PSample(.05);
	//initPosVals = initialPosSample();
	//initMomVals = initialMomentumSample();
	//psinotVals  = psinotSample();
	//dpsinotVals = dpsinotSample();
	
	//init deshi
	deshi::init();
	Render::UseDefaultViewProjMatrix();

	//gather initial data
	using namespace std;
	ifstream file;
	file.open("C:/Users/sushi/Dropbox/BHRR/initialdata/FSEinitialdata.csv");
	std::string line;


	forI(8) {
		getline(file, line, (i == 8 - 1) ? '\n' : ',');
		switch(i) {
			case 0: dt = stod(line); continue; break;
			case 1: Nt = stod(line); continue; break;
			case 2: Bi = stod(line); continue; break;
			case 3: Bf = stod(line); continue; break;
			case 4: L  = stod(line); continue; break;
			case 5: Np = stod(line); continue; break;
			case 6: {
				m   = stod(line);
				mm  = m*m;
				mmm = m*m*m;
			}break;
			case 7: {
				g   = stod(line);
				gg  = g*g;
				ggg = g*g*g;
			}break;
		}
	}

	//getline(file, line, ',');
	forI(6) {
		switch (i) {
			case 0: { //Plist
				forX(idx, Nt) {
					getline(file, line, (idx == Nt - 1) ? '\n' : ',');
					PList.add(stod(line));
				}
			}break;
			case 1: { //psizero
				forX(idx, 2 * Np + 1) {
					getline(file, line, (idx == 2 * Np) ? '\n' : ',');
					PsiZeroList.add(stod(line));
				}
			}break;
			case 2: { //dpsizero
				forX(idx, 2 * Np + 1) {
					getline(file, line, (idx == 2 * Np) ? '\n' : ',');
					dPsiZeroList.add(stod(line));
				}
			}break;
			case 3: {
				forX(idx, 2 * Np) {
					getline(file, line, (idx == 2 * Np - 1) ? '\n' : ',');
					ccList.add(stod(line));
				}
			}break;
			case 4: {
				forX(idx, 2 * Np) {
					getline(file, line, (idx == 2 * Np - 1) ? '\n' : ',');
					d2ccList.add(stod(line));
				}
			}break;
			case 5: {
				forX(idx, 2 * Np + 1) {
					Point nu;
					getline(file, line, ',');
					nu.pos = stod(line);
					getline(file, line, ',');
					nu.mom = stod(line);
					pointIter.add(nu);
					//idx++;

				}
			}break;
		}
	}

	
	file.close();

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
			static int latch = 0;

			//if (!latch || DeshInput->LMousePressed()){
				do_math();

				forI(nextPointIter.count) {
					Log("", "last ", TOSTRING(pointIter[i].pos), " next ", TOSTRING(nextPointIter[i].pos));
				}
				pointIter = nextPointIter;

				NN++;
				latch++;
		//	}
			ImGui::SetNextWindowPos(ImVec2(0, 0));
			ImGui::SetNextWindowSize(ImVec2(DeshWindow->width, DeshWindow->height));
			ImGui::Begin("graph", 0, ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize);
			ImGui::PlotLines("wave", &posshow[0], posshow.count, 0, (const char*)0, -10, 10, ImVec2(800, 800));
			ImGui::SameLine();
			ImGui::Text(TOSTRING(latch).str);
			ImGui::End();
			
		}

		UI::Update();
		Render::Update();                          //place imgui calls before this
		DeshTime->frameTime = TIMER_END(t_f); TIMER_RESET(t_f);
	}



	//cleanup deshi
	deshi::cleanup();
}