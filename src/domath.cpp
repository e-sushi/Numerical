//this file contains all the functions and variables involved in actually calcuating our solutions.
//we use 4th order Runge Kutta to approximate solutions to our equations
// entry point: do_math()
//        init: init_math()

//NOTE throughout the code we have comments explaing what each thing is doing and these comments reference things in our paper
//     however as the paper gets updated equation numbers may be wrong due to LaTeX's automatic numbering.
 

//@vars
//NOTE currently (3/2/2022) these varaibles are gathered from the initial data file we generate using Mathematica
// ideally we will generate initial data using C++ in the future, so these may need to be moved somewhere else later
f64 dt; //= 0.01; // time step value
u32 Nt; //= 50;   // how many time steps you're going to do 
f64 Bi; //= 0;    // initial velocity as a fraction of speed of light
f64 Bf; //= 0.9;  // final velocity
f64 L;  //= 10;   // length of world
u32 Np; //= 5;    // 2 * Np is the number of grid points in the world

f64 m;   //= 1;
f64 mm;  //= m*m;
f64 mmm; //= m*m*m;

f64 g;   //= 0.2;
f64 gg;  //= g*g;
f64 ggg; //= g*g*g;

u32 NN = 0;    // num of time steps taken

//initial data arrays
array<f64>   PList;
array<f64>   PsiZeroList;
array<f64>   dPsiZeroList;
array<f64>   ccList;
array<f64>   d2ccList;
array<Point> initialPointList;

//iterated arrays
array<Point> points; //zeta in the paper subsection 3.3
array<Point> nextPoints; //staging array

//staging arrays for tempdata required throughout RK
//they are always allocated to the same size as the world
//and their data is simply overwritten
array<Point> pointsStage;
array<Point> pointsDStage; //first derivative of point array
array<Point> pointsD2Stage; //NOTE the only thing we actually take the second derivative of is momentum, so maybe just store that instead of the whole point
array<f64>   NLTStage;

//@functions
//helpers
typedef f64 (*MathFunc)(f64 p);
//i dont think this will actually be used anywhere
f64 riemannSum(f64 start, f64 end, f64 dx, MathFunc func) {
	f64 sum = 0;
	for (f64 i = start; i < end; i += dx)
		sum += func(i) * dx;
	return sum;
}

f64 sech(f64 p) { return 1 / cosh(p); }
//non linear term helper func
inline f64 NLT(f64 curr) { return -mm * curr + 4 * gg * (curr * curr * curr); }

array<Point> do_math() {
	f64 dx = L / (2 *  Np);

	for(u32 rki = 0; rki < 4; rki++){
		
		{//setup point fields based on what iteration were on
			
			//Here we are defining the staging script F (\FF) that appears in the paper,
			//equation (3.25) in subsection 3.3. Referring to the paper, \zeta is the phase space vector,
			//whose top component is the momentum list and bottom component is the position list.

			//At the first step of Runge Kutta we want to take the argument of \FF to be the actual solution to the
			//previous time step, which is called \zeta^(0) in the paper. The result is then called \zeta^(1)
			if(!rki){
				pointsStage=points;
			}
			//In the second step we take the argument of \FF to be \zeta^(0) + 0.5 * \zeta^(1)
			//and the result is called \zeta^(2)
			//In the third step we take the argument of \FF to be \zeta^(0) + 0.5 * \zeta^(2)
			//and the result is called \zeta^(3)
			else if(rki<3){
				forI(points.count){
					pointsStage[i].pos = points[i].pos + 0.5 * points[i].posi[rki-1];
					pointsStage[i].mom = points[i].mom + 0.5 * points[i].momi[rki-1];
				}
			}
			//In the fourth iteration we take the argument of \FF to be \zeta^(0) + \zeta^(3) 
			//and the result is called \zeta^(4)
			else{
				forI(points.count){
					pointsStage[i].pos = points[i].pos + points[i].posi[rki-1];
					pointsStage[i].mom = points[i].mom + points[i].momi[rki-1];
				}
			}

		}

	
		{//fill out derivatives using SLAC 
		
			//Here we compute the spacial derivatives of the staging points list using
			//the SLAC derivative. This is equation (3.12) of subsection 3.1 in the paper for
			//the action of the SLAC derivative in terms of the c coefficients and equation (3.11)
			//for the value cs (or equation (3.17) in the case of the second derivative).

			//the first derivative
			for(u32 j = 0; j < 2*Np; j++){
				Point sum1,sum2;
				for(u32 kk = 0; kk < j; kk++){
					sum1.pos += ccList[j-1-kk]*pointsStage[kk].pos;
					sum1.mom += ccList[j-1-kk]*pointsStage[kk].mom;
				}
				for(u32 kk = j; kk < 2 * Np - 1; kk++){
					sum2.pos += ccList[2*Np-2-kk+j]*pointsStage[kk+1].pos;
					sum2.mom += ccList[2*Np-2-kk+j]*pointsStage[kk+1].mom;
				}
				sum1.pos-=sum2.pos;
				sum1.mom-=sum2.mom;

				pointsDStage[j] = sum1;
			}

			//the second derivative
			for(u32 j = 0; j < 2*Np; j++){
				Point sum1,sum2;
				for(u32 kk = 0; kk < j; kk++){
					sum1.pos += d2ccList[j-1-kk]*pointsStage[kk].pos;
				}
				for(u32 kk = j; kk < 2 * Np - 1; kk++){
					sum2.pos += d2ccList[2*Np-2-kk+j]*pointsStage[kk+1].pos;
				}
				sum1.pos+=(-(π * π)/(3.0*(4.0*Np*Np*dx*dx)) * (4.0*Np*Np-1)) * pointsStage[j].pos;
				sum1.pos-=sum2.pos;

				pointsD2Stage[j] = sum1;
			}

		}

		{//Non linear term list

			for(u32 j = 0; j < 2*Np; j++){
				NLTStage[j] = NLT(pointsStage[j].pos);
			}

		}

		f64 RS1=0,RS2=0,RS3=0,RS4=0; //Riemann Sums 

		{//evaluate Riemann Sums

			for(u32 j = 0; j < 2*Np; j++){
				RS1+=pointsStage[0].mom*pointsDStage[j].pos;
				RS2+=PsiZeroList[j]*pointsDStage[j].pos;
				RS3+=PsiZeroList[j]*(pointsD2Stage[j].pos-NLTStage[j]);
				RS4+=PsiZeroList[j]*pointsDStage[j].mom;
			}
			RS1*=dx; RS2*=dx; RS3*=dx; RS4*=dx;

		}

		//TODO Hsc and Hbar

		{//do final RK step

			//Here we evaluate the discretized Forced Soliton Equation. We apply
			//\FF to the current staged points list and the result is stored in the 
			//momi and posi arrays on the points. 

			for(u32 j = 0; j < 2*Np; j++){
				pointsStage[j].posi[rki] = dt * ((PList[NN+1]+RS1)/(RS2*RS2)*(pointsDStage[j].pos - RS2 * PsiZeroList[j]) + pointsStage[j].mom);
				pointsStage[j].momi[rki] = 
					dt * 
					((PList[NN + 1] + RS1) / (RS2 * RS2) * (pointsDStage[j].mom - RS4 * PsiZeroList[j]) - 
					(((PList[NN + 1] + RS1) * ((PList[NN + 1] + RS1)) / (RS2 * RS2 * RS2) * dPsiZeroList[j] +
										(pointsD2Stage[j].pos - NLTStage[j]) - RS3 * PsiZeroList[j])));
			}


		}
	}

	//Once we complete the 4 RK iterations we use them to construct the next time step.

	for(Point& curp : points){
		Point nu;
		curp.pos = curp.pos + (1.0/6.0)*(curp.posi[0]+2*curp.posi[1]+2*curp.posi[2]+curp.posi[3]);
		curp.mom = curp.mom + (1.0/6.0)*(curp.momi[0]+2*curp.momi[1]+2*curp.momi[2]+curp.momi[3]);

	}

	return points;

}

void init_math() {
	//TODO replace std io with File from deshi 
	using namespace std;
	ifstream file;
	file.open("src/initialdata/FSEinitialdata.csv");
	std::string line;
	
	//gather global constants used to make initial data
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

	//Plist
	forX(idx, Nt) {
		getline(file, line, (idx == Nt - 1) ? '\n' : ',');
		PList.add(stod(line));
	}
			
	//psizero
	forX(idx, 2 * Np) {
		getline(file, line, (idx == 2 * Np - 1) ? '\n' : ',');
		PsiZeroList.add(stod(line));
	}
			
	//dpsizero
	forX(idx, 2 * Np) {
		getline(file, line, (idx == 2 * Np - 1) ? '\n' : ',');
		dPsiZeroList.add(stod(line));
	}
			
	forX(idx, 2 * Np - 1) {
		getline(file, line, (idx == 2 * Np - 2) ? '\n' : ',');
		ccList.add(stod(line));
	}
			
	forX(idx, 2 * Np - 1) {
		getline(file, line, (idx == 2 * Np - 2) ? '\n' : ',');
		d2ccList.add(stod(line));
	}
			
	forX(idx, 2 * Np) {
		Point nu;
		getline(file, line, ',');
		nu.pos = stod(line);
		getline(file, line, ',');
		nu.mom = stod(line);
		points.add(nu);
	}

	pointsStage.resize(2*Np);
	pointsDStage.resize(2*Np);
	pointsD2Stage.resize(2*Np);

	file.close();
}
