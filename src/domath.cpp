//this file contains all the functions and variables involved in actually calcuating our solutions.
//we use 4th order Runge Kutta to approximate solutions to our equations
// entry point: do_math()


//@vars
//NOTE currently (3/2/2022) these varaibles are gathered from the initial data file we generate using Mathematica
// ideally we will generate initial data using C++ in the future, so these may need to be moved somewhere else later
f64 dt; //= 0.01; // time step value
u32 Nt; //= 50;   // how many time steps you're going to do 
f64 Bi; //= 0;    // initial velocity as a fraction of speed of light
f64 Bf; //= 0.9;  // final velocity
f64 L;  //= 10;   // length of world
u32 Np; //= 5;    // 2 * Np + 1 is the number of grid points

f64 m;   //= 1;
f64 mm;  //= m*m;
f64 mmm; //= m*m*m;

f64 g;   //= 0.2;
f64 gg;  //= g*g;
f64 ggg; //= g*g*g;

u32 NN = 0;    // num of time steps taken

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

void do_math() {

}
