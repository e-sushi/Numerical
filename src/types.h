#include "kigu/common.h"

//lattice point type
//represents position and momentum at some point in the world
//keeps track of it's RK4 iterations
struct Point {
	//in the paper K is momentum g is position

	f64 pos = 0, mom = 0;
	//iteration vars, so i dont have to keep calculating them for each point when i already have
	f64 posi[4] = { 0 }, momi[4] = { 0 };

	//copy operator cuz c-arrays werent being copied 
	void operator = (const Point& rhs) {
		pos = rhs.pos; mom = rhs.mom;
		memcpy(posi, rhs.posi, 4 * sizeof(f64));
		memcpy(momi, rhs.momi, 4 * sizeof(f64));
	}
};