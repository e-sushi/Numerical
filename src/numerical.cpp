#include "kigu/common.h"
#include "kigu/array.h"

#include <fstream>

//TODO maybe make a compiler flag for compiling with/without deshi
#include "deshi.h"

#include "types.h"
#include "domath.cpp"
#include "showmath.cpp"

#include <fstream>
#include <limits>

int main() {
	deshi::init();
	init_math();

	while(!deshi::shouldClose()){
		show_math(do_math());
	}

	deshi::cleanup();
}