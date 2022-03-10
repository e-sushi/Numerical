#include "kigu/common.h"
#include "kigu/array.h"
#include "core/io.h"

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
		DeshWindow->Update();
		DeshTime->Update();
		DeshInput->Update();
		DeshiImGui::NewFrame();

		show_math(do_math());

		Render::Update();
		memory_clear_temp();
		
	}

	deshi::cleanup();
}