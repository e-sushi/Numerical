#include "kigu/common.h"
#include "kigu/array.h"
#include "core/io.h"
#include "core/graphing.h"
#include "core/memory.h"

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

	array<array<Point>> frames;
	s32 frame_show = 0;
	frames.add(do_math());

	while(!deshi::shouldClose()){
		DeshWindow->Update();
		DeshTime->Update();
		DeshInput->Update();
		DeshiImGui::NewFrame();
		
		if(DeshInput->KeyDown(Key::SPACE))
			frames.add(do_math());
		if(DeshInput->KeyDown(Key::RIGHT)) 
			frame_show = Min(frame_show+1, frames.count-1);
		if(DeshInput->KeyDown(Key::LEFT)) 
			frame_show = Max(frame_show-1, 0);


		show_math(frames);

		UI::Update();
		Render::Update();
		memory_clear_temp();
		
	}

	deshi::cleanup();
}