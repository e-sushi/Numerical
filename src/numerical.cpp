#include "../kigu/common.h"
#include "../kigu/array.h"

//TODO maybe make a compiler flag for compiling with/without deshi
#define DESHI_DISABLE_IMGUI
#include "core/assets.h"
#include "core/commands.h"
#include "core/console.h"
#include "core/input.h"
#include "core/logger.h"
#include "core/memory.h"
#include "core/imgui.h"
#include "core/renderer.h"
#include "core/storage.h"
#include "core/threading.h"
#include "core/time.h"
#include "core/ui.h"
#include "core/window.h"
#include "core/io.h"

#include "domath.cpp"
#include "showmath.cpp"

#include <fstream>
#include <limits>

int main() {
	deshi::init();

	while(!deshi::shouldClose()){
		do_math();
		show_math();
	}

	deshi::cleanup();
}