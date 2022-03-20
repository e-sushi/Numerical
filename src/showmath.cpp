//this file contains all of the variables and functions related to displaying our solutions calculated in domath.cpp
// entry point: show_math()

void show_math(array<Point> in, u32 frame, u32 max_frames) {
    persist b32 init = 0;
	persist const u32 res = 300;
	persist Graph g;



	if(!in.count)return;
	if(!init){
		g.xAxisLabel = cstr("x");
		g.yAxisLabel = cstr("y");
		init=1;
	}
	else{
		UI::Begin("graphe", vec2::ONE, vec2(600,600), UIWindowFlags_NoScroll);
		u32 res = Min(DeshWinSize.x, UI::GetWindow()->width - UI::GetStyle().windowMargins.x*2);
		array<vec2g> temp;
        temp.resize(in.count);
        g.data={temp.data,temp.count};

        UI::Text(toStr("frame: ", frame, " / ", max_frames).str);
	
		forI(in.count){
			temp[i].x = i*L/(2*Np);
            temp[i].y = in[i].pos;
		}
		
		
		draw_graph(g, UI::GetWindowRemainingSpace());
		UIItem* gr = UI::GetLastItem();
		

		static vec2 mp;
		static vec2 gcp;
		if(UI::IsLastItemHovered() && DeshInput->LMousePressed()){
			UI::SetPreventInputs();
			mp = DeshInput->mousePos;
			gcp = g.cameraPosition;
		}
		if(mp!=vec2::ONE*FLT_MAX && DeshInput->LMouseDown()){
			g.cameraPosition = gcp - (DeshInput->mousePos - mp) /  (g.dimensions_per_unit_length*g.aspect_ratio);
		}
		if(DeshInput->LMouseReleased()){
			UI::SetAllowInputs();
			mp=vec2::ONE*FLT_MAX;
		}
		g.cameraZoom -= 0.2*g.cameraZoom*DeshInput->scrollY;
		
		UI::End();

	}

}