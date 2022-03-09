//this file contains all of the variables and functions related to displaying our solutions calculated in domath.cpp
// entry point: show_math()

void show_math(array<Point> in) {

    ImGui::SetNextWindowPos({0,0});
    ImGui::SetNextWindowSize({DeshWinSize.x, DeshWinSize.y});
    ImGui::Begin("dataplot");{
        ImGui::PlotLines("FSEPosSol", 
        [](void* data, int idx){
            return f32(((Point*)data+idx)->pos);
        }, 
        in.data, in.count);

        ImGui::PlotLines("FSEMomSol", 
        [](void* data, int idx){
            return f32(((Point*)data+idx)->mom);
        }, 
        in.data, in.count);
    }ImGui::End();

}