#include "sw.h"

int main(int argc, char** argv) { 
    SWFrame* sw = new SWFrame(); 

    sw->Alloc();
    sw->eval();
    return 0;
}
