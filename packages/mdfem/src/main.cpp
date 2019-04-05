#include "fem.h"

int main(int argc, char** argv) { 
    FEMFrame* fem = new FEMFrame(); 
    fem->n_bdy_nodes = 8;

    fem->Alloc();
    fem->eval();
    INFO_Printf("fem->n_bdy_nodes = %d\n", fem->n_bdy_nodes);
    return 0;
}
