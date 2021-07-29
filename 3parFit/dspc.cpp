#include "eclCovmatAlgorithm.h"

using namespace Belle2;

int main(){

    eclCovmatAlgorithm * eclFg = new eclCovmatAlgorithm();
    eclFg->calibrate();

    if(eclFg) delete eclFg;
    return 0;
}
