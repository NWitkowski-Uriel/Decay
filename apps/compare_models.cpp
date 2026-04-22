#include <iostream>
#include "Constants.h"
#include "DecayFunctions.h"
#include "RadialGrid.h"

int main() {
    RadialGrid::instance().setNumPoints(200);

    std::cout << "Comparing models at y=0 for proton" << std::endl;

    for (double mt = mp + 10; mt < mp + 2000; mt += 50) {
        double d = dN_dmt_proton_from_Delta_Dirac(mt, 0.0);
        double bw = dN_dmt_proton_from_Delta_BW(mt, 0.0);
        double ps = dN_dmt_proton_from_Delta_PS(mt, 0.0);

        double r_bw = (d > 0) ? bw / d : 0.0;
        double r_ps = (d > 0) ? ps / d : 0.0;

        std::cout << mt << "  BW/Dirac=" << r_bw << "  PS/Dirac=" << r_ps << std::endl;
    }

    return 0;
}
