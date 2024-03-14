
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <Kokkos_Core.hpp>

#include "matar.h"
#include "driver.h"


int main(int argc, char* argv[])
{


    if (argc < 2)
    {
        std::cout << "Fierro requires a mesh file input" << std::endl;
    }

    int num_solvers = 2;


    Kokkos::initialize();


    Driver driver(argv[1]);

    driver.initialize(num_solvers);
    driver.setup();
    driver.run();
    driver.finalize();

    Kokkos::finalize();

    std::cout << "**** End of main **** " << std::endl;
    return 0;
}
