#include "simulation.hpp"
#include "visualization.hpp"
#include "runner.hpp"

#include <omp.h>
#include <filesystem>

//#define FLOAT_DEBUG

void show_help();

#ifdef FLOAT_DEBUG
#include <fenv.h>
#endif


int main(int argc, char *argv[]) {
    omp_set_num_threads(6);

#ifdef FLOAT_DEBUG
    feenableexcept(FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW | FE_INVALID);
#endif

    std::string mode;

    if (argc < 2) {
        std::cout << "you didn't specify a program mode." << std::endl;
        mode = "help";
    } else {
        mode = std::string(argv[1]);
    }

    if (mode == "cylinder_flow") lb::runner::cylinder_flow();
    else if (mode == "shear_layer") lb::runner::shear_layer();
    else if (mode == "convergence_lbgk")
        lb::runner::convergence_test(17, 128, 0.25, 20, 0.05, 3e4, lb::CollisionType::LBGK);
    else if (mode == "convergence_kbc")
        lb::runner::convergence_test(17, 128, 0.25, 20, 0.05, 3e4, lb::CollisionType::KBC);
    else if (mode == "debug") lb::runner::debug();
    else if (mode == "help") show_help();
    else {
        std::cerr << "ERROR: " << mode << " is not a valid mode." << std::endl;
        return -1;
    }

    return 0;
}

void show_help() {
    std::cout << "Available Modes: \n"
              << "cylinder_flow : Simulates flow around a cylinder using the KBC-D collision model.\n"
              << "shear_layer : Simulates a periodic shear layer using the KBC-D collision model. \n"
              << "convergence_lbgk : taylor-green-vortex benchmark for lbgk. \n"
              << "convergence_kbc : taylor-green-vortex benchmark for kbc \n";
}

