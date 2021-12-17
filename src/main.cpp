
#include "simulation.hpp"
#include "visualization.hpp"
#include <omp.h>
#include <filesystem>

namespace lb {
    void visual_sim(int argc, char *argv[]) {
        auto *sim = new simulation(256, 256, 1000.0, 0.1);
        sim->doubly_periodic_shear_layer();
        std::cout << *sim << std::endl;

        try {
            lb::visualization::initialize(sim, argc, argv);
            lb::visualization::get_instance().run();
        }
        catch (std::runtime_error& e) {
            std::cerr << e.what() << std::endl;
        }
    }

    void convergence_test(int count, int start_size, float_type log_step_size, int iterations, float_type Vmax, float_type Re) {

        std::cout << "Simulation parameters:" << std::endl;

        std::cout << "- time steps: " << iterations << std::endl;
        std::cout << "- mean velocity: " << Vmax << std::endl;
        std::cout << "- Reynolds number: " << Re << std::endl;

        std::cout << "\n";

        std::cout << "L \tError (L2) \tConvergence Rate" << std::endl;

        std::vector<float_type> errors;
        std::vector<float_type> Ls;

        for (int n = 0; n < count; n++) {

            int L = start_size * powf(2, n * log_step_size);

            auto *sim = new simulation(L, L, Re, Vmax);
            sim->taylor_green();

            int nx = sim->l.nx;
            int ny = sim->l.ny;

            float_type Kx = 2 * M_PI / nx;
            float_type Ky = 2 * M_PI / ny;
            float_type Ksqr = Kx * Kx + Ky * Ky;
            float_type nu = sim->visc;
            float_type l2error = 0;
            float_type dA = 1.0 / L / L;

            for(int t = 0; t < iterations; t++) {
                sim->step();
                for(int i = 0; i < nx; i++) {
                    for(int j = 0; j < ny; j++) {
                        float_type u = - Vmax * cos(Kx * i) * sin(Ky * j) * exp(-nu*Ksqr * t);
                        float_type v = Vmax * cos(Ky * j) * sin(Kx * i) * exp(-nu*Ksqr * t);

                        auto& node = sim->l.get_node(i, j);

                        float_type error_u = abs(u - node.u());
                        float_type error_v = abs(v - node.v());
                        l2error += (error_u * error_u + error_v * error_v) * dA;
                    }
                }
            }


            std::cout << L << ", \t" << l2error;

            errors.push_back(l2error);
            Ls.push_back(L);

            if(n > 0) {
                std::cout << ",\t" << -log(errors[n] / errors[n-1]) / log(Ls[n] / (float_type) Ls[n-1]); // log-log graph slope
            }

            std::cout << std::endl;
        }
    }
}

int main(int argc, char *argv[])
{
	omp_set_num_threads(6);

    std::string mode;

    if (argc < 2) {
        mode = "visual";
    } else {
        mode = std::string(argv[1]);
    }

    if (mode == "visual") lb::visual_sim(argc, argv);
    if (mode == "taylor_green_convergence") lb::convergence_test(25, 128, 0.125, 100, 0.05, 1e4);
    else {
        std::cerr << "ERROR: " << mode << " is not a valid mode." << std::endl;
        return -1;
    }

	return 0;
}

