#include "simulation.hpp"
#include "visualization.hpp"
#include <omp.h>
#include <filesystem>



namespace lb {
    void visual_sim(int argc, char *argv[]) {
        auto sim = simulation(512, 256, 1e5, 0.05, lb::CollisionType::KBC);
        sim.doubly_periodic_shear_layer();
        std::cout << sim << std::endl;

        printf("beta = %.20f\n", sim.beta);

        try {
            lb::visualization::initialize(&sim, argc, argv);
            lb::visualization::get_instance().run();
        }
        catch (std::runtime_error& e) {
            std::cerr << e.what() << std::endl;
        }
    }

    void convergence_test(int count, int start_size, scalar_t log_step_size, int iterations, scalar_t Vmax, scalar_t Re, lb::CollisionType collision_type) {

        std::cout << "Simulation parameters:" << std::endl;
        std::cout << "- time steps: " << iterations << std::endl;
        std::cout << "- mean velocity: " << Vmax << std::endl;
        std::cout << "- Reynolds number: " << Re << std::endl;
        std::cout << "- Collision type: " << collision_type << std::endl;
        std::cout << "\n";

        std::cout
            << std::setw(20) << "domain size L"
            << std::setw(20) << "L2 Error"
            << std::setw(20) << "L2 rate"
            << std::setw(20) << "LInf Error"
            << std::setw(20) << "LInf rate"
            << std::setw(20) << "beta"
            << std::endl;

        std::vector<scalar_t> l2errors;
        std::vector<scalar_t> linferrors;
        std::vector<scalar_t> Ls;

        for (int n = 0; n < count; n++) {

            int L = start_size * powf(2, n * log_step_size);

            auto *sim = new simulation(L, L, Re, Vmax, collision_type);
            sim->taylor_green();

            int nx = sim->l.nx;
            int ny = sim->l.ny;

            scalar_t Kx = 2 * M_PI / nx;
            scalar_t Ky = 2 * M_PI / ny;
            scalar_t Ksqr = Kx * Kx + Ky * Ky;
            scalar_t nu = sim->visc;
            scalar_t l2error = 0;
            scalar_t dA = 1.0 / L / L;

            for(int t = 0; t < iterations; t++) {
                sim->step();
            }

            scalar_t linferror = 0;

            for(int i = 0; i < nx; i++) {
                for(int j = 0; j < ny; j++) {
                    scalar_t u = - Vmax * cos(Kx * i) * sin(Ky * j) * exp(-nu * Ksqr * iterations);
                    scalar_t v = Vmax * cos(Ky * j) * sin(Kx * i) * exp(-nu * Ksqr * iterations);

                    auto& node = sim->l.get_node(i, j);

                    scalar_t error_u = abs(u - node.u()) / L; // velocity measured relative to domain size
                    scalar_t error_v = abs(v - node.v()) / L;

                    scalar_t error = (error_u * error_u + error_v * error_v);
                    l2error += error * dA;
                    linferror = std::max(error_u, std::max(error_v, linferror));
                }
            }

            l2error = std::sqrt(l2error);
            l2errors.push_back(l2error);
            linferrors.push_back(linferror);
            Ls.push_back(L);

            if(n > 0) {
                std::cout
                    << std::setw(20) << L
                    << std::setw(20) << l2error
                    << std::setw(20) << -log(l2errors[n] / l2errors[n - 1]) / log(Ls[n] / (scalar_t) Ls[n - 1])
                    << std::setw(20) << linferror
                    << std::setw(20) << -log(linferrors[n] / linferrors[n - 1]) / log(Ls[n] / (scalar_t) Ls[n - 1]);

                printf(",\t %.20f\n", sim->beta);
            }
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
    if (mode == "taylor_green_convergence") lb::convergence_test(11, 128, 0.5, 100, 0.1, 1e20, lb::CollisionType::LBGK);
    else {
        std::cerr << "ERROR: " << mode << " is not a valid mode." << std::endl;
        return -1;
    }

	return 0;
}

