//
// Created by jonas on 12/22/21.
//

#ifndef LB2D_RUNNER_HPP
#define LB2D_RUNNER_HPP

#include "simulation.hpp"
#include "visualization.hpp"
#include <omp.h>
#include <filesystem>
#include <tuple>

namespace lb::runner {
    void cylinder_flow() {
        auto sim = simulation(600, 50, 1e5, 0.1, lb::CollisionType::KBC);

        sim.initialize([&](const simulation &sim, int i, int j) -> std::tuple<scalar_t, scalar_t, scalar_t, bool> {

            scalar_t D = sim.l.nx / 40.0;

            bool wall = std::pow(i - 10 * D, 2) + std::pow(j - (int) sim.l.ny / 2, 2) < std::pow(0.5 * D, 2);

            scalar_t u = std::max(sim.Vmax * (9.5 * D - i) / (9.5 * D), 0.0);

            return std::tuple{u, 0.0, 1.0, wall};
        });

        sim.periodic = false;
        std::cout << sim << std::endl;

        printf("beta = %.20f\n", sim.beta);

        try {
            lb::visualization::initialize(&sim, 0, nullptr);
            lb::visualization::get_instance().run();
        }
        catch (std::runtime_error &e) {
            std::cerr << e.what() << std::endl;
        }
    }


    void shear_layer() {
        auto sim = simulation(512, 512, 3000, 0.05, lb::CollisionType::KBC);

        scalar_t kappa = 80.0;
        scalar_t delta = 0.05;

        sim.initialize([&](const simulation &sim, int i, int j) -> std::tuple<scalar_t, scalar_t, scalar_t, bool> {
            scalar_t u;
            if (j <= sim.l.ny / 2) {
                u = sim.Vmax * tanh(kappa * (j / (scalar_t) sim.l.ny - 0.25));
            } else {
                u = sim.Vmax * tanh(kappa * (0.75 - j / (scalar_t) sim.l.ny));
            }

            scalar_t v = sim.Vmax * delta * sin(2 * M_PI * ((scalar_t) i / sim.l.nx + 0.25));
            scalar_t rho = 1.0;

            return std::tuple{u, v, rho, false};
        });

        std::cout << sim << std::endl;

        printf("beta = %.20f\n", sim.beta);

        try {
            lb::visualization::initialize(&sim, 0, nullptr);
            lb::visualization::get_instance().run();
        }
        catch (std::runtime_error &e) {
            std::cerr << e.what() << std::endl;
        }
    }

    void debug() {
        auto sim = simulation(100, 112, 0, 0, CollisionType::NONE);
        sim.debug_advection();

        try {
            lb::visualization::initialize(&sim, 0, nullptr);
            lb::visualization::get_instance().run();
        }
        catch (std::runtime_error &e) {
            std::cerr << e.what() << std::endl;
        }
    }

    void convergence_test(int count = 11, int start_size = 128, scalar_t log_step_size = 0.5, int iterations = 50,
                          scalar_t Vmax = 0.05, scalar_t Re = 3e4,
                          lb::CollisionType collision_type = CollisionType::LBGK) {

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

            auto sim = simulation(L, L, Re, Vmax, collision_type);

            scalar_t Kx = 2 * M_PI / L;
            scalar_t Ky = 2 * M_PI / L;
            scalar_t Ksqr = Kx * Kx + Ky * Ky;
            scalar_t nu = sim.visc;
            scalar_t l2error = 0;
            scalar_t dA = 1.0 / L / L;
            scalar_t mach = sim.Vmax / velocity_set().cs;

            sim.initialize([&](const simulation &sim, int i, int j) -> std::tuple<scalar_t, scalar_t, scalar_t, bool> {
                scalar_t u = -sim.Vmax * std::cos(Kx * i) * std::sin(Ky * j);
                scalar_t v = sim.Vmax * std::cos(Ky * j) * std::sin(Kx * i);
                scalar_t rho = 1 - mach * mach / (2 * Ksqr) *
                                   (Ky * Ky * std::cos(2 * Kx * i) + Kx * Kx * std::cos(2 * Ky * j));

                return std::tuple{u, v, rho, false};
            });

            int nx = sim.l.nx;
            int ny = sim.l.ny;

            for (int t = 0; t < iterations; t++) {
                sim.step();
            }

            scalar_t linferror = 0;

            for (int i = 0; i < nx; i++) {
                for (int j = 0; j < ny; j++) {
                    scalar_t u = -Vmax * std::cos(Kx * i) * std::sin(Ky * j) * std::exp(-nu * Ksqr * sim.time);
                    scalar_t v = Vmax * std::cos(Ky * j) * std::sin(Kx * i) * std::exp(-nu * Ksqr * sim.time);

                    auto &node = sim.l.get_node(i, j);

                    scalar_t error_u = u - node.u();
                    scalar_t error_v = v - node.v();

                    scalar_t error = (error_u * error_u + error_v * error_v);
                    l2error += error * dA;
                    linferror = std::max(sqrt(error), linferror);
                }
            }

            l2error = std::sqrt(l2error);
            l2errors.push_back(l2error);
            linferrors.push_back(linferror);
            Ls.push_back(L);

            if (n > 0) {
                std::cout
                        << std::setw(20) << L << ", "
                        << std::setw(20) << l2error << ", "
                        << std::setw(20) << -log(l2errors[n] / l2errors[n - 1]) / log(Ls[n] / (scalar_t) Ls[n - 1])
                        << ", "
                        << std::setw(20) << linferror << ", "
                        << std::setw(20) << -log(linferrors[n] / linferrors[n - 1]) / log(Ls[n] / (scalar_t) Ls[n - 1]);

                printf(",\t %.20f\n", sim.beta);
            }
        }
    }

}

#endif //LB2D_RUNNER_HPP
