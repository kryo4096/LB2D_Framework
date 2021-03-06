/** 
 *  @file
 *  @author Fabian Bösch
 *  @brief simulation
 */

#ifndef LB_SIMULATION_HPP_INCLUDED
#define LB_SIMULATION_HPP_INCLUDED

#define _USE_MATH_DEFINES

#include <cmath>

#include "H_root.hpp"
#include "lattice.hpp"
#include "timer.hpp"

#include <sstream>
#include <random>

#include <omp.h>
#include <cassert>
#include <complex>

namespace lb {

/**
 *  @brief Simulation class implementing LB
 * 
 *  This class holds a lattice as member (see @ref simulation::l) and 
 *  carries out the simulation steps on top of it. The main methods of 
 *  this class are @ref simulation::advect() and 
 *  @ref simulation::collide().
 */

    enum CollisionType {
        LBGK,
        KBC,
        NONE,
    };

    const std::string CollisionTypeNames[]{
            "LBGK",
            "KBC",
            "None"
    };

    std::ostream &operator<<(std::ostream &os, CollisionType type) {
        return os << CollisionTypeNames[type];
    }

    class simulation {
    public:
        lattice l;                 ///< lattice
        std::vector<int> shift;    ///< amount of nodes to shift each population in data structure during advection
        scalar_t Re;       ///< Reynolds number
        scalar_t Vmax;     ///< mean flow velocity
        scalar_t visc;     ///< viscosity
        scalar_t beta;     ///< LB parameter beta
        unsigned int time;         ///< simulation time
        bool file_output;          ///< flag whether to write files
        unsigned int output_freq;  ///< file output frequency
        unsigned int output_index; ///< index for file naming
        CollisionType collisionType;
        bool periodic = true;
		bool regularized = false;
		bool calculate_wall_force = false;
	    scalar_t wall_force_sum_x = 0;
	    scalar_t wall_force_sum_y = 0;

    public: // ctor


        /**
         *  @brief Construct from domain size and flow parameters
         *  @param[in] nx    extent in x direction
         *  @param[in] ny    extent in y direction
         *  @param[in] Re   Reynolds number
         *  @param[in] Vmax mean flow velocity
         */
        simulation(unsigned int nx, unsigned int ny, scalar_t Re, scalar_t Vmax,
                   CollisionType collisionType = CollisionType::LBGK)
                : l(nx, ny),
                  shift(velocity_set().size),
                  Re(Re),
                  Vmax(Vmax),
                  visc(Vmax * nx / Re),
                  beta(1 / (2 * visc / (velocity_set().cs * velocity_set().cs) + 1)),
                  time(0),
                  file_output(false), // set to true if you want to write files
                  output_freq(100),
                  output_index(0),
                  collisionType(collisionType) {
            // define amount to shift populations for advection
            for (unsigned int i = 0; i < velocity_set().size; ++i) {
                shift[i] = l.index(velocity_set().c[0][i], velocity_set().c[1][i]) - l.index(0, 0);
            }
        }

		void calculate_beta() {
			visc = Vmax * l.nx / Re;
			beta = 1 / (2 * visc / (velocity_set().cs * velocity_set().cs) + 1);
		}

        void debug_advection() {
            for (int i = 0; i < l.nx; i++) {
                for (int j = 0; j < l.ny; j++) {

                    auto &node = l.get_node(i, j);

                    if (std::abs<int>(i - l.nx / 2) < 10 && std::abs<int>(j - l.ny / 2) < 10)
                        for (int k = 0; k < 9; k++) {
                            node.f(k) = 1.0;
                        }
                }
            }
        }

        inline bool should_collide(int i, int j) {
            return periodic || i < l.nx - 1;
        }

        template<typename Func>
        void initialize(Func initializer) {

            for (int i = 0; i < l.nx; i++) {
                for (int j = 0; j < l.ny; j++) {

                    auto[u, v, rho, wall] = initializer(*this, i, j);

                    scalar_t f_eq[9];
                    auto node = l.get_node(i, j);

                    if (wall) {
                        l.add_wall(node.coord, node.coord);
                    }

                    v9::f_eq(f_eq, rho, u, v);

                    for (int p = 0; p < 9; p++) {

                        node.f(p) = f_eq[p];
                    }

                    node.u() = u;
                    node.v() = v;
                    node.rho() = rho;

                }
            }
        }

        /**
         *  @brief advect the populations
         *
         *  Include periodic boundary conditions here also
         */
        void advect() {
            // Velocity 1 : (1, 0)
#pragma omp for
            for (int j = 0; j < l.ny; j++) {
                for (int i = l.nx - 1; i >= 0; i--) {
                    int idx = l.index(i, j);
                    l.f[1][idx + shift[1]] = l.f[1][idx];
                }
            }

            // Velocity 2 : (0, 1)
#pragma omp for
            for (int i = 0; i < l.nx; i++) {
                for (int j = l.ny - 1; j >= 0; j--) {
                    int idx = l.index(i, j);
                    l.f[2][idx + shift[2]] = l.f[2][idx];
                }

            }

            // Velocity 3 : (-1, 0)
#pragma omp for
            for (int j = 0; j < l.ny; j++) {
                for (int i = 0; i < l.nx; i++) {
                    int idx = l.index(i, j);
                    l.f[3][idx + shift[3]] = l.f[3][idx];
                }
            }

            // Velocity 4: (0, -1)
#pragma omp for
            for (int i = 0; i < l.nx; i++) {
                for (int j = 0; j < l.ny; j++) {
                    int idx = l.index(i, j);
                    l.f[4][idx + shift[4]] = l.f[4][idx];
                }
            }

            //velocity 5: (1,1)
#pragma omp for
            for (int d = 0; d < l.nx; d++) {
                for (int i = d, j = l.ny - 1; i >= 0 && j >= 0; i--, j--) {
                    l.f[5][l.index(i + 1, j + 1)] = l.f[5][l.index(i, j)];
                }
            }

#pragma omp for
            for (int d = 0; d < l.ny - 1; d++) {
                for (int i = l.nx - 1, j = d; i >= 0 && j >= 0; i--, j--) {
                    l.f[5][l.index(i + 1, j + 1)] = l.f[5][l.index(i, j)];
                }
            }

            //velocity 6: (-1,1)
#pragma omp for
            for (int d = 0; d < l.ny - 1; d++) {
                for (int i = 0, j = d; i < l.nx && j >= 0; i++, j--) {
                    l.f[6][l.index(i, j) + shift[6]] = l.f[6][l.index(i, j)];
                }
            }

#pragma omp for
            for (int d = 0; d < l.nx; d++) {
                for (int i = d, j = l.ny - 1; i < l.nx && j >= 0; i++, j--) {
                    l.f[6][l.index(i, j) + shift[6]] = l.f[6][l.index(i, j)];
                }
            }

            // velocity 7: (-1, -1)
#pragma omp for
            for (int d = 1; d < l.ny; d++) {
                for (int i = 0, j = d; i < l.nx && j < l.ny; i++, j++) {
                    l.f[7][l.index(i, j) + shift[7]] = l.f[7][l.index(i, j)];
                }
            }

#pragma omp for
            for (int d = 0; d < l.nx; d++) {
                for (int i = d, j = 0; i < l.nx && j < l.ny; i++, j++) {
                    l.f[7][l.index(i, j) + shift[7]] = l.f[7][l.index(i, j)];
                }
            }

            // velocity 8: (1, -1)
#pragma omp for
            for (int d = 0; d < l.ny; d++) {
                for (int i = l.nx - 1, j = d; i >= 0 && j < l.ny; i--, j++) {
                    l.f[8][l.index(i, j) + shift[8]] = l.f[8][l.index(i, j)];
                }
            }

#pragma omp for
            for (int d = 0; d < l.nx - 1; d++) {
                for (int i = d, j = 0; i >= 0 && j < l.ny; i--, j++) {
                    l.f[8][l.index(i, j) + shift[8]] = l.f[8][l.index(i, j)];
                }
            }
            if (periodic) {
                auto apply_buffers = [&](int i, int j) {
                    int buf_index = l.index(i, j);

                    int new_i = (i + l.nx) % l.nx;
                    int new_j = (j + l.ny) % l.ny;

                    int new_index = l.index(new_i, new_j);

                    for (int p = 1; p < velocity_set().size; p++) {
                        int source_i = i - velocity_set().c[0][p];
                        int source_j = j - velocity_set().c[1][p];

                        if (source_i >= 0 && source_j >= 0 && source_i < l.nx && source_j < l.ny) {
                            l.f[p][new_index] = l.f[p][buf_index];
                        }
                    }
                };

#pragma omp for
                for (int i = 0; i < l.nx + 1; i++) apply_buffers(i, -1);
#pragma omp for
                for (int i = 0; i < l.nx + 1; i++) apply_buffers(i, l.ny);
#pragma omp for
                for (int j = 0; j < l.ny + 2; j++) apply_buffers(-1, j - 1);
#pragma omp for
                for (int j = 0; j < l.ny; j++) apply_buffers(l.nx, j);
            }
        }


        /**  @brief apply wall boundary conditions */
        void boundary_conditions() {

			scalar_t force_x = 0;
			scalar_t force_y = 0;

            for (auto &node: l.wall_nodes) {
                for (int k = 0; k < 9; k++) {
                    auto opp = velocity_set().opposites[k];
                    if (!l.get_node(node.index + shift[opp]).has_flag_property("wall")) {
                        l.f[opp][node.index + shift[opp]] = node.f(k);

						force_x += 2 * node.f(k) * v9::c[0][k];
	                    force_y += 2 * node.f(k) * v9::c[1][k];
                    }
                }
            }

			if(time % 10000 > 2000) {
				wall_force_sum_x += force_x;
				wall_force_sum_y += force_y;
			}

			if (calculate_wall_force && time % 10000 == 0) {
				std::cout << Vmax / v9::cs << ", " << wall_force_sum_x / 8000 << ", " << wall_force_sum_y / 8000 << "\n";
				Vmax += 0.01;
				calculate_beta();
				wall_force_sum_x = 0;
				wall_force_sum_y = 0;
			}

            if (!periodic) {
                for (int j = 0; j < (int) l.ny; j++) {
                    scalar_t u = Vmax;

                    scalar_t v = 0;
                    scalar_t rho = 1.0;

                    scalar_t f_eq[9];
                    velocity_set().f_eq(f_eq, rho, u, v);

                    auto node = l.get_node(0, j);
                    node.u() = u;
                    node.v() = v;
                    node.rho() = rho;

                    for (int p = 0; p < 9; p++) {
                        node.f(p) = f_eq[p];
                    }
                }

                for (int i = 0; i <= (int) l.nx; i++) {
                    l.get_node(i, 0).f(2) = l.get_node(i, -1).f(4);
                    l.get_node(i - 1, 0).f(5) = l.get_node(i, -1).f(8);
                    l.get_node(i + 1, 0).f(6) = l.get_node(i, -1).f(7);
                    l.get_node(i, l.ny - 1).f(4) = l.get_node(i, l.ny).f(2);
                    l.get_node(i - 1, l.ny - 1).f(8) = l.get_node(i, l.ny).f(5);
                    l.get_node(i + 1, l.ny - 1).f(7) = l.get_node(i, l.ny).f(6);
                }
            }
        }

        /** @brief collide the populations */
        void collide_lbgk() {
#pragma omp for
            for (int j = 0; j < static_cast<int>(l.ny); ++j) {
                for (int i = 0; i < static_cast<int>(l.nx); ++i) {

                    if (!should_collide(i, j)) continue;

                    scalar_t beta_ = beta;

                    if (!periodic && i > l.nx - 30) {
                        beta_ = 0.5;
                    }

                    int idx = l.index(i, j);

                    scalar_t f[9];
                    for (int k = 0; k < 9; k++) {
                        f[k] = l.f[k][idx];
                    }

                    scalar_t rho = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8];
                    scalar_t u = (f[1] - f[3] + f[5] - f[6] - f[7] + f[8]) / rho;
                    scalar_t v = (f[2] - f[4] + f[5] + f[6] - f[7] - f[8]) / rho;

                    scalar_t x_root = std::sqrt(1 + 3 * u * u);
                    scalar_t y_root = std::sqrt(1 + 3 * v * v);

                    scalar_t A = rho * (2 - x_root) * (2 - y_root);
                    scalar_t BX = (2 * u + x_root) / (1 - u);
                    scalar_t BY = (2 * v + y_root) / (1 - v);

                    l.f[0][idx] = f[0] + 2 * beta_ * (16.0 / 36.0 * A - f[0]);
                    l.f[1][idx] = f[1] + 2 * beta_ * (4.0 / 36.0 * A * BX - f[1]);
                    l.f[2][idx] = f[2] + 2 * beta_ * (4.0 / 36.0 * A * BY - f[2]);
                    l.f[3][idx] = f[3] + 2 * beta_ * (4.0 / 36.0 * A / BX - f[3]);
                    l.f[4][idx] = f[4] + 2 * beta_ * (4.0 / 36.0 * A / BY - f[4]);
                    l.f[5][idx] = f[5] + 2 * beta_ * (1.0 / 36.0 * A * BX * BY - f[5]);
                    l.f[6][idx] = f[6] + 2 * beta_ * (1.0 / 36.0 * A / BX * BY - f[6]);
                    l.f[7][idx] = f[7] + 2 * beta_ * (1.0 / 36.0 * A / BX / BY - f[7]);
                    l.f[8][idx] = f[8] + 2 * beta_ * (1.0 / 36.0 * A * BX / BY - f[8]);

                    l.u[idx] = u;
                    l.v[idx] = v;
                    l.rho[idx] = rho;
                }
            }

        }

        void collide_kbc() {
            scalar_t f[9];
            scalar_t f_eq[9];
            scalar_t delta_s[9];
            scalar_t delta_h[9];

#pragma omp for
            for (int j = 0; j < static_cast<int>(l.ny); ++j) {
                for (int i = 0; i < static_cast<int>(l.nx); ++i) {

                    if (!should_collide(i, j)) continue;

                    scalar_t beta_ = beta;

	                if (!periodic && i > l.nx - 30) {
		                beta_ = 0.5;
	                }

                    int idx = l.index(i, j);

                    for (int k = 0; k < 9; k++) {
                        f[k] = l.f[k][idx];
                    }

                    scalar_t rho = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8];
                    scalar_t u = (f[1] - f[3] + f[5] - f[6] - f[7] + f[8]) / rho;
                    scalar_t v = (f[2] - f[4] + f[5] + f[6] - f[7] - f[8]) / rho;

                    scalar_t x_root = std::sqrt(1 + 3 * u * u);
                    scalar_t y_root = std::sqrt(1 + 3 * v * v);

                    scalar_t A = rho * (2 - x_root) * (2 - y_root);
                    scalar_t BX = (2 * u + x_root) / (1 - u);
                    scalar_t BY = (2 * v + y_root) / (1 - v);

                    f_eq[0] = 16.0 / 36.0 * A;
                    f_eq[1] = 4.0 / 36.0 * A * BX;
                    f_eq[2] = 4.0 / 36.0 * A * BY;
                    f_eq[3] = 4.0 / 36.0 * A / BX;
                    f_eq[4] = 4.0 / 36.0 * A / BY;
                    f_eq[5] = 1.0 / 36.0 * A * BX * BY;
                    f_eq[6] = 1.0 / 36.0 * A / BX * BY;
                    f_eq[7] = 1.0 / 36.0 * A / BX / BY;
                    f_eq[8] = 1.0 / 36.0 * A * BX / BY;

                    scalar_t Pi = f[5] - f[6] + f[7] - f[8];
                    scalar_t N = f[1] - f[2] + f[3] - f[4];
                    scalar_t Pi_eq = f_eq[5] - f_eq[6] + f_eq[7] - f_eq[8];
                    scalar_t N_eq = f_eq[1] - f_eq[2] + f_eq[3] - f_eq[4];

                    delta_s[0] = 0;
                    delta_s[1] = 0.25f * (N - N_eq);
                    delta_s[2] = -0.25f * (N - N_eq);
                    delta_s[3] = 0.25f * (N - N_eq);
                    delta_s[4] = -0.25f * (N - N_eq);
                    delta_s[5] = 0.25f * (Pi - Pi_eq);
                    delta_s[6] = -0.25f * (Pi - Pi_eq);
                    delta_s[7] = 0.25f * (Pi - Pi_eq);
                    delta_s[8] = -0.25f * (Pi - Pi_eq);


                    scalar_t ds_dh = 0;
                    scalar_t dh_dh = 0;

                    for (int k = 0; k < 9; k++) {
                        delta_h[k] = f[k] - f_eq[k] - delta_s[k];
                        ds_dh += delta_h[k] * delta_s[k] / f_eq[k];
                        dh_dh += delta_h[k] * delta_h[k] / f_eq[k];
                    }

                    if (dh_dh == 0) dh_dh = 1e-20;

                    scalar_t gamma = 1 - (2 * beta_ - 1) * ds_dh / dh_dh;

					if(regularized) {
						gamma = 1;
					}

                    for (int k = 0; k < 9; k++) {
                        l.f[k][idx] -= beta_ * 2 * delta_s[k] + gamma * delta_h[k];
                    }

                    l.u[idx] = u;
                    l.v[idx] = v;
                    l.rho[idx] = gamma;
                }
            }
        }

        void calculate_moments() {
#pragma omp for
            for (int i = 0; i < l.nx; i++) {
                for (int j = 0; j < l.ny; j++) {
                    int idx = l.index(i, j);

                    scalar_t f[9];

                    for (int k = 0; k < 9; k++) {
                        f[k] = l.f[k][idx];
                    }

                    scalar_t rho = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8];
                    l.rho[idx] = rho;
                    l.u[idx] = (f[1] - f[3] + f[5] - f[6] - f[7] + f[8]) / rho;
                    l.v[idx] = (f[2] - f[4] + f[5] + f[6] - f[7] - f[8]) / rho;
                }
            }
        }

        /** @brief LB step */
        void step() {


#pragma omp parallel
            advect();
            boundary_conditions();

            switch (collisionType) {
                case CollisionType::LBGK:
#pragma omp parallel
                    collide_lbgk();
                    break;
                case CollisionType::KBC:
#pragma omp parallel
                    collide_kbc();
                    break;
                case CollisionType::NONE:
#pragma omp parallel
                    calculate_moments();
                    break;
            }


            // file io
            if (file_output && (((time + 1) % output_freq) == 0 || time == 0)) {
                write_fields();
                ++output_index;
            }

            ++time;
        }

        /** write macroscopic variables to ascii file */
        void write_fields() {
            std::stringstream fns;
            fns << "output/data_" << std::setfill('0') << std::setw(4) << output_index << ".txt";
            l.write_fields(fns.str());
        }

        /** print to output stream */
        friend std::ostream &operator<<(std::ostream &os, const simulation &sim) {
            os << "simulation parameters\n"
               << "---------------------\n";
            os << "collision type: " << sim.collisionType << "\n";
            os << "domain: " << sim.l.nx << " x " << sim.l.ny << "\n";
            os << "Re:     " << sim.Re << "\n";
            os << "Vmax:   " << sim.Vmax << " mach : "<< sim.Vmax / v9::cs <<"\n";
            os << "visc:   " << sim.visc << "\n";
            os << "beta:   " << sim.beta << "\n";
            return os;
        }


    };

} // lb

#endif // LB_SIMULATION_HPP_INCLUDED
