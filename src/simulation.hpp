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
    KBC
};

class simulation {
	public:
	    lattice l;                 ///< lattice
	    std::vector<int> shift;    ///< amount of nodes to shift each population in data structure during advection
	    const float_type Re;       ///< Reynolds number
	    const float_type Vmax;     ///< mean flow velocity
	    const float_type visc;     ///< viscosity
	    const float_type beta;     ///< LB parameter beta
	    unsigned int time;         ///< simulation time
	    bool file_output;          ///< flag whether to write files
	    unsigned int output_freq;  ///< file output frequency
	    unsigned int output_index; ///< index for file naming
	    CollisionType collisionType;

    public: // ctor


        /**
         *  @brief Construct from domain size and flow parameters
         *  @param[in] nx    extent in x direction
         *  @param[in] ny    extent in y direction
         *  @param[in] Re   Reynolds number
         *  @param[in] Vmax mean flow velocity
         */
        simulation(unsigned int nx, unsigned int ny, float_type Re, float_type Vmax, CollisionType collisionType=CollisionType::LBGK)
                : l(nx, ny),
                  shift(velocity_set().size),
                  Re(Re),
                  Vmax(Vmax),
                  visc(Vmax * ny / Re),
                  beta(1 / (2 * visc / (velocity_set().cs * velocity_set().cs) + 1)),
                  time(0),
                  file_output(false), // set to true if you want to write files
                  output_freq(100),
                  output_index(0),
                  collisionType(collisionType)
                  {
            // define amount to shift populations for advection
            for (unsigned int i = 0; i < velocity_set().size; ++i) {
                shift[i] = l.index(velocity_set().c[0][i], velocity_set().c[1][i]) - l.index(0, 0);
            }
        }

        /**
         *  @brief Initialize the flow field
         *
         *  Initialization includes defining initial density, velocity and
         *  populations. You can use Taylor-Green vortex flow conditions.
         */
        void taylor_green() {

            float_type mach = Vmax / velocity_set().cs;

            float_type Kx = 2 * M_PI / l.nx;
            float_type Ky = 2 * M_PI / l.ny;

            float_type Ksqr = Kx * Kx + Ky * Ky;

            for(int i = 0; i < l.nx; i++) {
                for(int j = 0; j < l.nx; j++) {
                    float_type u = - Vmax * cos(Kx * i) * sin(Ky * j);
                    float_type v = Vmax * cos(Ky * j) * sin(Kx * i);
                    float_type rho = 1 - mach * mach / (2 * Ksqr) * (Ky*Ky*cos(2 * Kx * i) + Kx*Kx*sin(2 * Ky * j));

                    float_type f_eq[9];
                    velocity_set().f_eq(f_eq, rho, u, v);

                    auto& node = l.get_node(i, j);

                    for(int p = 0; p < 9; p++) {
                        node.f(p) = f_eq[p];
                    }

                    node.u() = u;
                    node.v() = v;
                    node.rho() = rho;
                }
            }
        }

        void doubly_periodic_shear_layer() {
            float_type kappa = 80.0;
            float_type delta = 0.05;

            float_type Kx = 2 * M_PI / l.nx;


            for(int i = 0; i < l.nx; i++) {
                for(int j = 0; j < l.nx; j++) {

                    float_type u;
                    if (j <= l.ny / 2) {
                        u = Vmax * tanh(kappa * (j / (float_type) l.ny - 0.25));
                    } else {
                        u = Vmax * tanh(kappa * (0.75 - j / (float_type) l.ny));
                    }

                    float_type v = Vmax * delta * sin(2 * M_PI * ((float_type) i / l.nx + 0.25));
                    float_type rho = 1.0;

                    float_type f_eq[9];
                    velocity_set().f_eq(f_eq, rho, u, v);

                    auto& node = l.get_node(i, j);

                    for(int p = 0; p < 9; p++) {
                        node.f(p) = f_eq[p];
                    }

                    node.u() = u;
                    node.v() = v;
                    node.rho() = rho;
                }
            }
        }

        void initialize() {

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
            for (int d = 0; d < l.ny; d++) {
                for (int i = d, j = l.nx - 1; i >= 0 && j >= 0; i--, j--) {
                    l.f[5][l.index(i, j) + shift[5]] = l.f[5][l.index(i, j)];
                }
            }

			#pragma omp for
            for (int d = 0; d < l.nx - 1; d++) {
                for (int i = l.ny-1, j = d; i >= 0 && j >= 0; i--, j--) {
                    l.f[5][l.index(i, j) + shift[5]] = l.f[5][l.index(i, j)];
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

			auto apply_buffers = [&](int i, int j) {
				int buf_index = l.index(i, j);

				int new_i = (i + l.nx) % l.nx;
				int new_j = (j + l.ny) % l.ny;

				int new_index = l.index(new_i, new_j);

				for(int p = 1; p < velocity_set().size; p++) {
					int source_i = i - velocity_set().c[0][p];
					int source_j = j - velocity_set().c[1][p];

					if(source_i >= 0 && source_j >= 0 && source_i < l.nx && source_j < l.ny) {
						l.f[p][new_index] = l.f[p][buf_index];
					}
				}
			};

            #pragma omp barrier

			#pragma omp for
			for(int i = 0; i < l.nx + 1; i++) apply_buffers(i, -1);
			#pragma omp for
	        for(int i = 0; i < l.nx + 1; i++) apply_buffers(i, l.ny);
			#pragma omp for
	        for(int j = 0; j < l.ny + 2; j++) apply_buffers(-1, j-1);
			#pragma omp for
	        for(int j = 0; j < l.ny; j++) apply_buffers(l.nx, j);
        }


        /**  @brief apply wall boundary conditions */
        void wall_bc() {
            #pragma omp parallel for
            for (int i = 0; i < l.wall_nodes.size(); ++i) {
                // **************************
                // * fill in your code here *
                // **************************
            }
        }

        /** @brief collide the populations */
        void collide_lbgk() {

            static auto c = velocity_set().c;

            #pragma omp for schedule(static)
            for (int j = 0; j < static_cast<int>(l.ny); ++j) {
                for (int i = 0; i < static_cast<int>(l.nx); ++i) {
                    int idx = l.index(i, j);

                    float_type rho = 0;
                    float_type u = 0;
                    float_type v = 0;
                    float_type f[9];

                    for (int k = 0; k < velocity_set().size; k++) {
                        f[k] = l.f[k][idx];
                        rho += f[k];
                        u += f[k] * (float_type) c[0][k];
                        v += f[k] * (float_type) c[1][k];
                    }

                    u /= rho;
                    v /= rho;

                    l.u[idx] = u;
                    l.v[idx] = v;
                    l.rho[idx] = rho;

	                float_type f_eq[9];
	                velocity_set().f_eq(f_eq, rho, u, v);

					for(int k = 0; k < velocity_set().size; k++) {
						l.f[k][idx] = f[k] + 2 * beta * (f_eq[k] - f[k]);
					}
                }
            }

        }

        void collide_kbc() {
            #pragma omp for
            for (int j = 0; j < static_cast<int>(l.ny); ++j) {
                for (int i = 0; i < static_cast<int>(l.nx); ++i) {
                    int idx = l.index(i, j);

                    float_type rho = 0;
                    float_type u = 0;
                    float_type v = 0;
                    float_type f[9];

                    for (int k = 0; k < velocity_set().size; k++) {
                        f[k] = l.f[k][idx];
                        rho += f[k];
                        u += f[k] * (float_type) c[0][k];
                        v += f[k] * (float_type) c[1][k];
                    }

                    u /= rho;
                    v /= rho;

                    l.u[idx] = u;
                    l.v[idx] = v;
                    l.rho[idx] = rho;

                    float_type f_eq[9];
                    velocity_set().f_eq(f_eq, rho, u, v);

                    for(int k = 0; k < velocity_set().size; k++) {
                        l.f[k][idx] = f[k] + 2 * beta * (f_eq[k] - f[k]);
                    }


                }
            }
        }

        /** @brief LB step */
        void step() {

            Timer timer;

            timer.start();
            #pragma omp parallel
            advect();
            timer.stop();

            std::cout << "Advection took " << timer.elapsedMilliseconds() << " ms.\n";

            timer.start();
            switch ( collisionType ){
                case CollisionType::LBGK: {
                    #pragma omp parallel
                    collide_lbgk();
                    break;
                }
                case CollisionType::KBC: {
                    #pragma omp parallel
                    collide_kbc();
                    break;
                }
            }
            timer.stop();
            std::cout << "Collision took " << timer.elapsedMilliseconds() << " ms." << std::endl;



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
            os << "domain: " << sim.l.nx << " x " << sim.l.ny << "\n";
            os << "Re:     " << sim.Re << "\n";
            os << "Vmax:   " << sim.Vmax << "\n";
            os << "visc:   " << sim.visc << "\n";
            os << "beta:   " << sim.beta << "\n";
            return os;
        }
    };

} // lb

#endif // LB_SIMULATION_HPP_INCLUDED
