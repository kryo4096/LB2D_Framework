/** 
 *  @file
 *  @author Fabian Bösch
 *  @brief simulation
 */

#ifndef LB_SIMULATION_HPP_INCLUDED
#define LB_SIMULATION_HPP_INCLUDED

#include "H_root.hpp"
#include "lattice.hpp"
#include <sstream>

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

    public: // ctor


        /**
         *  @brief Construct from domain size and flow parameters
         *  @param[in] nx    extent in x direction
         *  @param[in] ny    extent in y direction
         *  @param[in] _Re   Reynolds number
         *  @param[in] _Vmax mean flow velocity
         */
        simulation(unsigned int nx, unsigned int ny, float_type _Re, float_type _Vmax)
                : l(nx, ny),
                  shift(velocity_set().size),
                  Re(_Re),
                  Vmax(_Vmax),
                  visc( /*fill in your code here*/ 0.001f),
                  beta( /*fill in your code here*/ 0.9f),
                  time(0),
                  file_output(false), // set to true if you want to write files
                  output_freq(100),
                  output_index(0) {
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
        void initialize() {
            // **************************
            // * fill in your code here *
            // * the lines below are    *
            // * just examples          *
            // **************************

			#pragma omp parallel
            for(int i = 0; i < l.nx; i++) {
	            for(int j = 0; j < l.nx; j++) {

					float u = 0;


		            float f_eq[9];
					velocity_set().f_eq(f_eq, 1.0, u, 0.0);

					for(int p = 0; p < 9; p++) {
						l.get_node(i, j).f(p) = f_eq[p];
					}


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
        void collide() {

            #pragma omp for
            for (int j = 0; j < static_cast<int>(l.ny); ++j) {
                for (int i = 0; i < static_cast<int>(l.nx); ++i) {
                    int idx = l.index(i, j);

                    l.rho[idx] = 0;
                    float_type u = 0;
                    float_type v = 0;

                    for (int k = 0; k < velocity_set().size; k++) {
                        l.rho[idx] += l.f[k][idx];
                        u += l.f[k][idx] * (float) velocity_set().c[0][k];
                        v += l.f[k][idx] * (float) velocity_set().c[1][k];
                    }

					l.u[idx] = u / l.rho[idx];
					l.v[idx] = v / l.rho[idx];

	                float f_eq[9];
	                velocity_set().f_eq(f_eq, l.rho[idx], l.u[idx], l.v[idx]);

					for(int k = 0; k < velocity_set().size; k++) {
						l.f[k][idx] += 2 * 0.99f * (f_eq[k] - l.f[k][idx]);
					}


                }
            }

        }

        /** @brief LB step */
        void step() {
            #pragma omp parallel
            {
                advect();
                #pragma omp barrier
                wall_bc();
                #pragma omp barrier
                collide();
                #pragma omp barrier
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
