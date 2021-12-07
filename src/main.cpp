
#include "simulation.hpp"
#ifdef USE_OPENGL_VISUALIZATION
#include "visualization.hpp"
#endif
#include <omp.h>
#include <filesystem>



int main(int argc, char *argv[])
{
	omp_set_num_threads(6);
	
	lb::simulation* sim = new lb::simulation(256,256,10000,0.05);
	sim->initialize();
	std::cout << *sim << std::endl;

	#ifdef USE_OPENGL_VISUALIZATION
        std::cout << std::filesystem::current_path() << std::endl;

		try {
			lb::visualization::initialize(sim, argc, argv);
			lb::visualization::get_instance().run();
		}
		catch (std::runtime_error e) {
			std::cerr << e.what() << std::endl;
		}
	
	#else
	
		// Here are some hints for getting aquainted with the lattice class
		// ================================================================
		
		// how to print the lattice:
		// -------------------------
		
		std::cout << sim->l << std::endl;
		
		// how to access the lattice:
		// --------------------------
		
		// 1) access via node proxy
		sim->l.get_node(1,0).f(0) = 2;
		
		// 2) access data directly (make sure you know what you're doing)
		sim->l.f[0][sim->l.index(2,0)] = 3;
		
		// 3) using iterators to nodes
		(sim->l.begin() + sim->l.index(0,0))->f(0) = 1;
		
		
		
		std::cout << sim->l << std::endl;
		
		
		
		// use a loop like this to run the simulation
		
		/*for (unsigned int i=0; i<500; ++i)
		{
			sim->step();
		}*/
	
	#endif
	
	return 0;
}
