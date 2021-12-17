/** 
 *  @file
 *  @author Fabian Bösch
 *  @brief velocity set
 */

#ifndef LB_VELOCITY_SET_HPP_INCLUDED
#define LB_VELOCITY_SET_HPP_INCLUDED

#include "global.hpp"
#include <array>
#include <cmath> 

namespace lb {

struct v9;                // forward declaration
const v9& velocity_set(); // forward declaration
	
/**
 *  @brief Lattice parameters for 9 velocity model.
 *  
 *  This class models a the singleton design pattern. That means there 
 *  exists only one single instance throughout the lifetime of the 
 *  program. To instantiate and access this object use the free function
 *  @ref velocity_set.
 * 
 *  This class holds parameters like lattice weights, molecular 
 *  velocities and speed of sound. It also exposes member functions to 
 *  compute the equilibrium populations.
 */
struct v9 // singleton
{
private:
	
	/** @brief Default constructor */
	v9(){};
	/** @brief Function for instantiating the singleton is a friend */
	friend const v9& lb::velocity_set();

public:
	
	v9(const v9&) = delete;
	v9& operator=(const v9&) = delete;

	
	//                                                     0,       1,       2,       3,       4,       5,       6,       7,       8
	const std::array<float_type, 9>         W =   {{ 16.0/36,  4.0/36,  4.0/36,  4.0/36,  4.0/36,  1.0/36,  1.0/36,  1.0/36,  1.0/36}};   ///< Lattice weights
	
	const std::array<std::array<int, 9>, 2> c = {{{{       0,       1,       0,      -1,       0,       1,      -1,      -1,       1}}, 
	                                              {{       0,       0,       1,       0,      -1,       1,       1,      -1,      -1}}}}; ///< Molecular velocities
	
	const float_type cs = 1.0/std::sqrt(3.0);   ///< Speed of sound
	
	const unsigned int size = 9;                ///< Number of velocities

	/** 
	 *  @brief Compute equilibrium.
	 * 
	 *  Compute f_eq from the locally conserved quantities rho, u and v 
	 *  (see also @ref v9::equilibrate).
	 *  @param[in,out] f_eq Pointer to an array of size 9 to store the computed values
	 *  @param[in]     rho  Local density
	 *  @param[in]     u    Local flow velocity in x-direction
	 *  @param[in]     v    Local flow velocity in y-direction
	 */
	inline void f_eq(float_type* f_eq, float_type rho, float_type u, float_type v) const
	{

		float x_root = sqrtf(1 + 3 * u * u);
		float y_root = sqrtf(1 + 3 * v * v);

		float A = rho * (2 - x_root) * (2 - y_root);
		float BX = (2 * u + x_root) / (1 - u);
		float BY = (2 * v + y_root) / (1 - v);

		for(int i = 0; i < velocity_set().size; i++) {
			float w = velocity_set().W[i];

			float cx = velocity_set().c[0][i];
			float cy = velocity_set().c[1][i];
			
			f_eq[i] = w * A * powf(BX, cx) * powf(BY, cy);
		}
	}

	/** 
	 *  @brief Equilibrate a node.
	 * 
	 *  Compute f_eq from the locally conserved quantities rho, u and v
	 *  and set the node's population to that equilibrium ( see also 
	 *  @ref v9::f_eq).
	 *  @tparam        Node A node type
	 *  @param[in,out] n    Reference to a Node object
	 *  @param[in]     rho  Local density
	 *  @param[in]     u    Local flow velocity in x-direction
	 *  @param[in]     v    Local flow velocity in y-direction
	 */
	template <typename Node>
	inline void equilibrate(Node& n, float_type rho, float_type u, float_type v) const
	{
		float_type f[9];
		f_eq(f, rho, u, v);
		n.f(0) = f[0];
		n.f(1) = f[1];
		n.f(2) = f[2];
		n.f(3) = f[3];
		n.f(4) = f[4];
		n.f(5) = f[5];
		n.f(6) = f[6];
		n.f(7) = f[7];
		n.f(8) = f[8];
	}

	/** 
	 *  @brief Equilibrate a node.
	 * 
	 *  Compute f_eq from the locally conserved quantities rho, u and v
	 *  and set the node's population to that equilibrium ( see also 
	 *  @ref v9::f_eq and v9::equilibrate). The locally conserved 
	 *  quantities are taken form the node object itself.
	 *  @tparam        Node A node type
	 *  @param[in,out] n    Reference to a Node object
	 */
	template <typename Node>
	inline void equilibrate(Node& n) const
	{
		return equilibrate(n, n.rho(), n.u(), n.v());
	}
};

/**
 *  @brief Get a reference single instance of the velocity set.
 *  @return 9-velocity set
 */
inline const v9& velocity_set()
{
	static v9 v_set;
	return v_set;
}

} // lb

#endif // LB_VELOCITY_SET_HPP_INCLUDED
