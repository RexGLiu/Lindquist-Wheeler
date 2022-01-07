#ifndef LATTICE_UNIVERSE_H
#define LATTICE_UNIVERSE_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include "Regge.h"
#include "redshift_data.h"
using namespace std;

#define err_tol DBL_EPSILON  //tolerance in fractional error.

//***************  Declaration of Lattice class  ***********************
class Lattice 
{
public:
	//*************************  Supporting data types  *********************
	
	//Used to store particle parameters.  Can also be used to initialise a Particle struct.
	//posn[] are in Schwarzschild coords, but traj[] is in block coords.
	//In both cases, coords are in (rad, ang, t) form.
	struct particle_params 
	{
		double posn[3], traj[3], D;		//Note: D = R/r_max
		
		particle_params(){}
		particle_params(const particle_params & in_params);
	};
	
	//Used to store parameters specific of the lattice cell.  Can also be used to initialise a Lattice class.
	//We only allow const functions to operate on the class.
	struct cell_params 
	{
		double dr, dt, R, rad_tol;
		long int N_ang;
		
		cell_params() : dr(0), N_ang(0), dt(0), R(0), rad_tol(0) {}
		cell_params(const cell_params & rhs) :	dr(rhs.dr), N_ang(rhs.N_ang), dt(rhs.dt), 
												R(rhs.R), rad_tol(rhs.rad_tol) {}
	};
    
    //**********************  Error-handling data types  *********************
	struct cell_err;
	enum sing_err {};
    class boundary_err;
    class propagation_err;
    
    class intercept_err : protected block::block_err
    {
    public:
        enum func_name {intercept1, intercept2};
        enum err_type {no_intercept, inner_particle_out, outer_particle_out, intercept_light};
        
        intercept_err(const block::block_err &err, const func_name in_name, const err_type in_type, const string &in_inner, const string &in_outer) :
        block::block_err(err), name(in_name), type(in_type), inner(in_inner), outer(in_outer)
        {}
        intercept_err(const func_name in_name, const err_type in_type, const string &in_inner, const string &in_outer) :
        block::block_err(), name(in_name), type(in_type), inner(in_inner), outer(in_outer)
        {}
        
        func_name get_func_name() const {return name;}
        err_type get_type() const {return type;}
        block_err::err_type get_block_err_type() const {return block::block_err::get_type();}
        string get_inner() const {return inner;}
        string get_outer() const {return outer;}
        double get_val1() const {return block::block_err::get_dbl1();}
        double get_val2() const {return block::block_err::get_dbl2();}
        double get_diff() const {return block::block_err::get_diff();}
        double get_tol() const {return block::block_err::get_tol();}
        double get_excess() const {return block::block_err::get_excess();}
        
    private:
        const func_name name;
        const err_type type;
        const string inner, outer;
    };
	
    //*************************  Lattice class constructor and reset func  *********************
	Lattice(const cell_params & in_cell, const particle_params & b, const particle_params & l);
	void reset(const double traj_l[]);
    void reset(const double traj_l[], const double init_freq);
    
    //Interface for private function next_redshift().
	void next_redshift(redshift_data &data);
	void next_redshift(double * const data);		//data[] must have length 3
    
    double photon_freq() const;

private:
    //********************** Private Nested Data Structures ************************************
	//For packaging all particle-related data into one entity.
	class Cell : public cell_params
	{
		friend class Lattice;
		
		public:
			Cell(const cell_params & parameters);
		
			//Const functions for calculating velocity of comoving observers at radius rad.
			double vel_b(const double rho, const block & curr_block, const bool sign=true) const;
			void vel_b(double* const traj, const double rho, const block & curr_block, const bool sign=true) const;
		
		private:
			const double r_sing;	//Supporting vars to help with lattice calculations involving cell parameters.
		
			//Const helper functions supporting for calculating velocities of comoving observers.
			//double vel_b(const double rad) const;
			//void vel_b(const double rad, double * const traj) const;
		
			//Const helper functions for error checking
			//void rad_ck(const double rad) const;
			inline void param_ck() const;
	};
	
	//Define a type to label type of particle.
	//"comoving" indicates particle of comoving type
	//"tau" indicates particle moves radially along surface of constant tau_CF
	enum particle_type {comoving, tau, other};
	
	struct Particle
	{
		double posn[3], traj[3], L; 		//all coords are block coords
		block block_p;
		block::Face face;
		const particle_type type;
        const string particle_name;
		
		Particle(const particle_params &particle, const Cell &cell, const string &name, const particle_type in_type = other);
		Particle(const particle_params &particle, const double in_D, const Cell &cell, const string &name, const particle_type in_type = other);
		Particle(const Particle &RHS);
        Particle(const Particle &RHS, const string &name);
        Particle(const Particle &RHS, const Cell &cell, const string &name, const particle_type in_type = other);
        void reset(const particle_params &particle, const Cell & cell);
        void reset(const double new_traj[], const Particle &particle, const Cell & cell);
		Particle& operator=(const Particle &RHS);

		void set_particle(const double in_posn[], const double in_traj[], const block &in_block, const Cell &cell);
        
        void set_traj(const Cell &cell, const bool sign=true);
        void set_E(const double E_new);
        void set_E(const double rad, const Cell &cell);
        void set_J(const double rad, const Cell &cell);
        
        double get_D() const {return D;}
        double get_E() const {return E;}
        double get_E2() const {return E2;}
        double get_J() const {return J;}
        double get_J2() const {return J2;}
        
        double get_t() const {return block_p.tau_to_t(posn[2], posn[0]);}
        double get_rad() const {return block_p.rho_to_rad(posn[0]);}
        
        private:
            double D, E, J, E2, J2;
                //D = 1-E2 = R/r_max.
                //E and J are particle's "energy" and "angular momentum".  E2 = E*E is the LW E.  J2 = J*J.

            inline void copy_particle(const Particle& RHS);
            inline void set_radial_timelike(double &traj0, double &traj2, const Cell &cell, const bool sign);
            void set_particle(const double in_traj[], const Cell &cell);
            inline void ck_E (const double rad, const Cell &cell) const;
            inline void vel_ck(const double vel) const;
	};
	
	//Function objects to serve as loop conditions for particle propagation
	struct compare_block_t;
	struct compare_r;
	struct compare_face;
    struct change_traj_dir;
    
    //********************** Private Data Members ***************************************************
	const Cell cell;								//Stores all cell parameters
	
	//The following stores init parameters of light and source.  Used for simulation resets.
	const Particle l_init, s_init;

	//Particles used in the simulation.
	Particle light, source, tau_surface;
    
    //Because of initialisation order, these particles declared last.
    //b_init stores init parameters of boundary and used for simulation resets.
    //boundary = boundary particle used in simulation
    const Particle b_init;
    Particle boundary;
	
	//a0 = cell size at time tau_CF when photon is emitted.
	//nu0 = initial photon frequency
	const double a0;
	double nu0;
	
    //********************** Private Functions ******************************************************
    //Determines initial position of boundary particle given init cell radius.
    //Returns result in a particle_params struct.
    Particle init_boundary(const double rad);
    inline void reach_init_boundary(double rad);
    
    //Propagates particle to next cell and determines z_LW, z_FLRW,
    //and r_obs which it stores in the three variables passed to it.
    void next_redshift(double &r_obs, double &z_FLRW, double &z_LW);
    
	//Returns radius of cell based on current time tau_CF of light.
	double cell_rad();
    
	//Particle manipulation functions
	void update(Particle & particle) const;
	void intercept1(Particle & particle1, Particle & particle2);
	void intercept2(Particle & particle1, Particle & particle2);
    void sync_tau(Particle & lagged, Particle & advanced) const;
	void cross_boundary(Particle & particle);
    inline void set_outer_inner(Particle* &p1, Particle* &p2) const;
	
    inline void at_boundary_ck(const Particle & particle) const;
    inline bool ck_intercept(Particle &inner, Particle &outer) const;
    inline bool ck_intercept(const Particle &p, const double L_tmp) const;
    inline void intercept_particles(Particle &p1, Particle &p2, const double L_tmp1, const double L_tmp2=0) const;
    inline void propagate_until_outward(Particle &p);
    
    inline void little_sync_tau(Particle & lagged, const Particle & advanced, const double L) const;
    inline void next_tau(Particle &p) const;
    inline bool next_block(Particle &p1, Particle &p2, const double tau) const;
    inline bool next_ang_block(Particle &p1, Particle &p2) const;
    
    //Inline helper functions that calculator factors used in photon_freq() and in
    //calculating new traj when crossing boundary
    inline double calc_factor2(const double rad, const Lattice::Particle &p) const
    {return cell.R-p.get_D()*rad;}  //get_factor2() calculates r*(R/r - R/r_max)
    
    //get_factor3() calculates r^3(E^2 - J^2/r^2(1-R/r))
    inline double calc_factor3(const double rad, const Lattice::Particle &p) const
    {return (p.get_E2()*rad*rad - p.get_J2())*rad + p.get_J2()*cell.R;}
    
    inline double calc_factor4(const double rad) const
    {return rad - cell.R;}
    
    //Const helper function returning Schwarzschild radius of a Regge block's rho-coord.
    double rho_to_rad(const double rho, const block & in_block) const;

	//Helper function for particle manipulation
	template <class T>
	inline void loop_propagate(Particle & particle, const T &condition) const;
	inline void sync_t(Particle &particle1, Particle &particle2) const;
    
    //Error-handling functions
    inline void init_boundary_ck(const double rad) const;
    void propagation_ck(const double inner, const double outer) const;
    void propagation_ck(const Particle *inner, const Particle *outer, intercept_err::func_name fn) const;
    static void process_const_err(const block::block_err &err);
    inline static void ck_block_err_type(const block::block_err &err);
    inline void rad_ck(const double rad) const;
};


//**************************** Definition of error data types *************************************
struct Lattice::cell_err
{
    const double val, ref;
    cell_err(const double in_val, const double in_ref) : val(in_val), ref(in_ref) {}
};

class Lattice::boundary_err
{
public:
    enum err_type{block_rad, block_time, rho_coord, time_coord};
    
    boundary_err(const double in_val1, const double in_val2, const err_type in_type) :
        val1(in_val1), val2(in_val2), type(in_type)
    {}
    double get_val1() const {return val1;}
    double get_val2() const {return val2;}
    double get_diff() const {return val1-val2;}
    err_type get_type() const {return type;}

private:
    const double val1, val2;
    const err_type type;
};

class Lattice::propagation_err
{
public:
    enum ck_type{tau_ck};
    
    double get_val1() const {return val1;}
    double get_val2() const {return val2;}
    double get_diff() const {return val1-val2;}
    double get_tol() const {return val1*err_tol;};
    double get_excess() const {return fabs(get_diff()) - get_tol();}
    ck_type get_type() const {return type;}
    
    propagation_err(const double in_val1, const double in_val2, const ck_type in_type) :
    val1(in_val1), val2(in_val2), type(in_type)
    {}
    
private:
    const double val1, val2;
    const ck_type type;
};

#endif