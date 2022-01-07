#include <cfloat>
#include "lattice_universe.h"
#include "loop_comparison.h"
using namespace std;

//*********************************************************************************
#define propagate(particle) (particle).face = (particle).block_p.propagate((particle).posn, \
    (particle).traj, (particle).L);

#define calc_next_tau(p) ((p).face==block::top ? \
    (p).block_p.tau_max((p).posn[0]+(p).L*(p).traj[0]) - (p).posn[2] : (p).L*(p).traj[2])

#define L_intersect(particle1, particle2) \
    ( ((particle1).posn[2]-(particle2).posn[2])/(particle2).traj[2] \
    - ((particle1).posn[0]-(particle2).posn[0])/(particle2).traj[0]) \
	/((particle1).traj[0]/(particle2).traj[0] - (particle1).traj[2]/(particle2).traj[2])


//********************************* Lattice class **********************************
//Constructor for Lattice class.
Lattice::Lattice(const cell_params & in_cell, const particle_params & b, const particle_params & l) :
cell(in_cell), l_init(l, cell, "light"), s_init(l, cell, "source", comoving),
        //set source particle using light block parameters and posn
b_init(init_boundary(b.posn[0]), cell, "boundary", comoving),
light(l_init), source(s_init), boundary(b_init),
tau_surface(l, b.D, cell, "tau_surface", tau),			//must init const members of tau_surface. other params arbitrary.
a0(b.posn[0])
{
	//Calculate initial photon frequency.
	nu0 = photon_freq();
}

//Reset Lattice for next simulation.
void Lattice::reset(const double traj_l[]) 
{
	//Reset boundary
	boundary = b_init;
	
	//Reset light
    light.reset(traj_l, l_init, cell);

	//Reset source
	source = s_init;
	
	//Reset initial photon frequency
	nu0 = photon_freq();
}

//Reset Lattice for next simulation and specify starting freq.
void Lattice::reset(const double traj_l[], const double init_freq)
{
    //Reset initial photon frequency
    nu0 = init_freq;

    //Reset boundary
    boundary = b_init;
    
    //Reset light
    light.reset(traj_l, l_init, cell);
    
    //Reset source
    source = s_init;
    
    //Reset light again w starting trajectory adjusted to match input frequency
    double tmp_freq = photon_freq(), new_traj[3];
    new_traj[0] = light.traj[0]*nu0/tmp_freq;
    new_traj[1] = light.traj[1]*nu0/tmp_freq;
    new_traj[2] = light.traj[2]*nu0/tmp_freq;
    light.reset(new_traj, l_init, cell);
}


Lattice::Particle Lattice::init_boundary(const double rad)
{
    //First ensure target boundary size is valid.
    init_boundary_ck(rad);
    
    //Prepare tau particle for propagation.
	tau_surface.set_particle(light.posn, light.traj, light.block_p, cell);

    //Propagate until reaching block containing target radius.
    reach_init_boundary(rad);
    
    return tau_surface;
}

//Checks whether target radius is greater than initial light radius.
inline void Lattice::init_boundary_ck(const double rad) const
{
    if (rad < l_init.get_rad() )
    {
        cell_err err(rad, l_init.get_rad());
        throw err;
    }
}

//Helper function for init_boundary().  Propagates tau_surface to target radius.
inline void Lattice::reach_init_boundary(double rad)
{
    propagate(tau_surface);
    
    while (tau_surface.block_p.get_rad1() <= rad) {
        update(tau_surface);
        propagate(tau_surface);
    }

    //Convert rad to block-rho coords.
    rad = tau_surface.block_p.rad_to_rho(rad);
    double L_tmp = (rad-tau_surface.posn[0])/tau_surface.traj[0];
    
    while (L_tmp > tau_surface.L) {
        if (tau_surface.face == block::right) {
            cerr << "Crossing right block in reach_init_boundary loop.\n";
            throw;
        }
        
        update(tau_surface);
        propagate(tau_surface);
        L_tmp = (rad-tau_surface.posn[0])/tau_surface.traj[0];
    }

    //Update tau_surface position to target position.
    //Note: posn[1] should still be 0 bc particle confined to tau-rho plane.
    tau_surface.posn[0] = rad;
    tau_surface.posn[2] += L_tmp*tau_surface.traj[2];
}

//Propagates photon to next cell, computes both z_LW and z_FRLW, and stores them in relevant vars.
void Lattice::next_redshift(redshift_data &data)
{next_redshift(data.r_obs, data.z_FLRW, data.z_LW);}

void Lattice::next_redshift(double * const data)
{next_redshift(data[0], data[1], data[2]);}

void Lattice::next_redshift(double &r_obs, double &z_FLRW, double &z_LW)
{
	intercept1(light, boundary);
	cross_boundary(light);
    propagate_until_outward(light);
    intercept1(light, source);
    
    source.set_traj(cell, source.traj[0] >= 0);
    light.set_traj(cell);
    
	//compute z_LW redshift
	double tmp = photon_freq();
	z_LW = (nu0 - tmp)/tmp;
	
	//compute current radius of observation
	r_obs = light.get_rad();
	
	//compute z_FLRW
	tmp = cell_rad();
	z_FLRW = (tmp - a0)/a0;
}

double Lattice::cell_rad()
{
	//Prepare tau particle for propagation.
	tau_surface.set_particle(light.posn, light.traj, light.block_p, cell);
	
	intercept1(tau_surface, boundary);
	
	return boundary.get_rad();
}

//Propagates a particle across boundary.
void Lattice::cross_boundary(Particle & particle)
{
    //First check to make sure particle at boundary
    at_boundary_ck(particle);
    particle.set_traj(cell);
    boundary.set_traj(cell);
    
    //Prepare repeated variables used in calculation of both new E and new traj[] components.
    double rad = boundary.get_rad(),
    factor1 = cell.R+(2*boundary.get_E2()-1)*rad,
    factor2 = boundary.get_E2()*calc_factor2(rad, boundary),
    factor3 = calc_factor3(rad, particle),
    factor4 = calc_factor4(rad),
    E_new = (particle.get_E()*factor1 - 2*sqrt(factor2*factor3)/rad)/factor4;  //calculation of new E

    //Calculation of new traj[] components.  Assumes that both particle.traj[0] and boundary.traj[0] are positive.
    particle.traj[0] = (2*particle.get_E()*sqrt(factor2/factor4)*rad - factor1*sqrt(factor3/factor4)/rad)/factor4;
    particle.traj[2] = (particle.get_E()*factor1*sqrt(rad/factor4) - 2*sqrt(factor2*factor3/(rad*factor4)))/factor4;
    
    particle.set_E(E_new);
}

//Function checks whether particle is at boundary or not.  If not, it throws an exception.
inline void Lattice::at_boundary_ck(const Particle & particle) const
{
    try {
        if (block::compare_rad(boundary.block_p, particle.block_p) != 0)
        {
            boundary_err err(boundary.block_p.get_rad0(), particle.block_p.get_rad0(), Lattice::boundary_err::block_rad);
            throw err;
        }
        if (block::compare_time(boundary.block_p, particle.block_p) != 0)
        {
            boundary_err err(boundary.block_p.get_t(), particle.block_p.get_t(), Lattice::boundary_err::block_time);
            throw err;
        }
        if (boundary.posn[0] != particle.posn[0])
        {
            boundary_err err(boundary.posn[0], particle.posn[0], Lattice::boundary_err::rho_coord);
            throw err;
        }
        if (boundary.posn[2] != particle.posn[2])
        {
            boundary_err err(boundary.posn[2], particle.posn[2], Lattice::boundary_err::time_coord);
            throw err;
        }
    } catch (const block::block_err &err) {
        switch (err.get_type()) {
            case block::block_err::const_Schw_dr:
                cerr << "Block error: Particle and boundary blocks have different Schw_dr.\n";
                cerr << "Boundary & particle Schw_dr: " << err.get_dbl1() << " " << err.get_dbl2() << "\n";
                cerr << "Difference: " << err.get_diff() << "\n";
                cerr << "Tolerance: " << err.get_tol() << endl;
                exit(EXIT_FAILURE);
            case block::block_err::const_Schw_dt:
                cerr << "Block error: Particle and boundary blocks have different Schw_dt.\n";
                cerr << "Boundary & particle Schw_dt: " << err.get_dbl1() << " " << err.get_dbl2() << "\n";
                cerr << "Difference: " << err.get_diff() << "\n";
                cerr << "Tolerance: " << err.get_tol() << endl;
                exit(EXIT_FAILURE);
            default:
                throw err;
        }
    }
}

//Function assumes light and source are at same spacetime point.  Otherwise result is meaningless.
double Lattice::photon_freq() const
{
    double rad = source.get_rad();
    return (rad*sqrt(light.get_E2()*source.get_E2()) - sqrt(calc_factor2(rad, source)*calc_factor3(rad, light))/rad)/calc_factor4(rad);
}

void Lattice::update(Particle & particle) const
{
	if (particle.face==block::right) 
		rad_ck(particle.block_p.get_rad0() + cell.dr);
	else if (particle.face==block::left) 
		rad_ck(particle.block_p.get_rad0());
	
	particle.block_p.update(particle.posn, particle.traj, particle.L, particle.face);

    //If particle crosses any of radially faces, a correction to traj[] is applied.
    if (particle.face==block::right)
        particle.set_traj(cell);
    else if (particle.face==block::left)
        particle.set_traj(cell, false);
}

//Propagate two particles until they intercept each other.  That is, Schwarzschild r and t (but not angle) coincide.
//This function is designed for two particles travelling in the same direction.
//Hence not suitable if any particles could possibly change directions at some point.
//Function currently assumes outer particle will always travel purely radially outwards.  Otherwise we need to pass
//a L_tmp2 argument into intercept_particles() and we need to call next_ang_block with 'outer' and 'inner' swapped.
void Lattice::intercept1(Particle & particle1, Particle & particle2) 
{
	double L_tmp;
	Particle *inner=&particle1, *outer=&particle2;
	
	//Synchronise the Schwarzschild t coordinate of the two particles' blocks.
    sync_t(particle1, particle2);

    //Next sync the tau-coords.
    if (particle1.posn[2] < particle2.posn[2])
        sync_tau(particle1, particle2);
    else sync_tau(particle2, particle1);
	
	//Make 'inner' point to 'inner particle' and 'outer' to 'outer particle'.  Hence, 'inner' must catch up to 'outer'.
    set_outer_inner(inner, outer);
    
    //No particle can ever catch up to light.
    if (outer->particle_name=="light")
    {
        intercept_err err(intercept_err::intercept1, intercept_err::intercept_light, inner->particle_name, "light");
        throw err;
    }

    //Functors to be used to check whether or not the two particles have propagated to the same block.
    compare_r rad_condition(*outer);
    compare_block_t t_condition(*inner);

	//Loop until particles meet.
	while (true) {
        //Error checking: make sure particles still within their respective blocks.  Also make sure particles haven't passed
        //each other w/o intercepting.
        propagation_ck(inner, outer, intercept_err::intercept1);
		
        //Propagate the two particles until they reach same block (ie same Schwarz t and r coords for their blocks.)
		while (rad_condition(*inner))
		{
            loop_propagate<compare_r>(*inner, rad_condition);
            loop_propagate<compare_block_t>(*outer, t_condition);
		}

        propagate(*inner);
        propagate(*outer);
		
		//If the two particles intercept bf any exit the block, then update particles' coords to be point of interception.
		L_tmp = L_intersect(*inner, *outer);
		if (ck_intercept(*inner, L_tmp))
		{
            intercept_particles(*inner, *outer, L_tmp);
			return;			//particles successfully intercepted, so exit function
		}

        //If particles don't intercept in current block, then propagate particles to next block and loop again.
        
        //If inner particle crosses an angular face before intercepting outer particle, it may intercept outer particle in next block.
        //So propagate inner particle to next block and check interception again.  In next block, separation between particles should be smaller too.
        if (next_ang_block(*inner, *outer))
            return;			//particles successfully intercepted, so exit function
        
        //Otherwise, both particles move into blocks of different Schwarz t or r coords.
        update(*inner);
        update(*outer);
        
        //Synchronise the Schwarzschild t coordinate of the two particles' blocks again.
        sync_t(*inner, *outer);
	}
}

//Propagate two particles until they intercept each other.  That is, Schwarzschild r and t (but not angle) coincide.
//Can be used for any two particles travelling in any direction.
//Function assumes that the two particles will definitely intercept.
//An infinite loop will ensue should this not be the case.
void Lattice::intercept2(Particle & particle1, Particle & particle2) 
{
	double tau_i, tau_o;
	Particle *inner=&particle1, *outer=&particle2;
	
	//Synchronise the Schwarzschild t coordinate of the two particles' blocks.
    sync_t(particle1, particle2);
    
    //Next sync the tau-coords.
    if (particle1.posn[2] < particle2.posn[2])
        sync_tau(particle1, particle2);
    else sync_tau(particle2, particle1);
	
	//Make 'inner' point to "inner particle" and 'outer' to "outer particle".  Hence, 'inner' must catch up to 'outer'.
    set_outer_inner(inner, outer);
    
	bool rho_loop;  //used for outer-loop condition
	//loop through the Regge tau-layers until 'inner' & 'outer' intersect
	while(true) {
        
		//When the two particles enter a new tau layer, tau coordinate of 'outer' will always
		//be more neg than tau coordinate of 'inner' because 'outer' has the greater radius.
		//Propagate 'outer' to make the tau coordinate of the two particles equal.
        sync_tau(*outer, *inner);
		
		//reset the face variables so that we can propagate through next tau layer.
		inner->face = outer->face = block::undef;

		//Increment the two particles until they reach the same block or until they both reach the next tau-layer
		//Each time, we increment by the smallest delta_tau necessary for one of the particles to reach a face.
        rho_loop = true;
		while (rho_loop)
		{
            if (inner->face == block::undef)
                propagate(*inner);
            if (outer->face == block::undef)
                propagate(*outer);
			
			//Error checking to ensure against propagation errors
            propagation_ck(inner, outer, intercept_err::intercept2);
            
            //Used to check whether light and source are not intercepting--need to verify correctness.
            if (outer->traj[2]>inner->traj[2]) {
                intercept_err err(intercept_err::intercept2, intercept_err::intercept_light, inner->particle_name, outer->particle_name);
                throw err;
            }
            
			//Check for intersection of particles: first check that particles are in same block
            if (ck_intercept(*inner, *outer))
                return;     //Particles have now intercepted each other. Therefore leave function.
            
			//If the two particles do not meet, then propagate them to next blocks.
            tau_i = calc_next_tau(*inner);
            tau_o = calc_next_tau(*outer);
            
			if (tau_i > tau_o)      //'outer' crosses a face before 'inner'.
			{
                //Note: if 'outer' exits its Regge block before 'inner' exits its block, then face crossed
                //cannot be the top face bc top face of 'outer' has tau-coords >= range of tau-coords
                //within Regge block of 'inner'.
                //However, due to finite floating point precision, this is not always the case.
                //If the fractional tau-difference between the two particles exceeds our
                //floating-point error tolerance, throw an exception.
                if (outer->face == block::top)
                    propagation_ck(tau_i, tau_o);
                
                rho_loop = next_block(*outer, *inner, tau_o);
			} else      //if (tau_i < tau_o)  //'inner' crosses a face before 'outer'.
                rho_loop = next_block(*inner, *outer, tau_i);
		}
	}
}

inline void Lattice::propagate_until_outward(Particle &p)
{loop_propagate<change_traj_dir>(p, change_traj_dir(change_traj_dir::pos));}

//This function also assumes blocks of both particles have the same time coordinate.
//Function propagates 'lagged' particle until either it catches up to 'advanced' particle or crosses into next tau layer.
//L_tmp must be calculated beforehands for 'lagged' particle and passed into function.  L_tmp will be updated as 'lagged' particle propagates.
void Lattice::sync_tau(Particle & lagged, Particle & advanced) const
{
    double L_tmp;
    
	while (lagged.posn[2] < advanced.posn[2])
	{
        //Assuming no deflection of lagged.traj[], L_tmp will denote the 'distance' needed for 'lagged' particle to travel to match
        //tau-coord of 'advanced' particle.
        propagate(lagged);

        //Whenever particle crosses rho-face, vel_b() will make a slight adjustment to traj[2].
        //Therefore, must re-calculate L after crossing each face.
        while ((L_tmp = (advanced.posn[2] - lagged.posn[2])/lagged.traj[2]) >= lagged.L && lagged.face!=block::top)
        {
            update(lagged);
            propagate(lagged);
        }

        //If loop stopped only bc 'lagged' crosses top face before catching up to 'advanced', we must propagate both particles to
		//next tau-layer and repeat process.
		if (L_tmp >= lagged.L && lagged.face==block::top)
		{
            update(lagged);
            next_tau(advanced);
		} else  //otherwise loop stopped bc 'lagged' can catch up in current block (ie L_tmp < lagged.L)
            little_sync_tau(lagged, advanced, L_tmp);
	}
}

inline void Lattice::little_sync_tau(Particle & lagged, const Particle & advanced, const double L) const
{
    lagged.posn[0] += L*lagged.traj[0];
    lagged.posn[2] = advanced.posn[2];
    if (lagged.type == other)
        lagged.posn[1] += L*lagged.traj[1];
}

inline void Lattice::next_tau(Particle &p) const
{
    do {
        propagate(p);
        update(p);
    } while (p.face!=block::top);
}

//Function propagates particles to next block.  p1 must be particle that crosses first into next block.
//tau = p1's final tau coordinate.
inline bool Lattice::next_block(Particle &p1, Particle &p2, const double tau) const
{
    //In general, only inner particle can cross a top face before an outer particle because top face has
    //increasing tau-coord with increasing radius.  However due to floating-point error, this may not necessarily
    //be the case, especially at extreme radii.  We must perform an error check before calling this function.
    
    //If inner particle crosses top face, then the two particles will not meet in current tau-layer.
    //Therefore propagate both particles to next tau-layer.
    if (p1.face == block::top)
    {
        //propagate both particles to the next tau-layer
        update(p1);
        update(p2);
        loop_propagate<compare_face>(p2, compare_face(block::top));
        
        return false;   //break out of rho-loop
    }
    
    //Otherwise p1 does not cross top face first.
    //First propagate p1 to its next block.
    update(p1);
    p1.face = block::undef;
    
    //Then synchronise p2's tau coord with p1.
    double L_tmp = tau/p2.traj[2];
    little_sync_tau(p2, p1, L_tmp);
    propagate(p2);
    
    return true;    //continue w rho-loop
}

//Helper function for intercept1().  Keeps propagating particle p1 to next angular block and checking for interception against p2.
//Function assumes p2 travels purely radially.  Otherwise, we need to pass a L_tmp2 argument into intercept_particles.
inline bool Lattice::next_ang_block(Particle &p1, Particle &p2) const
{
    double L_tmp;
    
    while (p1.face==block::back || p1.face==block::front)
    {
        update(p1);
        propagate(p1);
		
		//Checks for interception, and carries it out if possible.
		L_tmp = L_intersect(p1, p2);
		if (ck_intercept(p1, L_tmp))
		{
            intercept_particles(p1, p2, L_tmp);
			return true;			//indicate particles have intercepted
		}
    }
    
    return false;  //indicate particles did not intercept
}

//Function checks whether two particles will intercept.  If so, then function propagates the two particles to interception
//then returns true.  Otherwise returns false.
inline bool Lattice::ck_intercept(Particle &inner, Particle &outer) const
{
    if (!block::compare_rad(inner.block_p, outer.block_p))
    {
        double  L_tmp1 = L_intersect(inner, outer), L_tmp2 = L_intersect(outer, inner);
        
        //Given particles are in same block, check if they intersect and update accordingly if so.
        if (ck_intercept(inner, L_tmp1) && ck_intercept(outer, L_tmp2))
        {
            intercept_particles(inner, outer, L_tmp1, L_tmp2);
            return true;
        }
    }
    
    return false;
}

//Given a particle and L_tmp, the "distance" to interception with another particle, checks if interception will take place within the block.
inline bool Lattice::ck_intercept(const Particle &p, const double L_tmp) const
{
    return p.L > L_tmp && L_tmp >= 0;
}

inline void Lattice::intercept_particles(Particle &p1, Particle &p2, const double L_tmp1, const double L_tmp2) const
{
    p2.posn[0] = p1.posn[0] += L_tmp1*p1.traj[0];
    p2.posn[2] = p1.posn[2] += L_tmp1*p1.traj[2];
    
    if (p1.type == other)
        p1.posn[1] += L_tmp1*p1.traj[1];
    if (p2.type == other)
        p2.posn[1] += L_tmp2*p2.traj[1];
}

//***************  Definitions of inline helper functions  ***********************
//Helper function which propagates a particle until "condition" is satisfied.
template <class T>
inline void Lattice::loop_propagate(Particle & particle, const T &condition) const
{
	while (condition(particle)) 
	{
		propagate(particle);
		update(particle);
	}
}

//Helper function to synchronise the Schwarzschild t coordinate of the two particles' blocks.
inline void Lattice::sync_t(Particle &particle1, Particle &particle2) const
{
	if (particle1.block_p.get_t() < particle2.block_p.get_t())
        loop_propagate<compare_block_t>(particle1, compare_block_t(particle2));
	else
        loop_propagate<compare_block_t>(particle2, compare_block_t(particle1));
}

//Proceses any block const error thrown.  Outputs const which mismatch and terminates program.
void Lattice::process_const_err(const block::block_err &err)
{
    switch (err.get_type()) {
        case block::block_err::const_R:
            cerr << "Error: Cannot copy blocks.  Const R mismatch.\n";
            cerr << "R values: " << err.get_dbl1() << " " << err.get_dbl2() << "\n";
            cerr << "Difference: " << err.get_diff() << "\n";
            cerr << "Tolerance: " << err.get_tol() << endl;
            exit(EXIT_FAILURE);
            
        case block::block_err::const_Schw_dr:
            cerr << "Error: Cannot copy blocks.  Const Schw_dr mismatch.\n";
            cerr << "Schw_dr values: " << err.get_dbl1() << " " << err.get_dbl2() << "\n";
            cerr << "Difference: " << err.get_diff() << "\n";
            cerr << "Tolerance: " << err.get_tol() << endl;
            exit(EXIT_FAILURE);

        case block::block_err::const_Schw_dphi:
            cerr << "Error: Cannot copy blocks.  Const Schw_dphi mismatch.\n";
            cerr << "Schw_dphi values: " << err.get_dbl1() << " " << err.get_dbl2() << "\n";
            cerr << "Difference: " << err.get_diff() << "\n";
            cerr << "Tolerance: " << err.get_tol() << endl;
            exit(EXIT_FAILURE);

        case block::block_err::const_Schw_dt:
            cerr << "Error: Cannot copy blocks.  Const Schw_dt mismatch.\n";
            cerr << "Schw_dt values: " << err.get_dbl1() << " " << err.get_dbl2() << "\n";
            cerr << "Difference: " << err.get_diff() << "\n";
            cerr << "Tolerance: " << err.get_tol() << endl;
            exit(EXIT_FAILURE);
            
        case block::block_err::const_N_ang:
            cerr << "Error: Cannot copy blocks.  Const N_ang mismatch.\n";
            cerr << "N_ang values: " << err.get_dbl1() << " " << err.get_dbl2() << "\n";
            exit(EXIT_FAILURE);
            
        case block::block_err::const_rad_tol:
            cerr << "Error: Cannot copy blocks.  Const rad_tol mismatch.\n";
            cerr << "rad_tol values: " << err.get_dbl1() << " " << err.get_dbl2() << "\n";
            cerr << "Difference: " << err.get_diff() << "\n";
            cerr << "Tolerance: " << err.get_tol() << endl;
            exit(EXIT_FAILURE);
    }
}

void Lattice::propagation_ck(const double inner, const double outer) const
{
    if ((inner-outer) > err_tol*inner)
    {
        Lattice::propagation_err err(inner, outer, Lattice::propagation_err::tau_ck);
        
        cerr << "Error: outer particle crosses top boundary before inner particle.\n";
        cerr << "Tau params: " << err.get_val1() << " " << err.get_val2() << "\n";
        cerr << "Difference: " << err.get_diff() << "\n";
        cerr << "Tolerance: " << err.get_tol() << "\n";
        cerr << "Amount above tolerance: " << err.get_excess() << "\n";
        throw err;
    }
}

void Lattice::propagation_ck(const Particle *inner, const Particle *outer, Lattice::intercept_err::func_name fn) const
{
    intercept_err::err_type type;
    try {
        type = intercept_err::inner_particle_out;     //record particle type in case an error needs to be thrown.
        inner->block_p.posn_ck(inner->posn);
        
        type = intercept_err::outer_particle_out;     //record particle type in case an error needs to be thrown.
        outer->block_p.posn_ck(outer->posn);
    } catch (const block::block_err &bl_err) {
        ck_block_err_type(bl_err);
        intercept_err err(bl_err, fn, type, inner->particle_name, outer->particle_name);
        throw err;
    }
    
    if (block::compare_rad(inner->block_p, outer->block_p)>0)
    {
        type = intercept_err::no_intercept;     //record particle type to include in error to be thrown.
        intercept_err err(fn, type, inner->particle_name, outer->particle_name);
        throw err;
    }
}

inline void Lattice::ck_block_err_type(const block::block_err &err)
{
    switch (err.get_type()) {
        case block::block_err::invalid_rho:
        case block::block_err::invalid_psi:
        case block::block_err::invalid_tau:
            return;
        default:
            cerr << "Error: unknown block error during propagation.\n";
            throw err;
    }
}

//Sets p1 to point to the inner particle and p2 to the outer particle.
inline void Lattice::set_outer_inner(Particle* &p1, Particle* &p2) const
{
    double delta_rad = block::compare_rad(p1->block_p, p2->block_p);
	if (delta_rad > 0 || (delta_rad==0 && p1->posn[0] > p2->posn[0]))
    {
        Particle *p_tmp = p1;
		p1 = p2;
		p2 = p_tmp;
	}
}

//Checks if given Schwarz radius is lies outside horizon.
inline void Lattice::rad_ck(const double rad) const
{
	if (rad < cell.r_sing)
	{
        sing_err err;
		throw err;
	}
}

//Converts block's rho coord to Schwarzschild r coord.
inline double Lattice::rho_to_rad(const double rho, const block & in_block) const
{
	double rad = in_block.rho_to_rad(rho);
	rad_ck(rad);
	return rad;
}

//*********************  Copy constructor for particle_params  ***********************
Lattice::particle_params::particle_params(const particle_params & in_params)
{
	posn[0] = in_params.posn[0];
	posn[1] = in_params.posn[1];
	posn[2] = in_params.posn[2];
	
	traj[0] = in_params.traj[0];
	traj[1] = in_params.traj[1];
	traj[2] = in_params.traj[2];
}
