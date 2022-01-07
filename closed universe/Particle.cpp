//*********************  Particle manipulation functions  ***********************

#include "lattice_universe.h"

Lattice::Particle::Particle(const particle_params &particle, const Cell &cell, const string &name, const particle_type in_type) :
block_p(cell.dt, cell.N_ang, cell.dr, cell.R, particle.posn, posn, cell.rad_tol), type(in_type), particle_name(name), clock(0), D(particle.D)
{set_particle(particle.traj, cell);}

Lattice::Particle::Particle(const Particle &RHS) : block_p(RHS.block_p), type(RHS.type), particle_name(RHS.particle_name), clock(RHS.clock), D(RHS.D)
{copy_particle(RHS);}

Lattice::Particle::Particle(const Particle &RHS, const string &name) : block_p(RHS.block_p), type(RHS.type), particle_name(name), clock(RHS.clock), D(RHS.D)
{copy_particle(RHS);}

Lattice::Particle::Particle(const Particle &RHS, const Cell &cell, const string &name, const particle_type in_type) :
block_p(RHS.block_p), type(in_type), particle_name(name), clock(RHS.clock), D(RHS.D)
{
	posn[0] = RHS.posn[0];
	posn[1] = (type==other ? RHS.posn[1] : 0);
	posn[2] = RHS.posn[2];
    
    set_particle(RHS.traj, cell);
}

void Lattice::Particle::reset(const particle_params &particle, const Cell & cell)
{
    clock = 0;
    block_p.reset(particle.posn, posn);
    set_particle(particle.traj, cell);
}

void Lattice::Particle::reset(const double new_traj[], const Particle &particle, const Cell & cell)
{
    *this = particle;
    set_particle(new_traj, cell);
}

Lattice::Particle& Lattice::Particle::operator=(const Particle& RHS)
{
	if (this != &RHS)
	{
        if (type != RHS.type)
        {
            cerr << "Error: Particle types do not match.  Cannot copy.\n";
            exit(EXIT_FAILURE);
        }
        if (D != RHS.D)
        {
            cerr << "Error: Particle D parameters do not match.  Cannot copy.\n";
            exit(EXIT_FAILURE);
        }
        
		copy_particle(RHS);
		
        //Copy RHS.block_p into block_p.
        //If any constants mismatch, output which and terminate program.
        try {
            block_p = RHS.block_p;
        }
        catch (const block::block_err &err) {
            process_const_err(err);
        }
	}
	
	return *this;
}

//Copies all particle parameters from RHS except for block_p.
//block_p must be copied separately according to whether copy constructor or assignment is being used.
void Lattice::Particle::copy_particle(const Particle& RHS)
{
	posn[0] = RHS.posn[0];
	posn[1] = RHS.posn[1];
	posn[2] = RHS.posn[2];
	
	traj[0] = RHS.traj[0];
	traj[1] = RHS.traj[1];
	traj[2] = RHS.traj[2];
	
	L = RHS.L;
    
    D = RHS.D;
    E = RHS.E;
    E2 = RHS.E2;
    J = RHS.J;
    J2 = RHS.J2;
	
	face = RHS.face;
    
    clock = RHS.clock;
}

void Lattice::Particle::set_particle(const double in_traj[], const Cell &cell)
{
	L = 0;
	face = block::undef;  //this value would have no meaning for block::update function.
	
	//Set velocity of particle according to type.
	switch (type)
	{
		case comoving:
            //In comoving case, traj deduced from D.
            //comoving_vel(*this, cell);
            set_traj(cell);
            
            J=J2=0;
            E2 = 1-D;
            E = sqrt(E2);
            
            if (E2 < 0) cerr << "Warning: " << particle_name << " E is NaN.\n";
			break;
            
		default:
			traj[0] = in_traj[0];
			traj[1] = in_traj[1];
			traj[2] = in_traj[2];
            
            //In default case, E, E2, J, J2 are deduced from some initial traj.
            double rad = get_rad();
            set_E(rad, cell);
            set_J(rad, cell);
            D = 1-E2;
			break;
	}
    
    if (isnan(traj[0])) {
		cerr << "Error: radial velocity undefined.\n";
		exit(EXIT_FAILURE);
	}
}

//This function will also copy in_posn[] which are block coords of a particle in in_block.
inline void Lattice::Particle::set_particle(const double in_posn[], const double in_traj[],
											const block &in_block, const Cell &cell)
{
	block_p = in_block;
	
	posn[0] = in_posn[0];
	posn[1] = (type==other ? in_posn[1] : 0);
	posn[2] = in_posn[2];
	
	set_particle(in_traj, cell);
}

//'sign' is true if particle outgoing, otherwise it is false
void Lattice::Particle::set_traj(const Cell &cell, const bool sign)
{
    double rad = get_rad();
    
    switch (type) {
        case comoving:
        {
            double v1 = (cell.R - D*rad)/(rad - cell.R), v2 = rad*(1 - D)/(rad - cell.R);
            
            try {
                vel_ck(v1);
                vel_ck(v2);
            } catch (cell_err &err) {
                cerr << "Error: radius exceeds R/D for particle " << particle_name << ".\n";
                cerr << "Radii: " << rad << " " << cell.R/D << "\n";
                cerr << "Difference and tolerance: " << rad-cell.R/D << " " << cell.R/D*err_tol << "\n";
                cerr << "v1: " << v1 << "\n";
                cerr << "v2: " << v2 << "\n";
                
                exit(EXIT_FAILURE);
            }
            
            traj[0] = (sign ? 1 : -1) * sqrt(v1);
            traj[2] = sqrt(v2);
            //Assume traj[1] == 0, which should be case for comoving particles.

            break;
        }
        default:
        {
            double rad_ratio = rad/(rad-cell.R), v0 = E2*rad_ratio - J2/(rad*rad);
            
            if (rad_ratio < 0) {
                cerr << "Error: Particle inside horizon.\n";
                cerr << setprecision(10) << rad << " " << cell.R << " " << rad-cell.R << "\n";
                exit(EXIT_FAILURE);
            }
            if (v0 < 0) {
                if (J2*err_tol > -rad*rad*v0)
                    v0 = 0;
                else {
                    cerr << "Error: Radial velocity is NaN.\n";
                    cerr << setprecision(15) << E2 << " " << rad_ratio << " " << E2*rad_ratio << "\n";
                    cerr << J2 << " " << rad*rad << " " << J2/(rad*rad) << "\n";
                    cerr << v0 << " " << v0/(E2*rad_ratio) << " " << err_tol << " " << E2*rad_ratio*err_tol + v0 << " " << (E2*rad_ratio*err_tol > -v0) << "\n";
                    cerr << v0 << " " << v0/(J2/(rad*rad)) << " " << err_tol << " " << J2*err_tol + rad*rad*v0 << " " << (J2*err_tol > -rad*rad*v0) << "\n";
                    exit(EXIT_FAILURE);
                }
            }
            
            if ((sign && traj[0] < 0) || (!sign && traj[0] > 0)) {
                cerr << "Non-comoving traj have opposite signs: " << sign << " " << traj[0] << "\n";
                throw;
            }
            
            traj[2] = E*sqrt(rad_ratio);
            traj[1] = J/rad;
            traj[0] = (traj[0] < 0 ? -1 : 1)*sqrt(v0);
            
            break;
        }
    }
}

inline void Lattice::Particle::vel_ck(const double vel) const
{
    if (vel < 0) {
        cell_err err(vel, 0);
        throw err;
    }
}

void Lattice::Particle::set_E(const double E_new)
{
    E = E_new;
    E2 = E_new*E_new;
    
    if (E < 0) {
        cerr << "Error: E in particle " << particle_name << " is negative.\n";
        cerr << "E = " << E << "\n";
        throw;
    }
}

void Lattice::Particle::set_E(const double rad, const Cell &cell)
{
    E = traj[2]*sqrt((rad-cell.R)/rad);
    E2 = traj[2]*traj[2]*(rad-cell.R)/rad;
    
    ck_E(rad, cell);
}

void Lattice::Particle::set_J(const double rad, const Cell &cell)
{
    J = traj[1]*rad;
    J2 = J*J;
}

inline void Lattice::Particle::ck_E(const double rad, const Cell &cell) const
{
    if (isnan(E)) {
        cerr << "Error: E in particle " << particle_name << " is NaN.\n";
        cerr << "E2 = " << E2 << "\n";
        cerr << "rad = " << rad << "\n";
        cerr << "R = " << cell.R << "\n";
        cerr << "rad - R = " << rad - cell.R << "\n\n";
        throw;
    }
}