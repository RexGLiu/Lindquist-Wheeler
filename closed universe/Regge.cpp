#include <cstdlib>
#include "Regge.h"
using namespace std;

//const double two_pi = 2*3.14159265358979323846264338327;
const double two_pi = 2*acos(-1);
#define err_tol DBL_EPSILON

//*********************************************************************************
//Function macro computes range in tau coordinate for a given rho-coord.
#define tau_range(rho_p) (tau[avg] + tanhb*(rho_p))

//Function macro computes range in psi coordinate for a given rho-coord.
#define psi_range(rho_p) (psi[avg] + tana*(rho_p))

//*********************************************************************************
//d_WE = r_l/rt_l - R*log(sqrt(r_l/R) - sqrt(r_l0/R));
#define update_d0_WE(RParams) (RParams)[d0_WE] = (RParams)[side_r]*(RParams)[rt0_rt] \
    + R*log((RParams)[side_rt0] + (RParams)[side_rt])
    //+ const -- ignorable since we ultimately care only abt differences

#define get_stepsize() (fwd_prefactor*(RadParams[lwr][side_rt0]*RadParams[lwr][side_rt]+RadParams[lwr][side_r0]))

#define inv_stepsize() (inv_prefactor*(2*RadParams[uppr][side_r]-inv_const1 \
        +sqrt(4*RadParams[uppr][side_r]*(RadParams[uppr][side_r]-inv_const2)-inv_const3)))
        //r = upper radius from which block size will be deduced.

//Calculates a provisional value for the new forward block length.
#define provis_fwd_blockln() calc_blockln(get_stepsize())

#define back_blockln() calc_blockln(inv_stepsize())

//***************** Functions to compare tau and radius of two particles **********

//Returns radius of block A minus radius of block B.
//Schw_dr of A and B must equal
double block::compare_rad(const block &A, const block &B)
{
    compare_ck(A.Schw_dr, B.Schw_dr, block_err::const_Schw_dr);
	return long(A.rad-B.rad)*A.Schw_dr;
}

//Returns time of block A minus time of block B
//Schw_dt of A and B must equal
double block::compare_time(const block &A, const block &B)
{
    compare_ck(A.Schw_dt, B.Schw_dt, block_err::const_Schw_dt);
	return (A.t-B.t)*A.Schw_dt;
}

//Used by the above functions to check if valid comparison being made.
//A and B are values being compared, and ctype is type of comparison attempted.
void block::compare_ck(const double val1, const double val2, block_err::err_type type)
{
    if (fabs(val1-val2) > err_tol*val1)
	{
        block_err err(val1, val2, type);
		throw err;
	}
}

//********************* Functions specific to closed universe *************************
//Given a trajectory traj[], positions posn[], and target Schwarzschild time,
//function determines lambda-value L needed for trajectory to reach target time.
double block::t_to_L(const double time, const double posn[], const double traj[]) const
{
	double tmp = 2*(time - t*Schw_dt)-Schw_dt;
	return (tmp*tau_max(posn[0]) - posn[2]*Schw_dt)/(traj[2]*Schw_dt - tmp*traj[0]*tanhb);
}

//************************** block class constructor **********************************

//Initialise a new block given Schwarzschild (r,phi,t) coords in the block.
//S_coord[] stores the Schwarzschild coords, and corresponding block coords will
//be stored in b_coord[]
block::block(double delta_t, const unsigned long int ang_divs, double delta_r, double in_R,
			 const double S_coord[], double b_coord[], const double in_tol) :
N_ang(ang_divs), R(in_R), 
Schw_dr(delta_r), Schw_dphi(two_pi/ang_divs), Schw_dt(delta_t),
Schw_dr_2(Schw_dr/2), Schw_dphi_2(Schw_dphi/2), Schw_dphi2_4(Schw_dphi*Schw_dphi/4), Schw_dt_2(Schw_dt/2),
dtau_f(Schw_dt_2*Schw_dr*R), 
rad_tol(in_tol), 
fwd_prefactor(rad_tol/(R*Schw_dr)), inv_prefactor(rad_tol/(2*Schw_dr*(2*rad_tol+R))),
inv_const1(2*R+rad_tol), inv_const2(R-rad_tol), inv_const3(rad_tol*(4*R-rad_tol))
{
	if (0 >= R || R+Schw_dr > S_coord[0] || Schw_dr <= 0 || Schw_dt<=0)
	{
		cerr << "Error: invalid block parameters.\n";
        cerr << "R: " << R << "\n";
		cerr << "Cell Radius and Min radius allowed: " << S_coord[0] << " " << R+Schw_dr << "\n";
		cerr << "Difference: " << S_coord[0] - (R+Schw_dr) << "\n";
        cerr << "Schw_dr: " << Schw_dr << "\n";
		cerr << "Schw_dt: " << Schw_dt << "\n";
        
		exit(EXIT_FAILURE);
	}
	
	reset(S_coord, b_coord);
}

//Copy constructor
block::block(const block &RHS) : 
R(RHS.R), Schw_dr(RHS.Schw_dr), Schw_dphi(RHS.Schw_dphi), Schw_dt(RHS.Schw_dt), 
Schw_dr_2(RHS.Schw_dr_2), Schw_dphi_2(RHS.Schw_dphi_2), Schw_dphi2_4(RHS.Schw_dphi2_4), Schw_dt_2(RHS.Schw_dt_2),
dtau_f(RHS.dtau_f), N_ang(RHS.N_ang), 
rad_tol(RHS.rad_tol), 
fwd_prefactor(RHS.fwd_prefactor), inv_prefactor(RHS.inv_prefactor),
inv_const1(RHS.inv_const1), inv_const2(RHS.inv_const2), inv_const3(RHS.inv_const3)
{copy(RHS);}

//Assignment operator
const block& block::operator=(const block &RHS)
{
    consts_ck(RHS);
	
	if (this != &RHS)
		copy(RHS);
	
	return *this;
}

//Reset the block parameters using the given (r,phi,t) coordinates.
void block::reset(const double S_coord[], double b_coord[])
{
	if (R+Schw_dr > S_coord[0]) {
        block::radial_err err(R, Schw_dr, S_coord[0]);
		throw err;
	}
	
	//We will loop block parameters from r=R+Schw_dr until we reach block that contains radius r0.
	//Initialise parameters for loop.
	block_ln = 1;
	block_dr = Schw_dr;
	rad = 1;
	b_coord[0] = S_coord[0] - (R+Schw_dr);
	
	//Initial radial parameters of block.
    InitRadParams(lwr, Schw_dr);
    InitRadParams(uppr, 2*Schw_dr);
    
	//Loop and update radial params until we reach correct block.
	while (b_coord[0] > block_dr) {
		rad += block_ln;
		//b_coord[0] -= block_dr;
        b_coord[0] = S_coord[0] - (rad*Schw_dr + R);
		
		//First update all lower params
		copy_RadParams(RadParams[lwr], RadParams[uppr]);
			//Note: this function call involves copying d0_WE which hasn't been defined yet.
		
        fwd_blockln();
    }
    
	//Calculate remaining params.
	update_d0_WE(RadParams[lwr]);
	update_d0_WE(RadParams[uppr]);
	d_WE = RadParams[uppr][d0_WE] - RadParams[lwr][d0_WE];	//Note: this d is twice the d defined in Williams-Ellis
	
    //First normalise S_coord[1] to be in range [0, 2pi).
    double tmp_ang = fmod(S_coord[1], two_pi);
    if (tmp_ang < 0) tmp_ang += two_pi;
    
    //Compute ang.
	ang = tmp_ang/Schw_dphi;
    b_coord[1] = fmod(tmp_ang, Schw_dphi);
    if ((Schw_dphi - b_coord[1]) <= err_tol*b_coord[1])
    {
        ang++;
        if (ang == N_ang) ang = 0;
        
        b_coord[1] = 0;
    }
	
    //Compute t based on S_coord[2].  If S_coord[2]==0, then t=0 and tau = -tau_max.
	t = floor(S_coord[2]/Schw_dt);
    double t_remainder = S_coord[2] - t*Schw_dt;
	
	UpdateNonRadParams();
	
	//Convert S_coord[] into block coords
	b_coord[0] = b_coord[0]*2*rho/block_dr - rho;
    b_coord[1] = (S_coord[1]==0 ? 0 : b_coord[1]*S_coord[0] - psi_max(b_coord[0]));
    double max = tau_max(b_coord[0]);
	b_coord[2] = 2*max*t_remainder/Schw_dt - max;
}

//First sets RadParams[s][side_r0] to r0.  Then updates RadParams[s][side_r] accordingly.
//Finally, if update_all is set, function updates remaining RadParams based on RadParams[s][side_r0].
inline void block::InitRadParams(const side s, const double r0, bool update_all)
{
    RadParams[s][side_r0] = r0;
    RadParams[s][side_r] = RadParams[s][side_r0]+R;
    
    if (update_all)
        UpdateRemainingRadParams(s);
}

//Function updates remaining RadParams based on RadParams[s][side_r0].
inline void block::UpdateRemainingRadParams(const side s)
{
    RadParams[s][side_rt0] = sqrt(RadParams[s][side_r0]);
    RadParams[s][side_rt] = sqrt(RadParams[s][side_r]);
    RadParams[s][rt0_rt] = RadParams[s][side_rt0]/RadParams[s][side_rt];
}

//Light propagation without a stopping condition. Light is propagated to the closest point where 
//any of its (rho, psi, tau) coordinates match any of the (rho, psi, tau) coordinates in 'stop'.
//If such a point is outside the block, then light is propagated to the closest face instead.
//The function returns a type "face" indicating the face exited or the stopping condition satisfied.
block::Face block::propagate(double posn[], double traj[], double &L) const 
{
	Face face = top;
	L = (tau_range(posn[0]) - posn[2]) / (traj[2] - traj[0]*tanhb);
    
	register double L_temp, tmp = traj[0]*tana, denom = traj[1] - tmp;
	if (denom > 0)    //necessary condition to intercept back face
	{
		L_temp = (psi_range(posn[0]) - posn[1]) / denom;
        shortest_L(L, L_temp, face, back);
	}
    
    denom = traj[1] + tmp;
    if (denom < 0)    //necessary condition to intercept front face
    {
        L_temp = (psi_range(posn[0]) + posn[1]) / -denom;
        shortest_L(L, L_temp, face, front);
    }
	
	if (traj[0]>0)    //necessary condition to intercept right face
	{
		L_temp = (rho - posn[0]) / traj[0];
        shortest_L(L, L_temp, face, right);
	}
	else if (traj[0]<0)    //necessary condition to intercept left face
	{	
		L_temp = (rho + posn[0]) / -traj[0];
        shortest_L(L, L_temp, face, left);
	}
	
	return face;
}

//Compares L and L_tmp.  If L_tmp is positive but smaller than L, then L gets updated to L_tmp,
//and face gets updated to face_tmp.
inline void block::shortest_L(double &L, const double L_tmp, block::Face &face, const block::Face face_tmp) const
{
    if (L_tmp >= 0 && L > L_tmp)
    {
        L = L_tmp;
        face = face_tmp;
    }
}

//Same light propagation function as above but with stopping condition.
block::Face block::propagate(double posn[], double traj[], double &L, const double stop[]) const 
{
	Face face = propagate(posn, traj, L);
	
	double L_temp = (stop[0] - posn[0])/traj[0];
    shortest_L(L, L_temp, face, rho_stop);

	L_temp = (stop[1] - posn[1])/traj[1];
    shortest_L(L, L_temp, face, psi_stop);
	
	L_temp = (stop[2] - posn[2])/traj[2];
    shortest_L(L, L_temp, face, tau_stop);
	
	return face;
}

//Propagate into next block and update block parameters, posn, and traj
void block::update(double posn[], double traj[], const double L, const Face face)
{
	switch (face)
	{
		case right:
			//update posn[] and block parameters
			rad+=block_ln;
            
			//First update all lower params
			copy_RadParams(RadParams[lwr], RadParams[uppr]);
            
            fwd_blockln();
			
            update_d0_WE(RadParams[uppr]);
			d_WE = RadParams[uppr][d0_WE] - RadParams[lwr][d0_WE];
            //Note: this d is twice the d defined in Williams-Ellis
			
			UpdateNonRadParams(right);
            
			posn[0] = -rho;
			posn[1] += L*traj[1];
			posn[2] += L*traj[2];
			
			return;
		case left:
			//First update all upper params
			copy_RadParams(RadParams[uppr], RadParams[lwr]);
			
			//Then determine the new block_ln based on new RadParams[uppr][side_r].
			block_ln = back_blockln();
			
			rad -= block_ln;
			block_dr = block_ln*Schw_dr;
			
			//Update lower params
            InitRadParams(lwr, rad*Schw_dr);
			update_d0_WE(RadParams[lwr]);
			d_WE = RadParams[uppr][d0_WE] - RadParams[lwr][d0_WE];	//Note: this d is twice the d defined in Williams-Ellis
            
			UpdateNonRadParams(left);
            
			posn[0] = rho;
			posn[1] += L*traj[1];
			posn[2] += L*traj[2];
			
			return;
		case top:
			//update posn
            double tmp_posn[3];
            tmp_posn[0] = posn[0];
            tmp_posn[1] = posn[1];
            tmp_posn[2] = posn[2];
            
			posn[0] += L*traj[0];
			posn[1] += L*traj[1];
            posn[2] = -tau_max(posn[0]);	//This computation is slower than using L*traj[0]
                            //but reduces accumulation of floating-point error.
			
			//update traj
            trig_denom = (d_WE - dpsi)*(d_WE + dpsi)/2;
            
            traj0 = traj[0];
			traj[0] += dtau*(dtau*traj0 - rho_x_2*traj[2])/trig_denom;
            traj[2] += dtau*(dtau*traj[2] - rho_x_2*traj0)/trig_denom;
            
			//update corresponding Schwarz coords
			t++;
            
			return;
		case back:
			//update posn
			posn[0] += L*traj[0];
            posn[1] = -psi_max(posn[0]);  	//This computation is slower than using L*traj[1]
                            //but reduces accumulation of floating-point error.
			posn[2] += L*traj[2];
			
			//update traj
            sin2a = rho_x_2*dpsi;
            trig_denom = (d_WE*d_WE + dtau*dtau)/2;
            
			traj0 = traj[0];
            traj[0] += (sin2a*traj[1] - dpsi2*traj0)/trig_denom;
			traj[1] -= (dpsi2*traj[1] + sin2a*traj0)/trig_denom;
			
			//update corresponding Schwarz coords
			ang++;
            if (ang == N_ang) ang = 0;
            
			return;
		case front:
			//update posn
			posn[0] += L*traj[0];
            posn[1] = psi_max(posn[0]);     //This computation is slower than using L*traj[1]
                            //but reduces accumulation of floating-point error.
			posn[2] += L*traj[2];
			
			//update traj
            sin2a = rho_x_2*dpsi;
            trig_denom = (d_WE*d_WE + dtau*dtau)/2;
            
			traj0 = traj[0];
            traj[0] -= (dpsi2*traj0 + sin2a*traj[1])/trig_denom;
			traj[1] += (sin2a*traj0 - dpsi2*traj[1])/trig_denom;
            
			//update corresponding Schwarz coords
            if (ang == 0) ang = N_ang-1;
			else ang--;
			
			return;
		case rho_stop:
		case psi_stop:
		case tau_stop:
			posn[0] += L*traj[0];
			posn[1] += L*traj[1];
			posn[2] += L*traj[2];
			return;
		default:
			cerr << "Error: Face is undefined.  Unable to propagate." << endl;
			throw face;
			return;
	}
}

//Updates position of particle to exit point but doesn't enter next block.
//Afterwards, L set to 0.
void block::update_p(double posn[], const double traj[], double &L, const Face face) const
{
	switch (face)
	{
		case right:
		{
			posn[0] = rho;
			posn[1] += L*traj[1];
			posn[2] += L*traj[2];
			break;
		}
		case left:
		{
			posn[0] = -rho;
			posn[1] += L*traj[1];
			posn[2] += L*traj[2];
			break;
		}
		case top:
			posn[0] += L*traj[0];
			posn[1] += L*traj[1];
			posn[2] = tau_max(posn[0]);	//This computation is slower than using L*traj[0]
                    //but reduces accumulation of floating-point error.
			break;
		case back:
			posn[0] += L*traj[0];
			posn[1] = psi_max(posn[0]);  	//This computation is slower than using L*traj[1]
                        //but reduces accumulation of floating-point error.
			posn[2] += L*traj[2];
			break;
		case front:
			posn[0] += L*traj[0];
			posn[1] = -psi_max(posn[0]);	//This computation is slower than using L*traj[1]
                        //but reduces accumulation of floating-point error.
			posn[2] += L*traj[2];
			break;
		case rho_stop:
		case psi_stop:
		case tau_stop:
        case other:
			posn[0] += L*traj[0];
			posn[1] += L*traj[1];
			posn[2] += L*traj[2];
			break;
		default:
			cerr << "Error: Face is undefined.  Unable to propagate." << endl;
			throw face;
			break;
	}
    
    L = 0;
}

//Given a (rho, psi, tau) coordinate, function returns true if position within block, otherwise returns false.
void block::posn_ck(const double posn[]) const
{
	const double max_tau = tau_range(posn[0]), max_psi = psi_range(posn[0]);
    in_coord_range(posn[0], rho, block_err::invalid_rho);
    in_coord_range(posn[1], max_psi, block_err::invalid_psi);
    in_coord_range(posn[2], max_tau, block_err::invalid_tau);
}

//Given a rho position, function returns the block's range in tau-coordinate at that point.
//For example, if function returns value for TAU, then range is [-TAU, TAU].
//Function throws a block_err exception if rho_p is outside block.
double block::tau_max(const double rho_p) const 
{
    in_coord_range(rho_p, rho, block_err::invalid_rho);
	return tau_range(rho_p);
}

//Given a rho position, function returns the block's range in psi-coordinate at that point.
//For example, if function returns value for PSI, then range is [-PSI, PSI].
//Function throws a block_err exception if rho_p is outside block.
double block::psi_max(const double rho_p) const 
{
    in_coord_range(rho_p, rho, block_err::invalid_rho);
	return psi_range(rho_p);
}

//Checks if 'coord' is within +/-range.  If not, block_err thrown with err_type indicated by 'type'.
//This function is intended to check coordinate ranges, and 'type' should indicate the corresponding coordinate.
inline void block::in_coord_range(double coord, const double range, const block::block_err::err_type type) const
{
    if (range <= 0)
    {
        cerr << "Error: invalid range for in_coord_range()" << endl;
        exit(EXIT_FAILURE);
    }
    
    switch (type) {
        case block_err::invalid_rho:
        case block_err::invalid_psi:
        case block_err::invalid_tau:
            coord = fabs(coord);
            if ((coord - range) > err_tol*range)
            {
                block::block_err err(range, coord, type);
                throw err;
            }
            break;
            
        default:
            cerr << "Error: incorrect use of in_coord_range function.\n";
            cerr << "Type: " << type << endl;
            exit(EXIT_FAILURE);
    }
}

//Function converts rho-coordinate in a block to an approximate Schwarzschild r-coordinate by linearly interpolating
//between the left and right block boundaries.
double block::rho_to_rad(const double rho_p) const 
{return rad*Schw_dr + R + block_dr*(rho_p+rho)/(2*rho);}

//Function converts a Schwarzschild r coordinate to block-rho coordinate by inverse of above function.
double block::rad_to_rho(const double r) const
{return 2*rho*(r - rad*Schw_dr - R)/block_dr - rho;}

//Given a rho-coordinate rho_p, function converts tau-coordinate in block to an approximate Schwarzschild t-coordinate.
double block::tau_to_t(const double tau_p, const double rho_p) const 
{
	double max = tau_max(rho_p);
	return t*Schw_dt + Schw_dt*(tau_p+max)/(2*max);
}

//Given a rho-coordinate rho_p, function converts psi-coordinate in block to an approximate Schwarzschild phi-coordinate.
double block::psi_to_phi(const double psi_p, const double rho_p) const 
{
	double phi_p = ang*Schw_dphi + (psi_p+psi_max(rho_p))/rho_to_rad(rho_p);
	return (phi_p >= two_pi ? phi_p - two_pi : phi_p);
}

//*****************Supporting inline functions**************************************
inline void block::err_ck() const 
{
    if (isnan(rho))
    {
        cerr << "Error: invalid block dimensions.  'rho' is NaN." << endl;
        exit(EXIT_FAILURE);
    }
    
    if (psi[uppr]<psi[lwr])
    {
        cerr << "Error: invalid block dimensions.  Upper 'psi' length is smaller than lower length.\n";
        cerr << "Lengths: " << psi[uppr] << " " << psi[lwr] << "\n";
        cerr << "Differences: " << psi[uppr]-psi[lwr] << endl;
        exit(EXIT_FAILURE);
    }

    if (tau[uppr]<tau[lwr])
    {
        cerr << "Error: invalid block dimensions.  Upper 'tau' length is smaller than lower length.\n";
        cerr << "Lengths: " << tau[uppr] << " " << tau[lwr] << "\n";
        cerr << "Differences: " << tau[uppr]-tau[lwr] << endl;
        exit(EXIT_FAILURE);
    }

    if (rad==0)
    {
        cerr << "Error: 'rad' cannot be zero." << endl;
        exit(EXIT_FAILURE);
    }
    
    if (block_ln < 1)
    {
        cerr << "Error: step size must be positive.\n";
        cerr << "step size: " << block_ln << endl;
        exit(EXIT_FAILURE);
    }
}

inline void block::rad_ck() const
{
	if (rad+block_ln < rad)
	{
		cerr << "Error: 'rad' overflows.\n";
		exit(EXIT_FAILURE);
	}
}

inline void block::copy(const block & RHS)
{
	tau[lwr] = RHS.tau[lwr];
	tau[uppr] = RHS.tau[uppr];
    tau[avg] = RHS.tau[avg];
	psi[lwr] = RHS.psi[lwr];
	psi[uppr] = RHS.psi[uppr];
    psi[avg] = RHS.psi[avg];

	rho = RHS.rho;
    rho2_x_4 = RHS.rho2_x_4;
    rho_x_2 = RHS.rho_x_2;
	tanhb = RHS.tanhb;
	tana = RHS.tana;
	rad = RHS.rad;
	t = RHS.t;
	ang = RHS.ang;
    
    d_WE = RHS.d_WE;
    dpsi = RHS.dpsi;
    dpsi2 = RHS.dpsi2;
    dtau = RHS.dtau;
	
	copy_RadParams(RadParams[lwr], RHS.RadParams[lwr]);
	copy_RadParams(RadParams[uppr], RHS.RadParams[uppr]);
	
	block_ln = RHS.block_ln;
	block_dr = RHS.block_dr;
}

inline void block::copy_RadParams(double LHS[], const double RHS[])
{
	for (int i=0; i<key_tot; i++)
		LHS[i] = RHS[i];
}

inline void block::UpdateNonRadParams(const Face f)
{
	update_tau(f);
	update_psi(f);
	
    rho2_x_4 = (d_WE-dpsi)*(d_WE+dpsi) + dtau*dtau;
    rho_x_2 = sqrt(rho2_x_4);
	rho = rho_x_2/2;
	
	err_ck();
	
	tanhb= dtau/(2*rho);
	tana = dpsi/(2*rho);
}

inline void block::update_tau(const Face f)
{
	switch (f) {
		case right:
			tau[lwr] = tau[uppr];
			tau[uppr] = RadParams[uppr][rt0_rt]*Schw_dt_2;
			break;
		case left:
			tau[uppr] = tau[lwr];
			tau[lwr] = RadParams[lwr][rt0_rt]*Schw_dt_2;
			break;
		default:
			tau[lwr] = RadParams[lwr][rt0_rt]*Schw_dt_2;
			tau[uppr] = RadParams[uppr][rt0_rt]*Schw_dt_2;
			break;
	}
	
    tau[avg] = (tau[lwr]+tau[uppr])/2;

	dtau = dtau_f*block_ln/
			((RadParams[uppr][rt0_rt]+RadParams[lwr][rt0_rt])
			 *RadParams[lwr][side_r]*RadParams[uppr][side_r]);
}

inline void block::update_psi(const Face f)
{
	switch (f) {
		case right:
			psi[lwr] = psi[uppr];
			psi[uppr] = RadParams[uppr][side_r]*Schw_dphi_2;
			break;
		case left:
			psi[uppr] = psi[lwr];
			psi[lwr] = RadParams[lwr][side_r]*Schw_dphi_2;
			break;
		default:
			psi[lwr] = RadParams[lwr][side_r]*Schw_dphi_2;
			psi[uppr] = RadParams[uppr][side_r]*Schw_dphi_2;
			break;
	}
    
    psi[avg] = (psi[lwr]+psi[uppr])/2;
	
	dpsi = Schw_dphi_2*block_dr;
    dpsi2 = block_dr*block_dr*Schw_dphi2_4;
}

inline unsigned long int block::calc_blockln(const double step) const
{
	unsigned long int len = step;
	if (len+1 - step < err_tol) len++;

	//Block length must be at least 1.
	return (len > 1 ? len : 1);
}

//Determines the new forward block length and updates block_ln accordingly.
void block::fwd_blockln()
{
    //Determine the new block_ln
    block_ln = provis_fwd_blockln();
    
    //Check adding block_ln will not overflow rad.
    rad_ck();
    
    //Check that reversing will still give same block length.  If not, then decrement block length.
    //First update upper params needed for check.
    InitRadParams(uppr, (rad+block_ln)*Schw_dr, false);
    
    if (block_ln>1)
    {
        unsigned long int temp_ln;
        
        //If reversing doesn't give same block_ln, then decrease to match,
        //re-calculate upper params, and check again.
        while ((temp_ln = back_blockln()) < block_ln)
        {
            block_ln = temp_ln;
            InitRadParams(uppr, (rad+block_ln)*Schw_dr, false);
        }
    }
    
    block_dr = block_ln*Schw_dr;
    
    //Update remaining upper param
    UpdateRemainingRadParams(uppr);
}

//********************************** Error handling functions ********************************

//Function checks that all constants of RHS and *this match.
//Upon first mismatch found, it notes which const and throws and error.
void block::consts_ck(const block &RHS) const
{
    if (fabs(R-RHS.R) > err_tol*R)
    {
		block_err err(R, RHS.R, block_err::const_R);
		throw err;
    }
    
    if (fabs(Schw_dr-RHS.Schw_dr) > err_tol*Schw_dr)
    {
        block_err err(Schw_dr, RHS.Schw_dr, block_err::const_Schw_dr);
		throw err;
    }
    
    if (fabs(Schw_dphi-RHS.Schw_dphi) > err_tol*Schw_dphi)
    {
        block_err err(Schw_dphi, RHS.Schw_dphi, block_err::const_Schw_dphi);
		throw err;
    }

    if (fabs(Schw_dt-RHS.Schw_dt) > err_tol*Schw_dt)
    {
        block_err err(Schw_dt, RHS.Schw_dt, block_err::const_Schw_dt);
		throw err;
    }

    if (N_ang != RHS.N_ang)
    {
        block_err err(N_ang, RHS.N_ang, block_err::const_N_ang);
		throw err;
    }

    if (rad_tol != RHS.rad_tol)
    {
		block_err err(rad_tol, RHS.rad_tol, block_err::const_rad_tol);
		throw err;
    }
}

//****************************** Methods for error handling classes *********************
double block::block_err::get_diff() const
{
    ck_dbl();
    return dbl1-dbl2;
}

double block::block_err::get_tol() const
{
    ck_dbl();
    return dbl1*err_tol;
}

double block::block_err::get_excess() const
{
    ck_dbl();
    return fabs(dbl1-dbl2) - dbl1*err_tol;
}

double block::block_err::get_dbl1() const
{
    ck_dbl();
    return dbl1;
}

double block::block_err::get_dbl2() const
{
    ck_dbl();
    return dbl2;
}

unsigned long int block::block_err::get_int1() const
{
    ck_int();
    return int1;
}

unsigned long int block::block_err::get_int2() const
{
    ck_int();
    return int2;
}

void block::block_err::ck_dbl() const
{
    switch (type) {
        case const_N_ang:
            cerr << "Error: cannot access double values for an integer error type." << endl;
            exit(EXIT_FAILURE);
    }
}

void block::block_err::ck_int() const
{
    switch (type) {
        case const_R:
        case const_Schw_dr:
        case const_Schw_dphi:
        case const_Schw_dt:
        case const_rad_tol:
        case invalid_rho:
            cerr << "Error: cannot access int values for a double error type." << endl;
            exit(EXIT_FAILURE);
    }
}