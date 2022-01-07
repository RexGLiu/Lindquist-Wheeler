#ifndef REGGE_H
#define REGGE_H
#include <iostream>
#include <cmath>
#include <cfloat>
using namespace std;
//const double two_pi = 2*3.14159265358979323846264338327;

class block
{
public:
    enum rad_err {};
    
    class radial_err;
    class block_err
    {
    public:
        enum err_type {const_R, const_Schw_dr, const_Schw_dphi, const_Schw_dt, const_N_ang, const_rad_tol,
            invalid_rho, invalid_psi, invalid_tau, dummy};
        // used to indicate which type of comparison attempted
        
        block_err(const double in1, const double in2, err_type in_type)
        : dbl1(in1), dbl2(in2), type(in_type),
          int1(0), int2(0)
        {}
        block_err(const unsigned long int in1, const unsigned long int in2, err_type in_type)
        : int1(in1), int2(in2), type(in_type),
          dbl1(0), dbl2(0)
        {}
        
        unsigned long int get_int1() const;
        unsigned long int get_int2() const;
        
        double get_dbl1() const;
        double get_dbl2() const;
        double get_diff() const;
        double get_tol() const;
        double get_excess() const;
        err_type get_type() const {return type;}

    protected:
        //Constructor used to create dummy block_err.  Primarily used to create derived class objects when the block_err component is irrelevant.
        block_err() : int1(0), int2(0), type(dummy), dbl1(0), dbl2(0) {}
        
    private:
        const err_type type;
        // indicates which type of const mismatched
        
        const double dbl1, dbl2;  //used if const is a double
        const unsigned long int int1, int2;  //used if const is a long int
        
        void ck_dbl() const;
        void ck_int() const;

    };
	
	//Define an enumeration to indicate which face is exited by a particle or which stopping condition satisfied.
	//0 = top; 1 = back; -1 = front; 2 = right; -2 = left;
	//3 = rho stopping condition; 4 = psi stopping condition; 5 = tau stopping condition.
	//"undef" indicated none of the above -- eg because face variable not yet properly initialised.
	enum Face {undef=-3, left, front, top, back, right, rho_stop, psi_stop, tau_stop, other};
	
	block(double delta_t, const unsigned long int ang_divs, double delta_r, double in_R,
		  const double S_coord[], double b_coord[], const double in_tol=0);
	block(const block &RHS);
	const block& operator=(const block &RHS);
	void reset(const double S_coord[], double b_coord[]);
	
	double get_tau0() const {return tau[lwr];}
	double get_tau1() const {return tau[uppr];}
	double get_psi0() const {return psi[lwr];}
	double get_psi1() const {return psi[uppr];}
	double get_rho() const {return rho;}
	double get_R() const {return R;}
	double get_rad0() const {return rad*Schw_dr+R;}
    double get_rad1() const {return (rad+1)*Schw_dr+R;}
	double get_phi() const {return ang*Schw_dphi;}
	double get_t() const {return t*Schw_dt;}
	
	static double compare_rad(const block &A, const block &B);
	static double compare_time(const block &A, const block &B);
    static void compare_ck(const double val1, const double val2, block::block_err::err_type type);
	
	Face propagate(double posn[], double traj[], double &L) const;
	Face propagate(double posn[], double traj[], double &L, const double stop[]) const;
	void update(double posn[], double traj[], const double L, const Face face);
    void update_p(double posn[], const double traj[], double &L, const Face face) const;
	
	void posn_ck(const double posn[]) const;
	double tau_max(const double rho_p) const;
	double psi_max(const double rho_p) const;
	double rho_to_rad(const double rho_p) const;
	double rad_to_rho(const double r) const;
	double tau_to_t(const double tau_p, const double rho_p) const;
	double psi_to_phi(const double psi_p, const double rho_p) const;
    
private:
	//The following enums are used to label array indices for the different *Params[].
	//The last entry '*_tot' is only used to keep track of the total number of elements in the array.
	enum side{lwr, uppr, avg, side_tot};		//lwr=lower radial face; uppr=upper radial face; avg=average of upper and lower
	enum key{side_r0, side_r, side_rt0, side_rt, rt0_rt, d0_WE, key_tot};  //only relevant for RadParams[].
	
	const double Schw_dr, Schw_dphi, Schw_dt, R,
				 Schw_dr_2, Schw_dphi_2, Schw_dphi2_4, Schw_dt_2, 
				 dtau_f,	//constant factor appearing in every dtau calculation
				 //Following constants used in every forward or reverse calculation of step sizes.
				 rad_tol, fwd_prefactor, inv_prefactor, inv_const1, inv_const2, inv_const3;
	const unsigned long int N_ang;           //Number of angular divisions into which 2*Pi is divided.
	double	RadParams[side_tot-1][key_tot],  //Stores parameters associated with the radial faces.
                                             //'side_tot-1' since we won't store averages for RadParams
			tau[side_tot],
			psi[side_tot];
	double	rho, rho_x_2, rho2_x_4,
			tanhb, tana, 
			block_dr,
			//the vars below are updated only when used
			//they were moved here from update function to reduce overhead, as update fn gets looped a lot
			cosh2b, sinh2b, tanh2, traj0, traj1,
			tan2, cos2a, sin2a, trig_denom, trig_numer_diff, trig_numer_diff2, traj_tmp1, traj_tmp2,
			d_WE, dtau, dpsi, dpsi2;
	unsigned long int ang;	//"ang" stores which angular division, from 0 to N_ang-1, is block located.
	long int t;				//"t" stores block's time in units of Schw_dt.
    unsigned long int rad;	//"rad" stores block's radius in units of Schw_dr, with zero corresponding to horizon.
	//If 'rad' overflows, we can change it back to double.  Then introduce unsigned long int N to number blocks,
	//with block 0 being block R+Schw_dr to R+Schw_dr+Schw_dr*block_ln.
	//When we check if two blocks are same, we check N rather than radii.
	
	inline void copy(const block & RHS);
	inline void copy_RadParams(double LHS[], const double RHS[]);

	inline void err_ck() const;
	inline void rad_ck() const;
    inline void in_coord_range(double coord, const double range, const block_err::err_type type) const;
    
    void consts_ck(const block &RHS) const;
	
    inline void InitRadParams(const side s, const double r0, bool update_all=true);
    inline void UpdateRemainingRadParams(const side s);
	inline void UpdateNonRadParams(const Face f = undef);
			//Give a default value to call default behaviour in update_tau and update_psi.
	inline void update_tau(const Face f);
	inline void update_psi(const Face f);
    inline void shortest_L(double &L, const double L_tmp, block::Face &face, const block::Face face_tmp) const;
};

class block::radial_err
{
public:
    radial_err(const double in_R, const double in_Schw_dr, const double in_Schw_rad) :
        R(in_R), Schw_dr(in_Schw_dr), Schw_rad(in_Schw_rad)
    {}
    
    double get_R() const {return R;}
    double get_Schw_dr() const {return Schw_dr;}
    double get_min_rad() const {return R+Schw_dr;}
    double get_Schw_rad() const {return Schw_rad;}
    double get_diff() const {return R+Schw_dr - Schw_rad;}

private:
    const double R, Schw_dr, Schw_rad;
};

#endif