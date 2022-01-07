//****************************	 Functions of Cell class  *******************************

#include "lattice_universe.h"


Lattice::Cell::Cell(const cell_params & parameters) :
cell_params(parameters), r_sing(R + dr + (dr >= err_tol*R ? dr : err_tol*R)) {param_ck();}

//Checks if parameters for cell are valid.
inline void Lattice::Cell::param_ck() const
{
	if (R<=0 || dt<=0 || dr<=0 || N_ang<0)
	{
		cerr << "Error: Invalid cell parameters.\n";
		cerr << "R: " << R << "\n";
        cerr << "dt: " << dt << "\n";
        cerr << "dr: " << dr << "\n";
        cerr << "N_ang: " << N_ang << endl;
		exit(EXIT_FAILURE);
	}
}
