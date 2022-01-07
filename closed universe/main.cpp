#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cfloat>
#include <cstdlib>
#include "lattice_universe.h"
#include "redshift_data.h"
using namespace std;
//const double pi = 3.14159265358979323846264338327;
const double pi = acos(-1);

double input_params();
string get_name();
string get_address();
double*** alloc_array();
void dealloc_array(double*** data);
void output_open(ofstream &file, const string & filename, bool append=true);
void datatofile(double *** data, const double a0, ofstream &outdata);
void process_compare_err(const block::block_err &err);
void process_intercept_err(const Lattice::intercept_err &err);
void process_boundary_err(const Lattice::boundary_err &err);
void process_cell_err(const Lattice::cell_err &err);

Lattice& init_lattice(Lattice::particle_params & boundary, Lattice::particle_params & light, Lattice::cell_params & cell);
void simulation(Lattice &universe, Lattice::particle_params &light, double *** data);

const int tot_cells = input_params(), N_ang = input_params();
const double ang_max = pi*input_params(), ang_min = pi*input_params(), ang_inc = (ang_max - ang_min)/N_ang;

int main () {
	//Array to store simulation results
	//static double data[N_ang][tot_cells][3]={0};
	double*** data = alloc_array();

	string	address = get_address(),
			filename = get_name();

	ofstream outdata;
	output_open(outdata, address+filename);

	//Set up lattice universe
	Lattice::particle_params boundary;
	Lattice::particle_params light;
	Lattice::cell_params cell;
	Lattice &universe = init_lattice(boundary, light, cell);
    const double a0 = universe.get_a0();

	//Check angles correctly entered
	if (ang_inc <= 0) {
		cerr << "Error: Maximum angle is smaller than minimum angle.\n";
        exit(EXIT_FAILURE);
	}

	//Output all input parameters to screen.
	cout << "Parameters entered are:\n\
		N_ang = " << N_ang << "\n\
		cell.N_ang = " << cell.N_ang << "\n\
		ang_max = " << ang_max << "\n\
		ang_min = " << ang_min << "\n\
		R = " << cell.R << "\n\
		dr = " << cell.dr << "\n\
		dt = " << cell.dt << "\n\
		rad_tol = " << cell.rad_tol << "\n\
		boundary.D = " << boundary.D << "\n\
		light.posn[0] = " << light.posn[0] << "\n\
		light.posn[1] = " << light.posn[1] << "\n\
		light.posn[2] = " << light.posn[2] << "\n\
		source.D = " << light.D << "\n";   //Note: D for light's source should be stored in light.D.
	
	simulation(universe, light, data);
	datatofile(data, a0, outdata);
	
	outdata.close();
	dealloc_array(data);

	return 0;
}


//Possibly group file-handling together into new class
void output_open(ofstream &file, const string & filename, bool append) {
	if (append)
		file.open(filename.c_str(), ios_base::app);
	else
		file.open(filename.c_str());

	if (!file.is_open()) {
		cerr << "Error: output file unopened." << endl;
        exit(EXIT_FAILURE);
	}
}

//Creates, initialises, and returns a lattice universe to be used for duration of entire program.
//Also initialises boundary, light, and cell parameters.
Lattice& init_lattice(Lattice::particle_params & boundary, Lattice::particle_params & light, Lattice::cell_params & cell) {
	//Set cell parameters.
	cell.R = input_params();
	cell.dr = input_params()*cell.R;
	cell.N_ang = input_params();
	cell.dt = input_params()*cell.dr;
	cell.rad_tol = input_params();      //Note: this line is included for closed universe, rather than automatically computed.

	//Initialise boundary parameters.
	//For boundary, traj[] is superfluous.
	boundary.traj[0] = boundary.traj[1] = boundary.traj[2] = 0;
	
	//Initial position of boundary irrelevant.  It will be determined based on boundary's D and source's init params.
	boundary.posn[0] = 0;
	boundary.posn[1] = 0;
	boundary.posn[2] = 0;
	boundary.D = input_params();

	//Initialise light parameters.
	//Start w light moving completely in tangential direction.
	light.traj[0] = 0;
	light.traj[1] = light.traj[2] = 1;
	
	//Set initial position of light
	light.posn[0] = input_params()*cell.R;
	light.posn[1] = input_params();
	light.posn[2] = input_params();
    light.D = input_params();  //This D gets used by source, not light.

	//Create, initialise, return lattice
	static Lattice universe(cell, boundary, light);
	
	return universe;
}

//Main simulation function.
void simulation(Lattice &universe, Lattice::particle_params &light, double *** data) {
	double ang;
	//Title heading("");
	//redshift_data data;
	
	//Loop through simulation of different trajectories.
	for (int i=0; i<N_ang; i++) {
		//Set new trajectory for light
		ang = ang_max - i*ang_inc;
		light.traj[0] = light.traj[2]*cos(ang);
		light.traj[1] = light.traj[2]*sin(ang);
		
		//Store angle for current trajectory in redshift data package.
		//data.ang_i = ang;
		
		//Reset lattice for next simulation
		universe.reset(light.traj);
		
		//Output heading for next trajectory data
		//outdata << heading;

		for (int j=1; j<=tot_cells; j++) {
			try {
				universe.next_redshift(data[i][j-1]);
				//outdata << data << flush;
			}
			catch (Lattice::sing_err) {
                cerr << "Error: radius crosses horizon." << endl;
				break;  //end current trajectory and move on to next
			}
            catch (const block::block_err &err) {
                //Program attempted invalid comparison of tau or rad
                //between particles in two different blocks
                process_compare_err(err);
                exit(EXIT_FAILURE);
            }
            catch (const Lattice::intercept_err &err) {
                process_intercept_err(err);
                exit(EXIT_FAILURE);
            }
            catch (const Lattice::boundary_err &err) {
                process_boundary_err(err);
                exit(EXIT_FAILURE);
            }
            catch (const Lattice::cell_err &err)
            {
                process_cell_err(err);
                exit(EXIT_FAILURE);
            }
            catch (const Lattice::no_redshift &err) {
                cerr << "Warning: Light fails to return to source.  No redshift reported.  Stopping trajectory.\n";
                break;   //continue onto next cell
            }
		}
	}
}

//Output all simulation results to file
void datatofile(double*** data, const double a0, ofstream &outdata) {
	Title heading("");
	double ang;

    outdata  << "Initial boundary radius: " << setprecision(15) << a0 << "\n";
	
	for (int i=0; i<N_ang; i++) {
		ang = ang_max - i*ang_inc;
		
		//Output heading for next trajectory data
		outdata << heading(ang);

		for (int j=1; j<=tot_cells; j++)
			outdata << setprecision(15) << data[i][j-1];
		
		//Double linebreak to separate data series so gnuplot can process.
		outdata << "\n\n";
	}
}

//Dynamically creates the array to store the simulation results.
double*** alloc_array() {
	double *** array = new double**[N_ang];
	
	for (int i=0; i<N_ang; i++) {
		array[i] = new double*[tot_cells];
		
		for (int j=0; j<tot_cells; j++) {
			array[i][j] = new double[3];
			
			for (int k=0; k<3; k++)
				array[i][j][k] = 0;
		}
	}

	return array;
}

//Dynamically de-allocates memory held by data.
void dealloc_array(double*** data) {
	for (int i=0; i<N_ang; i++)
		for (int j=0; j<tot_cells; j++)
			delete [] data[i][j];
	
	for (int i=0; i<N_ang; i++)
		delete [] data[i];
	
	delete [] data;
}

double input_params()
{
	char ch;
	string tmp;
	double x;
	
	while (cin >> ch) {
		//Character '#' signifies rest of line is a comment.  Therefore discard rest of line.
		if (ch=='#')
			getline(cin, tmp);
		
		//Found a digit or start of a decimal point, ie a possible parameter
		else if (isdigit(ch) || ch=='+' || ch=='-' || ch=='.') {
			cin.unget();
			break;
		}
		
		//Not a digit, and not a comment.  Can't process input.
		else {	
			cerr << "Error: Character \'" << ch << "\' entered is invalid.\n";
            exit(EXIT_FAILURE);
		}
	}
	
	cin >> x;
	
	return x;
}

string get_address()
{
	string address, tmp;
	char ch;
	
	//Keep discarding comments until start of filename found.
	while (cin >> ch) {
		if (ch=='#')
			getline(cin, tmp);
		else {
			cin.unget();
			cin >> address;
			break;
		}
	}

	return address;
}

string get_name()
{
	string name = "", tmp;
	char ch;
	
	//Keep discarding comments until start of filename found.
	while (cin >> ch) {
		if (ch=='#')
			getline(cin, tmp);
		else {
			cin.unget();
			break;
		}
	}

	//Parse out the entire line for processing
	getline(cin, tmp);
	istringstream line(tmp);
	
	//Append words to filename until comment or end of line reached
	line >> name;
	while (line >> ch) {
		//Reached a comment, therefore stop.
		if (ch=='#')
			break;
		else line.unget();

		line >> tmp;
		name += " " + tmp;
	}
	
	return name;
}

//******************** Error handling functions *******************
void process_compare_err(const block::block_err &err)
{
    //Print error message according to type of error
    switch (err.get_type())
    {
        case block::block_err::const_Schw_dr :
            cerr << "Error: compare_rad cannot be used if blocks have different Schw_dr.\n";
            break;
        case block::block_err::const_Schw_dt :
            cerr << "Error: compare_time cannot be used if blocks have different Schw_dt.\n";
            break;
    }
    
    //Output Schw_dt or Schw_dr to screen accordingly.
    cerr << "Compare type: " << err.get_type() << "\n";
    cerr << "Parameters: " << err.get_dbl1() << " " << err.get_dbl2() << "\n";
    cerr << "Difference: " << err.get_diff() << "\n";
    cerr << "Tolerated difference: " << err.get_tol() << endl;
}

void process_intercept_err(const Lattice::intercept_err &err)
{
    cerr << "Function: " << (err.get_func_name()==Lattice::intercept_err::intercept1 ? "intercept1\n" : "intercept2\n" );
    cerr << "Intercept error type: " << err.get_type() << endl;
    if (err.get_type() != Lattice::intercept_err::no_intercept)
    {
        cerr << "Block error type: " << err.get_block_err_type() << "\n";
        cerr << "Inner: " << err.get_inner() << "\n";
        cerr << "Outer: " << err.get_outer() << "\n";
        cerr << "Val1: " << err.get_val1() << "\n";
        cerr << "Val2: " << err.get_val2() << "\n";
        cerr << "Diff: " << err.get_diff() << "\n";
        cerr << "Tol: " << err.get_tol() << "\n";
        cerr << "Excess: " << err.get_excess() << "\n";
        cerr << "err_tol: " << DBL_EPSILON << endl;
    }
}

void process_boundary_err(const Lattice::boundary_err &err)
{
    cerr << "Boundary Error:\n";
    cerr << "Error type: " << err.get_type() << "\n";
    cerr << "Val1: " << err.get_val1() << "\n";
    cerr << "Val2: " << err.get_val2() << "\n";
    cerr << "Diff: " << err.get_diff() << endl;
}

void process_cell_err(const Lattice::cell_err &err)
{
    cerr << "Val: " << err.val << "\n";
    cerr << "Ref: " << err.ref << "\n";
    cerr << "Diff: " << err.get_diff() << "\n";
    cerr << "Tol: " << err.get_tol() << "\n";
    cerr << "Excess: " << fabs(err.get_diff()) - err.get_tol() << endl;
}