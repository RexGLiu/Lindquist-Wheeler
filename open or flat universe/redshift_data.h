#ifndef REDSHIFT_DATA_H
#define REDSHIFT_DATA_H

#include <iostream>
#include <string>
using namespace std;

//Used for setting & formatting title for output.
struct Title
{
public:
	Title(const string & pref) : prefix(pref == "" ? "" : prefix+"\n") {}
	
	string get_title() const {return prefix + string("Init angle\tr_obs\tz_FLRW\tz_LW\n");}
	string operator()(const double ang) const;
private:
	const string prefix;
};

//Package in one place all redshift data to be outputted.
struct redshift_data {double ang_i, r_obs, z_FLRW, z_LW;};

ostream& operator<<(ostream& output, const Title &T);
ostream& operator<<(ostream& output, const redshift_data &data);
ostream& operator<<(ostream& output, const double * const data);

#endif /* defined(REDSHIFT_DATA_H) */
