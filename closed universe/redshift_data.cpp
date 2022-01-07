#include "redshift_data.h"
#include <iomanip>

ostream& operator<<(ostream& output, const Title &T)
{
	return output << T.get_title();
}
ostream& operator<<(ostream& output, const redshift_data &data)
{
	return output	<< setprecision(5) << data.ang_i << "\t" << data.r_obs << "\t"
	<< setprecision(15) << data.z_FLRW << "\t" << setprecision(15) << data.z_LW << "\n";
}
ostream& operator<<(ostream& output, const double * const data)
{
	return output	<< setprecision(5) << data[0] << "\t"
	<< setprecision(15) << data[1] << "\t" << setprecision(15) << data[2] << "\n";
}

string Title::operator()(double ang) const
{
    char tmp [35];
    sprintf(tmp, "%15.10f\nr_obs\tz_FLRW\tz_LW\n", ang);
    return prefix + string(tmp);
}
