/* *** 
   From: http://wiki.openstreetmap.org/wiki/Mercator#Spherical_Mercator
*** */

#include <math.h>
#define deg2rad(d) (((d)*M_PI)/180)
#define rad2deg(d) (((d)*180)/M_PI)
#define earth_radius 6378137
 
/* The following functions take or return there results in degrees */
 
// double y2lat_d(double y) { return rad2deg(2 * atan(exp(  deg2rad(y) ) ) - M_PI/2); }
// double x2lon_d(double x) { return x; }
// double lat2y_d(double lat) { return rad2deg(log(tan(M_PI/4+ deg2rad(lat)/2))); }
// double lon2x_d(double lon) { return lon; }
 
/* The following functions take or return there results in something close to meters, along the equator */
 
// double y2lat_m(double y) { return rad2deg(2 * atan(exp( (y / earth_radius ) )) - M_PI/2); }
// double x2lon_m(double x) { return rad2deg(x / earth_radius); }
// double lat2y_m(double lat) { return earth_radius * log(tan(M_PI/4+ deg2rad(lat)/2)); }
// double lon2x_m(double lon) { return deg2rad(lon) * earth_radius; }

#define y2lat_d(y) (rad2deg(2 * atan(exp(  deg2rad(y) ) ) - M_PI/2))
#define x2lon_d(x) (x)
#define lat2y_d(lat) (rad2deg(log(tan(M_PI/4+ deg2rad(lat)/2))))
#define lon2x_d(lon) (lon)

#define y2lat_m(y) (rad2deg(2 * atan(exp( (y / earth_radius ) )) - M_PI/2))
#define x2lon_m(x) (rad2deg(x / earth_radius))
#define lat2y_m(lat) (earth_radius * log(tan(M_PI/4+ deg2rad(lat)/2)))
#define lon2x_m(lon) (deg2rad(lon) * earth_radius)

