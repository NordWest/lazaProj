#ifndef TYPES_DEFINED
#include <ephem_types.h>
#endif

#ifndef EPHEM_UTIL
#include <ephem_util.h>
#endif
#include <ephem_read.h>

#include <QString>
#include <QList>
#include <stdio.h>
#include <math.h>

#ifndef DELE_H

#define EKV 0.409092804
#define FNLEN 255

#define DE405_CAT 192

#define MERCURY_NUM MERCURY
#define VENUS_NUM VENUS
#define EARTH_NUM EARTH
#define MARS_NUM MARS
#define JUPITER_NUM JUPITER
#define SATURN_NUM SATURN
#define URANUS_NUM URANUS
#define NEPTUNE_NUM NEPTUNE
#define PLUTO_NUM PLUTO
#define MOON_NUM MOON
#define SUN_NUM SUN
#define NUTATIONS_NUM NUTATIONS
#define LIBRATIONS_NUM LIBRATIONS
#define GEOCENTR_NUM 15


#define CENTER_BARY 0
#define CENTER_SUN 1

#define SK_EKVATOR 0
#define SK_ECLIPTIC 1


#define SCALE_KMSEC 1.377534792908385547227E+6
#define SCALE_AUDAY 2.3737097724907275533844871479929e+9

int planet_num(char *pname);
double det_planet_H(int pl_num);

class dele
{
public:

    char *fileName;

	dele();
    dele(const char *fname_bin);
	~dele();

    int init(const char *jpl_bin_name);

	int detR(double *x, double *y, double *z, double Time, int nplanet, int proizv, int centr, int sk);
	int detR(double *x, double *y, double *z, double Time, char *planet, int proizv, int centr, int sk);
	int detRtt(double *x, double *y, double *z, double Time, int nplanet, int centr, int sk);
    int detState(double *x, double *y, double *z, double *vx, double *vy, double *vz, double Time, int nplanet, int centr, int sk);

    int headParam(QString name, double &res);
    double headParam(QString name);

	/*
		proiz: 0 - положения, 1 - скорости
		centr: 0 - barycenter SS, 1 - geliocentr
		sk: 0 - ekvator, 1 - ectiptic
	*/

private:

};


#define DELE_H
#endif
