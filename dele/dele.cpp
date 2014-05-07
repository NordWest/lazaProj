//#include "stdafx.h"
#include "dele.h"

int planet_num(char *pname)
{
    if(strstr(pname, "Mercury")) return MERCURY_NUM;
    if(strstr(pname, "Venus")) return VENUS_NUM;
    if(strstr(pname, "Earth")) return GEOCENTR_NUM;
    if(strstr(pname, "EMB")) return EARTH_NUM;
    if(strstr(pname, "Mars")) return MARS_NUM;
    if(strstr(pname, "Jupiter")) return JUPITER_NUM;
    if(strstr(pname, "Saturn")) return SATURN_NUM;
    if(strstr(pname, "Uranus")) return URANUS_NUM;
    if(strstr(pname, "Neptune")) return NEPTUNE_NUM;
    if(strstr(pname, "Pluto")) return PLUTO_NUM;
    if(strstr(pname, "Moon")) return MOON_NUM;
    if(strstr(pname, "Sun")) return SUN_NUM;
    if(strstr(pname, "Sol")) return SUN_NUM;
    if(strstr(pname, "Nutations")) return NUTATIONS_NUM;
    if(strstr(pname, "Librations")) return LIBRATIONS_NUM;
    if(strstr(pname, "Geocentr")) return GEOCENTR_NUM;

    return -1;
}

double det_planet_H(int pl_num)
{
	switch(pl_num)
	{
	case MERCURY_NUM:
		return 0.162;
		break;
	case VENUS_NUM:
		return -1.472;
		break;
	case EARTH_NUM:
		return 3.000;
		break;
	case MARS_NUM:
		return -3.501;
		break;
	case JUPITER_NUM:
		return -5.970;
		break;
	case SATURN_NUM:
		return -5.306;
		break;
	case URANUS_NUM:
		return -0.621;
		break;
	case NEPTUNE_NUM:
		return 0.209;
		break;
	case PLUTO_NUM:
		return 7.118;
		break;
	default:
		return 3.0;
		break;
	}
}


///////////////////////////////////////////

dele::dele()
{
    setlocale(LC_NUMERIC, "C");
    this->fileName = new char[FNLEN];
    strcpy(this->fileName, "");
}

dele::dele(const char *fname_bin)
{
    dele();

    fileName = new char[FNLEN];
    strncpy(fileName, fname_bin, FNLEN);

    init(fileName);
}

dele::~dele()
{

}

int dele::init(const char *jpl_name)
{
        strcpy(this->fileName, jpl_name);

        int res = Initialize_Ephemeris(fileName);

        GetParams(H1, H2, R1);

        return res;
}

int dele::detR(double *x, double *y, double *z, double Time, char *planet, int proizv, int centr, int sk)
{
	int nplanet = planet_num(planet);
	if(nplanet<0) return 1;

	return(detR(x, y, z, Time, nplanet, proizv, centr, sk));
}

int dele::detRtt(double *x, double *y, double *z, double Time, int nplanet, int centr, int sk)
{
    //Initialize_Ephemeris(fileName);
    //List4 *cs;
    QList <double*> cs;
    double *xk, *xk1;
	int npl;

	int i;
	for(i=0; i<11; i++)
	{
        xk = new double[7];
		detR(&xk[0], &xk[1], &xk[2], Time, i, 0, centr, sk);
		detR(&xk[3], &xk[4], &xk[5], Time, i, 1, centr, sk);
                //this->g1041->getElem(rc, 8+i);
                xk[6] = H2.data.constValue[8+i];//rc->value;
        cs << xk;
	}

    xk = new double[7];
    xk1 = new double[7];
	detR(&xk[0], &xk[1], &xk[2], Time, nplanet, 0, centr, sk);
	detR(&xk[3], &xk[4], &xk[5], Time, nplanet, 1, centr, sk);

	npl = nplanet;
    if(nplanet==GEOCENTR_NUM) npl = 2;
        //this->g1041->getElem(rc, 8+npl);
        xk[6] = H2.data.constValue[8+npl];//rc->value;
    cs << xk;

    xk = cs[9];
    xk1 = cs[11];

	double rmod1;

	rmod1 = sqrt(xk1[0]*xk1[0] + xk1[1]*xk1[1] + xk1[2]*xk1[2]);

    *x = -(xk[6] + xk1[6])*xk1[0]/pow(rmod1, 3.0);
	*y = -(xk[6] + xk1[6])*xk1[1]/pow(rmod1, 3.0);
	*z = -(xk[6] + xk1[6])*xk1[2]/pow(rmod1, 3.0);

	double A[3];
	double c;
    c = H2.data.constValue[6];//rc->value;

	double vmod1 = sqrt(xk1[3]*xk1[3] + xk1[4]*xk1[4] + xk1[5]*xk1[5]);
	
	for(i=0; i<3; i++) A[i] = xk[6]*((4.0*xk[6]/rmod1 - pow(vmod1, 2.0))*xk1[i+3] + 4.0*(xk1[i]*xk1[i+3])*xk1[i+3])/(c*c*pow(rmod1, 3.0));

	*x+=A[0];
	*y+=A[1];
	*z+=A[2];

	int k;

	double rmodij, rmodj;


	for(i=0; i<3; i++)
	{
        A[i] = 0.0;

		for(k=0; k<9; k++)
		{
			if(k==npl) continue;
            xk = cs[k];
			rmodij = sqrt((xk[0]-xk1[0])*(xk[0]-xk1[0]) + (xk[1]-xk1[1])*(xk[1]-xk1[1]) + (xk[2]-xk1[2])*(xk[2]-xk1[2]));
			rmodj = sqrt(xk[0]*xk[0] + xk[1]*xk[1] + xk[2]*xk[2]);
			A[i] += xk[6]*((xk[i]-xk1[i])/pow(rmodij, 3.0)-(xk[i])/pow(rmodj, 3.0));
		}
	}

	*x+=A[0];
	*y+=A[1];
	*z+=A[2];

	return 0;
}

int dele::headParam(QString name, double &res)
{
    int i;
    QString param;
    char *oneParam = new char[6];

    for(i=0; i<400; i++)
    {
        strncpy(oneParam, &R1.constName[i][0], 6);
        strcpy(&oneParam[5], "\0");
        param = QString(oneParam).simplified();
        if(QString().compare(name, param)==0)
        {
            res = H2.data.constValue[i];
            return 0;
        }
    }
    return 1;
}

double dele::headParam(QString name)
{
    double res;
    headParam(name, res);
    return res;
}

int dele::detState(double *x, double *y, double *z, double *vx, double *vy, double *vz, double Time, int nplanet, int centr, int sk)
{
    double xt, yt, zt;
    double vxt, vyt, vzt;
    double Em;
    int npl = 0;
    stateData State;

    if((nplanet==GEOCENTR_NUM))
    {
        npl = 1;
        nplanet = EARTH_NUM;
    }

    if(nplanet==MOON_NUM)
    {
        npl = 2;
        nplanet = EARTH_NUM;
    }


    Interpolate_State( Time , nplanet , &State );

    *vx = State.Velocity[0];
    *vy = State.Velocity[1];
    *vz = State.Velocity[2];
    *vx = *vx/H1.data.AU*86400.0;
    *vy = *vy/H1.data.AU*86400.0;
    *vz = *vz/H1.data.AU*86400.0;

    *x = State.Position[0];
    *y = State.Position[1];
    *z = State.Position[2];
    *x = *x/H1.data.AU;
    *y = *y/H1.data.AU;
    *z = *z/H1.data.AU;

    if(npl)
    {
        xt = yt = zt = 0.0;
        vxt = vyt = vzt = 0.0;

        Interpolate_State( Time , MOON_NUM , &State );
        vxt = State.Velocity[0];
        vyt = State.Velocity[1];
        vzt = State.Velocity[2];
        vxt = vxt/H1.data.AU*86400.0;
        vyt = vyt/H1.data.AU*86400.0;
        vzt = vzt/H1.data.AU*86400.0;

        xt = State.Position[0];
        yt = State.Position[1];
        zt = State.Position[2];
        xt = xt/H1.data.AU;
        yt = yt/H1.data.AU;
        zt = zt/H1.data.AU;
        Em = H1.data.EMRAT;

        if(npl==1)
        {
            *x = *x - (1.0/(1.0+Em))*xt;
            *y = *y - (1.0/(1.0+Em))*yt;
            *z = *z - (1.0/(1.0+Em))*zt;

            *vx = *vx - (1.0/(1.0+Em))*vxt;
            *vy = *vy - (1.0/(1.0+Em))*vyt;
            *vz = *vz - (1.0/(1.0+Em))*vzt;


        }

        if(npl==2)
        {
            *x = *x + (Em/(1.0+Em))*xt;
            *y = *y + (Em/(1.0+Em))*yt;
            *z = *z + (Em/(1.0+Em))*zt;

            *vx = *vx + (Em/(1.0+Em))*vxt;
            *vy = *vy + (Em/(1.0+Em))*vyt;
            *vz = *vz + (Em/(1.0+Em))*vzt;
        }
    }

    if(centr)
    {
        Interpolate_State( Time , SUN_NUM , &State );
        *x -= State.Position[0]/H1.data.AU;
        *y -= State.Position[1]/H1.data.AU;
        *z -= State.Position[2]/H1.data.AU;

        *vx -= State.Velocity[0]/H1.data.AU*86400.0;
        *vy -= State.Velocity[1]/H1.data.AU*86400.0;
        *vz -= State.Velocity[2]/H1.data.AU*86400.0;
    }

    if(sk)
    {
        xt = *x;
        yt = *y;
        zt = *z;
        *y = cos(EKV)*yt + sin(EKV)*zt;
        *z = -sin(EKV)*yt + cos(EKV)*zt;

        xt = *vx;
        yt = *vy;
        zt = *vz;
        *vy = cos(EKV)*yt + sin(EKV)*zt;
        *vz = -sin(EKV)*yt + cos(EKV)*zt;
    }

    return 0;
}


int dele::detR(double *x, double *y, double *z, double Time, int nplanet, int proizv, int centr, int sk)
{
    double xt, yt, zt;
    double Em;
    int npl = 0;
    stateData State;
    if((nplanet==GEOCENTR_NUM))
    {
            npl = 1;
            nplanet = EARTH_NUM;
    }

    if(nplanet==MOON_NUM)
    {
            npl = 2;
            nplanet = EARTH_NUM;
    }

    Interpolate_State( Time , nplanet , &State );

    if(proizv)
    {
        *x = State.Velocity[0];
        *y = State.Velocity[1];
        *z = State.Velocity[2];

        *x = *x/H1.data.AU*86400.0;
        *y = *y/H1.data.AU*86400.0;
        *z = *z/H1.data.AU*86400.0;
    }
    else
    {
        *x = State.Position[0];
        *y = State.Position[1];
        *z = State.Position[2];

        *x = *x/H1.data.AU;
        *y = *y/H1.data.AU;
        *z = *z/H1.data.AU;
    }

    if(npl)
    {
        xt = yt = zt = 0.0;

        Interpolate_State( Time , MOON_NUM , &State );
        if(proizv)
        {
            xt = State.Velocity[0];
            yt = State.Velocity[1];
            zt = State.Velocity[2];
            xt = xt/H1.data.AU*86400.0;
            yt = yt/H1.data.AU*86400.0;
            zt = zt/H1.data.AU*86400.0;
        }
        else
        {
            xt = State.Position[0];
            yt = State.Position[1];
            zt = State.Position[2];
            xt = xt/H1.data.AU;
            yt = yt/H1.data.AU;
            zt = zt/H1.data.AU;
        }

        Em = H1.data.EMRAT;

        if(npl==1)
        {
            *x = *x - (1.0/(1.0+Em))*xt;
            *y = *y - (1.0/(1.0+Em))*yt;
            *z = *z - (1.0/(1.0+Em))*zt;
        }

        if(npl==2)
        {
            *x = *x + (Em/(1.0+Em))*xt;
            *y = *y + (Em/(1.0+Em))*yt;
            *z = *z + (Em/(1.0+Em))*zt;
        }

    }

    if(centr)
    {
        this->detR(&xt, &yt, &zt, Time, SUN_NUM, proizv, 0, 0);

        *x -= xt;
        *y -= yt;
        *z -= zt;
    }

    if(sk)
    {
        xt = *x;
        yt = *y;
        zt = *z;
        *y = cos(EKV)*yt + sin(EKV)*zt;
        *z = -sin(EKV)*yt + cos(EKV)*zt;
    }

    return 0;
}
