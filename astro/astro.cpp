#include "astro.h"



double mas_to_rad(double angle)
{
  return acos(-1)*angle/(180*3600000);
};

double mas_to_grad(double angle)
{
	return angle/MAS_IN_GRAD;
}

double ang_b_dir(double a1, double d1, double a2, double d2)
{
  double r1[3];
        r1[0] = cos(mas_to_rad(a1))*cos(mas_to_rad(d1)); 
		r1[1] = sin(mas_to_rad(a1))*cos(mas_to_rad(d1));
		r1[2] = sin(mas_to_rad(d1));
 double r2[3];
        r2[0] = cos(mas_to_rad(a2))*cos(mas_to_rad(d2)); 
		r2[1] = sin(mas_to_rad(a2))*cos(mas_to_rad(d2));
		r2[2] = sin(mas_to_rad(d2));
    return rad_to_mas(acos(r1[0]*r2[0]+r1[1]*r2[1]+r1[2]*r2[2]));
};

double hms_to_mas(QString ra, QString spl_symb)
{
 ra = ra.trimmed(); 
 QStringList ra_list = ra.split(spl_symb, QString::SkipEmptyParts);
 QString nums = ra_list.at(0);
 double h = nums.toDouble();
 nums = ra_list.at(1);
 double m = nums.toDouble();
 nums = ra_list.at(2);
 double s = nums.toDouble();
 return (h*3600.0+m*60.0+s)*15000.0;
};

double damas_to_mas(QString de, QString spl_symb)
{
 de = de.trimmed();
 double sign;
 if (de[0]=='-') sign = -1; else sign = 1; 
 QStringList de_list = de.split(spl_symb, QString::SkipEmptyParts);
 QString nums = de_list.at(0);
 double deg = nums.toDouble();
 deg = fabs(deg);		
 nums = de_list.at(1);
 double arcmin = nums.toDouble();
 nums = de_list.at(2);
 double arcsec = nums.toDouble();
 return sign*(deg*3600+arcmin*60+arcsec)*1000;
};

QString mas_to_hms(double ra, QString spl_symb, int ndig)
{
 int h = (int)(ra/54000000);
 ra = ra-h*54000000;
 int m = (int)(ra/900000);
 ra = ra-m*900000;
 double s = ra/15000;
 QString str = QString( "%1" ).arg( h,2,10,QLatin1Char( '0' ))+spl_symb;
   str = str + QString( "%1" ).arg( m,2,10,QLatin1Char( '0' ))+spl_symb;
   str = str + QString( "%1" ).arg( s,3+ndig,'f',ndig,QLatin1Char( '0' ));
  return str;
};


QString mas_to_damas(double de, QString spl_symb, int ndig)
{
int sign = 0;
if (de<0) sign = 1;
de = fabs(de);
int d = (int)(de/3600000);
de = de - d*3600000;
int am = (int)(de/60000);
de = de - am*60000;
double as = de/1000;
QString str = QString( "%1" ).arg( d,2,10,QLatin1Char( '0' ))+spl_symb;
  str = str + QString( "%1" ).arg( am,2,10,QLatin1Char( '0' ))+spl_symb;
  str = str + QString( "%1" ).arg( as,3+ndig,'f',ndig,QLatin1Char( '0' ));
  if (sign==1)
    str = "-"+str;
  else
    str = "+"+str;
  return str;
};

DATEOBS timeFromStrNA(QString strT)
{
  DATEOBS date_obs;
  QList<QString> ml;
  ml << "Jan" << "Feb" << "Mar" << "Apr" << "May" << "Jun" << "Jul" << "Aug" << "Sep" << "Oct" << "Nov" << "Dec";
  QString s;
  //get month
  s = strT.mid(15,3);
  for (int i=0;i<12;i++)
  {
    if (s == ml.at(i))
	{
	  date_obs.month = i+1;
	  break;
	}
  }
  //get day
  s = strT.mid(19,2);
  date_obs.day = s.toInt();  
  //get year
  s = strT.mid(35,4);
  date_obs.year = s.toInt();
  //get time
  QString sect;
  s = strT.mid(22,12);
  sect = s.section(':',0,0);
  date_obs.hour = sect.toInt();
  sect = s.section(':',1,1);
  date_obs.min = sect.toInt();
  sect = s.section(':',2,2);
  date_obs.sec = sect.toDouble();
  date_obs.pday = date_obs.day+(date_obs.hour*3600+date_obs.min*60+date_obs.sec)/86400;
  return date_obs;
};

DATEOBS timeFromStrFTN(QString strT)
{
DATEOBS date_obs;
  QString s;
  //get month
  s = strT.mid(5,2);
  date_obs.month = s.toInt();
  //get day
  s = strT.mid(8,2);
  date_obs.day = s.toInt();  
  //get year
  s = strT.mid(0,4);
  date_obs.year = s.toInt();
  //get time
  QString sect;
  s = strT.mid(11,12);
  sect = s.section(':',0,0);
  date_obs.hour = sect.toInt();
  sect = s.section(':',1,1);
  date_obs.min = sect.toInt();
  sect = s.section(':',2,2);
  date_obs.sec = sect.toDouble();
  date_obs.pday = date_obs.day+(date_obs.hour*3600+date_obs.min*60+date_obs.sec)/86400;
  return date_obs;

};

int timeFromStrFTN(DATEOBS *date_obs, QString strT)
{
  QStringList pDateStr;
  pDateStr = strT.split("T");
  if(pDateStr.size()!=2) return 1;
  QStringList dParts;
  dParts = pDateStr.at(0).split("-");
  if(dParts.size()!=3) return 1;
  //s = strT.mid(5,2);
  date_obs->month = dParts.at(1).toInt();
  //get day
  //s = strT.mid(8,2);
  date_obs->day = dParts.at(2).toInt();
  //get year
  //s = strT.mid(0,4);
  date_obs->year = dParts.at(0).toInt();
  //get time
  //QString sect;
  dParts.clear();
  dParts = pDateStr.at(1).split(":");
  if(dParts.size()!=3) return 1;
  //s = strT.mid(11,12);
  //sect = s.section(':',0,0);
  date_obs->hour = dParts.at(0).toInt();
  //sect = s.section(':',1,1);
  date_obs->min = dParts.at(1).toInt();
  //sect = s.section(':',2,2);
  date_obs->sec = dParts.at(2).toDouble();
  date_obs->pday = date_obs->day+(date_obs->hour*3600+date_obs->min*60+date_obs->sec)/86400;
  return 0;

};


DATEOBS timeFromStrSDSS(QString strT1, QString strT2)
{
	DATEOBS date_obs;
	QString s;
	s = strT1.section('/',0,0);
	date_obs.day = s.toInt();
	s = strT1.section('/',1,1);
	date_obs.month = s.toInt();
	s = strT1.section('/',2,2);
	date_obs.year = s.toInt();
	if (date_obs.year<10) date_obs.year = date_obs.year+ 2000;
	else date_obs.year = date_obs.year+ 1900;
	///
	s = strT2.section(':',0,0);
	date_obs.hour = s.toInt();
	s = strT2.section(':',1,1);
	date_obs.min = s.toInt();
	s = strT2.section(':',2,2);
	date_obs.sec = s.toDouble();
	date_obs.pday = date_obs.day+(date_obs.hour*3600+date_obs.min*60+date_obs.sec)/86400;
	return date_obs;
}

double getMJDfromYMD(QString strT)
{
	DATEOBS date_obs;
	QString elS = strT.section(' ',0,0); date_obs.year = elS.toInt();
	elS = strT.section(' ',1,1); date_obs.month = elS.toInt();
	elS = strT.section(' ',2,2); date_obs.pday = elS.toDouble();
	date_obs.day = (int)floor(date_obs.pday);
	double hh =   (date_obs.pday - date_obs.day)*24;
	date_obs.hour = (int)floor(hh);
	double mm = (hh - date_obs.hour)*60;
	date_obs.min = (int)floor(mm);
	date_obs.sec = (mm - date_obs.min)*60;
	double mjd = getJD(date_obs) - 2400000.5;
	return mjd;
}

double jd2mjd(double jd)
{
    return(jd-2400000.5);
}
double mjd2jd(double mjd)
{
    return(mjd+2400000.5);
}



double getJD(DATEOBS date_obs)
{
double d_day = date_obs.day;
double d_hour = date_obs.hour;
double d_min = date_obs.min;
double d_sec = date_obs.sec;
d_day = d_day + (d_hour+(d_min+d_sec/60.0)/60.0)/24.0;
int y;
int m;
if (date_obs.month<3)
{
	 y =date_obs.year-1;
	 m =date_obs.month+12;
}
else
{
	 y =date_obs.year;
	 m =date_obs.month;
}

int A = y/100;
double B = 2-A+floor(A/4);
double jd =floor(365.25*y)+floor(30.6001*(m+1))+d_day+1720994.5+B;
/*int A = y + 4716 -(int)((14-m)/12);
double B = fmod(m-3, 12);
double C = date_obs.day - 1;
double g = (int)((int)((A+184)/100)-38);
double jd = int(1461*A/4) + (int)((153*B+2)/5) + C - 1401 -0.5 - g + (d_hour+(d_min+d_sec/60.0)/60.0)/24.0;*/
return jd;
};

double getMJDfromStrNA(QString strT, double exp_time)
{
  DATEOBS date_obs = timeFromStrNA(strT);

  double mjd = getJD(date_obs) - 2400000.5 - exp_time/(2*86400);
  return mjd;
};

double getMJDfromStrFTN(QString strT, double exp_time)
{
  DATEOBS date_obs = timeFromStrFTN(strT);
  double mjd = getJD(date_obs) - 2400000.5 + exp_time/(2*86400);
  return mjd;
};

int getMJDfromStrFTN(double *mjd, QString strT, double exp_time)
{
  DATEOBS date_obs;
  if(timeFromStrFTN(&date_obs, strT)) return 1;
  *mjd = getJD(date_obs) - 2400000.5 + exp_time/(2*86400);
  return 0;
}

double getMJDfromStrT(QString strT)
{
  DATEOBS date_obs = timeFromStrFTN(strT);
  return(jd2mjd(getJD(date_obs)));
};

int getStrFTNfromMJD(QString *strT, double mjd, double exp_time)
{
    DATEOBS date_obs;
    getDATEOBSfromMJD(&date_obs, mjd - (exp_time/2.0)/86400.0);
    getStrFromDATEOBS(strT, date_obs, "", 2, 3);
    return 0;
}

int getStrTfromMJD(QString *strT, double mjd)
{
    DATEOBS date_obs;
    getDATEOBSfromMJD(&date_obs, mjd);
    getStrFromDATEOBS(strT, date_obs, "", 2, 3);
    return 0;
}

double getMJDfromStrSDSS(QString strT1, QString strT2, double exp_time)
{
  DATEOBS date_obs = timeFromStrSDSS(strT1, strT2);
  double mjd = getJD(date_obs) - 2400000.5 + exp_time/(2*86400);
  return mjd;
};

double getMJDfromStrSDSS(QString strT, double exp_time)
{
  DATEOBS date_obs = timeFromStrFTN(strT);
  double mjd = getJD(date_obs) - 2400000.5 + exp_time/(2*86400);
  return mjd;
};

double getYearFromMJD(double mjd)
{
   return 1858.87885010267 + mjd/365.25;
};

double getMJDFromYear(double year)
{
   return (year-1858.87885010267)*365.25;
};

DATEOBS getDATEOBSfromMJD(double mjd, int rnsec)
{
   DATEOBS date_obs;

    getDATEOBSfromMJD(&date_obs, mjd);

    roundDATEOBS(&date_obs, rnsec);

   return date_obs;
};

int roundDATEOBS(DATEOBS *date_obs, int nsec)
{
    int dayM;
    double val;

    val = 60.0 - pow(10, -nsec);

    if(date_obs->sec*pow(10, nsec)>=floor(val*pow(10, nsec)))
    {
        date_obs->sec = 0.0;
        date_obs->min +=1;
    }
    if(date_obs->min>=60)
    {
        date_obs->min -= 60;
        date_obs->hour++;
    }
    if(date_obs->hour>=24)
    {
        date_obs->hour -= 24;
        date_obs->day++;
    }
    dayM = dinm(date_obs->month, isVes(date_obs->year));
    if(date_obs->day>dayM)
    {
        date_obs->day -= dayM;
        date_obs->month++;
    }
    if(date_obs->month>12)
    {
        date_obs->year++;
        date_obs->month -= 12;
    }

}

void getDATEOBSfromMJD(DATEOBS *date_obs, double mjd, int rnsec)
{
    double parth_of_day = mjd - floor(mjd);
    double JD = 2400000.5+mjd;
    double x = JD+0.5;
    int Z = (int)floor(x);
    double F = x - Z;
    int q = (int)floor((Z - 1867216.25)/36524.25);
    int A = Z + 1 + q - (int)(q/4);
    int B = A + 1524;
    int C = (int)floor((B - 122.1)/365.25);
    int D = (int)floor(365.25*C);
    int E = (int)floor((B - D)/30.6001);
    int Day = B - D - (int)floor(30.6001*E) + F;
    int Month;
    if (E<13.5) Month = E - 1;
    else Month = E - 13;
    int Year;
    if (Month<2.5) Year = C - 4715;
    else Year = C - 4716;

    date_obs->year = Year;
    date_obs->month = Month;
    date_obs->day = Day;
    date_obs->pday = Day+parth_of_day;
    date_obs->hour = (int)floor(parth_of_day*24);
    date_obs->min = (int)floor(parth_of_day*1440) - date_obs->hour*60;
    date_obs->sec = parth_of_day*86400 - date_obs->hour*3600 - date_obs->min*60;

    roundDATEOBS(date_obs, rnsec);
}

QString getStrFromDATEOBS(DATEOBS date_obs, QString spl_symb, int format, int ndig)
{
  QString str;
  switch(format)
  {
    case 0:
	{
	  str = QString( "%1" ).arg( date_obs.year,4,10,QLatin1Char( '0' ))+spl_symb;
	  str = str + QString( "%1" ).arg( date_obs.month,2,10,QLatin1Char( '0' ))+spl_symb;
	  str = str + QString( "%1" ).arg( date_obs.pday,3+ndig,'f',ndig,QLatin1Char( '0' ));
	  break;	  
    }
	case 1:
	{
	  str = QString( "%1" ).arg( date_obs.year,4,10,QLatin1Char( '0' ))+spl_symb;
	  str = str + QString( "%1" ).arg( date_obs.month,2,10,QLatin1Char( '0' ))+spl_symb;
	  str = str + QString( "%1" ).arg( date_obs.day,2,10,QLatin1Char( '0' ))+" ";
	  str = str + QString( "%1" ).arg( date_obs.hour,2,10,QLatin1Char( '0' ))+spl_symb;
	  str = str + QString( "%1" ).arg( date_obs.min,2,10,QLatin1Char( '0' ))+spl_symb;
	  str = str + QString( "%1" ).arg( date_obs.sec,3+ndig,'f',ndig,QLatin1Char( '0' ));
	  break;
	}
	case 2:
	{
	  str = QString( "%1" ).arg( date_obs.year,4,10,QLatin1Char( '0' ))+"-";
	  str = str + QString( "%1" ).arg( date_obs.month,2,10,QLatin1Char( '0' ))+"-";
	  str = str + QString( "%1" ).arg( date_obs.day,2,10,QLatin1Char( '0' ))+"T";
	  str = str + QString( "%1" ).arg( date_obs.hour,2,10,QLatin1Char( '0' ))+":";
	  str = str + QString( "%1" ).arg( date_obs.min,2,10,QLatin1Char( '0' ))+":";
	  str = str + QString( "%1" ).arg( date_obs.sec,3+ndig,'f',ndig,QLatin1Char( '0' ));
	  break;
	}
  }
  return str;
};

void getStrFromDATEOBS(QString *rstr, DATEOBS date_obs, QString spl_symb, int format, int ndig)
{
  QString str;
  switch(format)
  {
    case 0:
        {
          str = QString( "%1" ).arg( date_obs.year,4,10,QLatin1Char( '0' ))+spl_symb;
          str = str + QString( "%1" ).arg( date_obs.month,2,10,QLatin1Char( '0' ))+spl_symb;
          str = str + QString( "%1" ).arg( date_obs.pday,3+ndig,'f',ndig,QLatin1Char( '0' ));
          break;
    }
        case 1:
        {
          str = QString( "%1" ).arg( date_obs.year,4,10,QLatin1Char( '0' ))+spl_symb;
          str = str + QString( "%1" ).arg( date_obs.month,2,10,QLatin1Char( '0' ))+spl_symb;
          str = str + QString( "%1" ).arg( date_obs.day,2,10,QLatin1Char( '0' ))+" ";
          str = str + QString( "%1" ).arg( date_obs.hour,2,10,QLatin1Char( '0' ))+spl_symb;
          str = str + QString( "%1" ).arg( date_obs.min,2,10,QLatin1Char( '0' ))+spl_symb;
          str = str + QString( "%1" ).arg( date_obs.sec,3+ndig,'f',ndig,QLatin1Char( '0' ));
          break;
        }
        case 2:
        {
          str = QString( "%1" ).arg( date_obs.year,4,10,QLatin1Char( '0' ))+"-";
          str = str + QString( "%1" ).arg( date_obs.month,2,10,QLatin1Char( '0' ))+"-";
          str = str + QString( "%1" ).arg( date_obs.day,2,10,QLatin1Char( '0' ))+"T";
          str = str + QString( "%1" ).arg( date_obs.hour,2,10,QLatin1Char( '0' ))+":";
          str = str + QString( "%1" ).arg( date_obs.min,2,10,QLatin1Char( '0' ))+":";
          str = str + QString( "%1" ).arg( date_obs.sec,3+ndig,'f',ndig,QLatin1Char( '0' ));
          break;
        }
  }
  rstr->clear();
  rstr->append(str);
};

double getAngleDeg(double cos_x, double sin_x)
{
  double ang = 180*atan2(cos_x, sin_x)/PI;
  if (ang<0)ang = 360 + ang;
  return ang;
};


double* getTangToRaDe(double ksi, double eta, double ra_c, double de_c)
{
    double rd_vector[2];
	ra_c = mas_to_rad(ra_c);
	de_c = mas_to_rad(de_c);
	ksi = mas_to_rad(ksi*1000);
        eta = mas_to_rad(eta*1000);
	double secD = 1/cos(de_c);
	double commP = 1-eta*tan(de_c);
	double ra = atan((ksi*secD)/(commP))+ra_c;
	double de = atan((eta+tan(de_c))*cos(ra - ra_c)/commP);
	rd_vector[0] = rad_to_mas(ra);
	rd_vector[1] = rad_to_mas(de);
	return rd_vector;
};

void getTangToRaDe1(double *rd0, double *rd1, double ksi, double eta, double ra_c, double de_c)
{
	ra_c = mas_to_rad(ra_c);
	de_c = mas_to_rad(de_c);
	ksi = mas_to_rad(ksi*1000);
        eta = mas_to_rad(eta*1000);
	double secD = 1/cos(de_c);
	double commP = 1-eta*tan(de_c);
	double ra = atan((ksi*secD)/(commP))+ra_c;
	double de = atan((eta+tan(de_c))*cos(ra - ra_c)/commP);
	*rd0 = rad_to_mas(ra);
	*rd1 = rad_to_mas(de);
};

void getTangToRaDe(double ksi, double eta, double ra_c, double de_c, double *const ep)
{
	ra_c = mas_to_rad(ra_c);
	de_c = mas_to_rad(de_c);
	ksi = mas_to_rad(ksi*1000);
    eta = mas_to_rad(eta*1000);
	double secD = 1/cos(de_c);
	double commP = 1-eta*tan(de_c);
	double ra = atan((ksi*secD)/(commP))+ra_c;
	double de = atan((eta+tan(de_c))*cos(ra - ra_c)/commP);
	ep[0] = rad_to_mas(ra);
	ep[1] = rad_to_mas(de);
};

double* getRaDeToTang(double ra, double de, double ra_c, double de_c)
{
double rd_vector[2];
ra_c = mas_to_rad(ra_c);
de_c = mas_to_rad(de_c);
ra = mas_to_rad(ra);
de = mas_to_rad(de);
double ksi = cos(de)*sin(ra-ra_c)/(sin(de)*sin(de_c)+cos(de)*cos(de_c)*cos(ra - ra_c));
double eta = (sin(de)*cos(de_c)-cos(de)*sin(de_c)*cos(ra - ra_c))/(sin(de)*sin(de_c)+cos(de)*cos(de_c)*cos(ra - ra_c));
rd_vector[0] = rad_to_mas(ksi)/1000;
rd_vector[1] = rad_to_mas(eta)/1000;
return rd_vector;
};

void getRaDeToTang1(double *ksi, double *eta, double ra, double de, double ra_c, double de_c)
{
//double rd_vector[2];
ra_c = mas_to_rad(ra_c);
de_c = mas_to_rad(de_c);
ra = mas_to_rad(ra);
de = mas_to_rad(de);
*ksi = cos(de)*sin(ra-ra_c)/(sin(de)*sin(de_c)+cos(de)*cos(de_c)*cos(ra - ra_c));
*eta = (sin(de)*cos(de_c)-cos(de)*sin(de_c)*cos(ra - ra_c))/(sin(de)*sin(de_c)+cos(de)*cos(de_c)*cos(ra - ra_c));
*ksi = rad_to_mas(*ksi)/1000.0;
*eta = rad_to_mas(*eta)/1000.0;
//return rd_vector;
};

void getRaDeToTang(double ra, double de, double ra_c, double de_c, double *const tang)
{
      ra_c = mas_to_rad(ra_c);
      de_c = mas_to_rad(de_c);
      ra = mas_to_rad(ra);
      de = mas_to_rad(de);
      double ksi = cos(de)*sin(ra-ra_c)/(sin(de)*sin(de_c)+cos(de)*cos(de_c)*cos(ra - ra_c));
      double eta = (sin(de)*cos(de_c)-cos(de)*sin(de_c)*cos(ra - ra_c))/(sin(de)*sin(de_c)+cos(de)*cos(de_c)*cos(ra - ra_c));
      tang[0] = rad_to_mas(ksi)/1000;
      tang[1] = rad_to_mas(eta)/1000;
};


double* getCposFromWCS(double x, double y, double *WCSD)
{
   double ksi = (x - WCSD[0])*WCSD[8] + (y - WCSD[1])*WCSD[9];
   double eta = (x - WCSD[0])*WCSD[10] + (y - WCSD[1])*WCSD[11];
   double* posC = getTangToRaDe(ksi*3600, eta*3600, WCSD[2]*3600000, WCSD[3]*3600000);
   return posC;
};

void getCposFromWCS1(double *ra, double *de, double x, double y, double *WCSD)
{
   double ksi = (x - WCSD[0])*WCSD[8] + (y - WCSD[1])*WCSD[9];
   double eta = (x - WCSD[0])*WCSD[10] + (y - WCSD[1])*WCSD[11];
   getTangToRaDe1(ra, de, ksi*3600, eta*3600, WCSD[2]*3600000, WCSD[3]*3600000);
}
/*
double* getPixPosFromWCS(double ra, double de, double *WCSD)
{
	double* tPos;
	tPos = getRaDeToTang(ra, de, WCSD[2]*3600000, WCSD[3]*3600000);
	double t1 = tPos[0]/3600; double t2 = tPos[1]/3600;
	double a1 = WCSD[8];double b1 = WCSD[9];
	double a2 = WCSD[10];double b2 = WCSD[11];
	double rPix[2];
	rPix[0] = (t2-t1*b2/b1)/(a2-a1*b2/b1);
	rPix[1] = (t1-a1*rPix[0])/b1;
	rPix[0] = rPix[0]+WCSD[0];
	rPix[1] = rPix[1]+WCSD[1];
	return rPix;
};
*/
void getPixPosFromWCS(double *x, double *y, double ra, double de, double *WCSD)
{
	double* tPos;
	tPos = new double[2];
	getRaDeToTang1(&tPos[0], &tPos[1], ra, de, WCSD[2]*3600000, WCSD[3]*3600000);

	double t1 = tPos[0]/3600.0; double t2 = tPos[1]/3600.0;
	double a1 = WCSD[8];double b1 = WCSD[9];
	double a2 = WCSD[10];double b2 = WCSD[11];
	double rPix[2];

	rPix[0] = (t2-t1*b2/b1)/(a2-a1*b2/b1);
	rPix[1] = (t1-a1*rPix[0])/b1;

	rPix[0] = rPix[0]+WCSD[0];
	rPix[1] = rPix[1]+WCSD[1];
	*x = rPix[0];
	*y = rPix[1];
};



void getPixPosFromWCS(double ra, double de, double *WCSD, double *const pixp)
{
    double tPos[2];
    getRaDeToTang(ra, de, WCSD[2]*3600000, WCSD[3]*3600000, tPos);
    double t1 = tPos[0]/3600; double t2 = tPos[1]/3600;
    double a1 = WCSD[8];double b1 = WCSD[9];
    double a2 = WCSD[10];double b2 = WCSD[11];
    pixp[0] = (t2-t1*b2/b1)/(a2-a1*b2/b1);
    pixp[1] = (t1-a1*pixp[0])/b1;
    pixp[0] = pixp[0]+WCSD[0];
    pixp[1] = pixp[1]+WCSD[1];
};

void mjdDateCode(QString *dateCode, double mJD)
{
    DATEOBS dObs;
    dObs = getDATEOBSfromMJD(mJD, 1);

    dateCode->clear();

    dateCode->append(QString("%1").arg((int)dObs.year, 4, 10, QLatin1Char( '0' )));
    dateCode->append(QString("%1").arg((int)dObs.month, 2, 10, QLatin1Char( '0' )));
    dateCode->append(QString("%1").arg((int)dObs.day, 2, 10, QLatin1Char( '0' )));
    dateCode->append(QString("%1").arg((int)dObs.hour, 2, 10, QLatin1Char( '0' )));
    dateCode->append(QString("%1").arg((int)dObs.min, 2, 10, QLatin1Char( '0' )));
    dateCode->append(QString("%1").arg((int)floor(dObs.sec*10+0.5), 3, 10, QLatin1Char( '0' )));
}

void mjdDateCode_file(QString *dateCode, double mJD)
{
    DATEOBS dObs;
    dObs = getDATEOBSfromMJD(mJD, 1);
    dateCode->clear();
    dateCode->append(QString("%1").arg((int)dObs.year, 4, 10, QLatin1Char( '0' )));
    dateCode->append(QString("%1").arg((int)dObs.month, 2, 10, QLatin1Char( '0' )));
    dateCode->append(QString("%1-").arg((int)dObs.day, 2, 10, QLatin1Char( '0' )));
    dateCode->append(QString("%1").arg((int)dObs.hour, 2, 10, QLatin1Char( '0' )));
    dateCode->append(QString("%1").arg((int)dObs.min, 2, 10, QLatin1Char( '0' )));
    //dateCode->append(QString("%1").arg((int)dObs.sec, 2, 10, QLatin1Char( '0' )));
    dateCode->append(QString("%1").arg((int)floor(dObs.sec*1000 + 0.5), 5, 10, QLatin1Char( '0' )));
}

//My
int isVes(int year)
{
    int k1, k2;
    int d1, d2, d3, d4, ym;

    k1 = k2 = 0;

    k1 = ((year/4.0)==(year/4));

    d4 = year;
    d1 = year/1000;
    d4 -= 1000*d1;
    d2 = d4/100;
    d4 -= 100*d2;
    d3 = d4/10;
    d4 -= 10*d3;

    ym = d1*10 + d2;

    k2 = k1;
    if((d4==0)&&(d3==0))
    {
        k2 = (ym/4.0)==(ym/4);
    }

    return((k1&&k2));
}

double dinm(int mounth, int isves)
{
    //mounth--;
    double m[12];
    m[0] = 31.0;
    m[1] = 28.0;
    m[2] = 31.0;
    m[3] = 30.0;
    m[4] = 31.0;
    m[5] = 30.0;
    m[6] = 31.0;
    m[7] = 31.0;
    m[8] = 30.0;
    m[9] = 31.0;
    m[10] = 30.0;
    m[11] = 31.0;
//	m[12] = m[0];

    if(mounth!=2) return(m[mounth-1]);
    else return(m[mounth-1] + isves);
}

double delta_ra(double a1, double a2)
{
	return ang_b_dir(a1,0,a2,0);
}

double grad_to_mas(double angle)
{
        return angle*MAS_IN_GRAD;
}

double grad_to_rad(double angle)
{
        return angle*atan(1)*4/180.0;
}

void rot2D(double *r, double ang)
{
        double x, y;
        x = r[0];
        y = r[1];

        r[0] = cos(ang)*x + sin(ang)*y;
        r[1] = -sin(ang)*x + cos(ang)*y;

}


//Degree

void getDegToTang(double *ksi, double *eta, double ra, double de, double ra_c, double de_c)
{
//double rd_vector[2];
ra_c = grad_to_rad(ra_c);
de_c = grad_to_rad(de_c);
ra = grad_to_rad(ra);
de = grad_to_rad(de);
*ksi = cos(de)*sin(ra-ra_c)/(sin(de)*sin(de_c)+cos(de)*cos(de_c)*cos(ra - ra_c));
*eta = (sin(de)*cos(de_c)-cos(de)*sin(de_c)*cos(ra - ra_c))/(sin(de)*sin(de_c)+cos(de)*cos(de_c)*cos(ra - ra_c));
*ksi = rad2grad(*ksi);
*eta = rad2grad(*eta);
};

void getTangToDeg(double *rd0, double *rd1, double ksi, double eta, double ra_c, double de_c)
{
        ra_c = grad_to_rad(ra_c);
        de_c = grad_to_rad(de_c);
        ksi = grad2rad(ksi);
        eta = grad2rad(eta);

        double secD = 1/cos(de_c);
        double commP = 1-eta*tan(de_c);
        double ra = atan2(ksi*secD, commP)+ra_c;
        double de = atan2((eta+tan(de_c))*cos(ra - ra_c),commP);
        *rd0 = rad2grad(ra);
        *rd1 = rad2grad(de);
};

void getPixPosToDegWCS(double *ra, double *de, double x, double y, double *WCSD)
{
   double ksi = (x - WCSD[0])*WCSD[8] + (y - WCSD[1])*WCSD[9];
   double eta = (x - WCSD[0])*WCSD[10] + (y - WCSD[1])*WCSD[11];
   getTangToDeg(ra, de, ksi, eta, WCSD[2], WCSD[3]);
}

void getDegToPixPosWCS(double *x, double *y, double ra, double de, double *WCSD)
{
        double* tPos;
        tPos = new double[2];
        getDegToTang(&tPos[0], &tPos[1], ra, de, WCSD[2], WCSD[3]);

        double t1 = tPos[0]; double t2 = tPos[1];
        double a1 = WCSD[8];double b1 = WCSD[9];
        double a2 = WCSD[10];double b2 = WCSD[11];
        double rPix[2];
        rPix[0] = (t2-t1*b2/b1)/(a2-a1*b2/b1);
        rPix[1] = (t1-a1*rPix[0])/b1;
        rPix[0] = rPix[0]+WCSD[0];
        rPix[1] = rPix[1]+WCSD[1];

        *x = rPix[0];
        *y = rPix[1];
};

double hms_to_deg(QString ra, QString spl_symb)
{
 ra = ra.trimmed();
 QString nums;
 nums = ra.section(spl_symb, 0, 0);
 double h = nums.toDouble();
 nums = ra.section(spl_symb, 1, 1);
 double m = nums.toDouble();
 nums = ra.section(spl_symb, 2, 2);
 double s = nums.toDouble();
 return (h+m/60.0+s/3600.0)*15.0;
};

int hms_to_deg(double *raDeg, QString ra, QString spl_symb)
{
     ra = ra.trimmed();
     QStringList ra_list = ra.split(spl_symb);
     if(ra_list.size()!=3) return 1;
     QString nums;
     nums = ra_list.at(0);
     double h = nums.toDouble();
     nums = ra_list.at(1);
     double m = nums.toDouble();
     nums = ra_list.at(2);
     double s = nums.toDouble();
     *raDeg = (h+m/60.0+s/3600.0)*15.0;
     return 0;
}

double damas_to_deg(QString de, QString spl_symb)
{
 de = de.trimmed();
 double sign;
 if (de[0]=='-') sign = -1; else sign = 1;
 QString nums;// = de_list.at(0);
 nums = de.section(spl_symb, 0, 0);
 double deg = nums.toDouble();
 deg = fabs(deg);
 nums = de.section(spl_symb, 1, 1);//de_list.at(1);
 double arcmin = nums.toDouble();
 nums = de.section(spl_symb, 2, 2);//de_list.at(2);
 double arcsec = nums.toDouble();
 return sign*(deg+arcmin/60+arcsec/3600.0);
};

int damas_to_deg(double *deDeg, QString de, QString spl_symb)
{
     de = de.trimmed();
     double sign;
     if (de[0]=='-') sign = -1; else sign = 1;
     QStringList de_list = de.split(spl_symb);
     if(de_list.size()!=3) return 1;
     QString nums;// = de_list.at(0);
     nums = de_list.at(0);
     double deg = nums.toDouble();
     deg = fabs(deg);
     nums = de_list.at(1);
     double arcmin = nums.toDouble();
     nums = de_list.at(2);
     double arcsec = nums.toDouble();
     *deDeg = sign*(deg+arcmin/60+arcsec/3600.0);
     return 0;
};

QString deg_to_hms(double ra, QString spl_symb, int ndig)
{
 int h = (int)(ra/15);
 ra = (ra/15.0-h)*60;
 int m = (int)(ra);
 ra = (ra-m)*60;
 double s = ra;
 QString str = QString( "%1" ).arg( h,2,10,QLatin1Char( '0' ))+spl_symb;
   str = str + QString( "%1" ).arg( m,2,10,QLatin1Char( '0' ))+spl_symb;
   str = str + QString( "%1" ).arg( s,3+ndig,'f',ndig,QLatin1Char( '0' ));
  return str;
};

void deg_to_hms(QString *str, double ra, QString spl_symb, int ndig)
{
     str->clear();
     int h = (int)(ra/15);
     ra = (ra/15.0-h)*60;
     int m = (int)(ra);
     ra = (ra-m)*60;
     double s = ra;
     str->append(QString( "%1" ).arg( h,2,10,QLatin1Char( '0' ))+spl_symb);
     str->append(QString( "%1" ).arg( m,2,10,QLatin1Char( '0' ))+spl_symb);
     str->append(QString( "%1" ).arg( s,3+ndig,'f',ndig,QLatin1Char( '0' )));
};

QString deg_to_damas(double de, QString spl_symb, int ndig)
{
int sign = 0;
if (de<0) sign = 1;
de = fabs(de);
int d = (int)(de);
de = (de - d)*60;
int am = (int)(de);
de = (de - am)*60;
double as = de;
QString str = QString( "%1" ).arg( d,2,10,QLatin1Char( '0' ))+spl_symb;
  str = str + QString( "%1" ).arg( am,2,10,QLatin1Char( '0' ))+spl_symb;
  str = str + QString( "%1" ).arg( as,3+ndig,'f',ndig,QLatin1Char( '0' ));
  if (sign==1)
    str = "-"+str;
  else
    str = "+"+str;
  return str;
};

void deg_to_damas(QString *str, double de, QString spl_symb, int ndig)
{
    str->clear();
    int sign = 0;
    if (de<0) sign = 1;
    de = fabs(de);
    int d = (int)(de);
    de = (de - d)*60;
    int am = (int)(de);
    de = (de - am)*60;
    double as = de;
    str->append(QString( "%1" ).arg( d,2,10,QLatin1Char( '0' ))+spl_symb);
    str->append(QString( "%1" ).arg( am,2,10,QLatin1Char( '0' ))+spl_symb);
    str->append(QString( "%1" ).arg( as,3+ndig,'f',ndig,QLatin1Char( '0' )));

    if(sign==1)str->prepend("-");
    else str->prepend("+");
};

double grad2rad(double grad)
{
    return(grad*PI/180.0);
}

double rad2grad(double rad)
{
    return(rad*180.0/PI);
}

double rad_to_mas(double angle)
{
        return angle*(180*3600000)/PI;
};

double rad2mas(double rad)
{
    return(rad*SECINRAD*1000.0);
}

double mas2rad(double mas)
{
    return(mas/1000.0/SECINRAD);
}

double getS(DATEOBS date_obs)
{
        double d_day = date_obs.day;
        double d_hour = date_obs.hour;
        double d_min = date_obs.min;
        double d_sec = date_obs.sec;
        d_day = d_hour*3600+d_min*60+d_sec;
        DATEOBS gm = date_obs;
        gm.hour = 0; gm.min = 0; gm.sec = 0;
        double JD = getJD(gm);
        double T=(JD-2451545)/36525;
        double S0=6*3600+41*60+50.54841+236.555367908*(JD-2451545)+0.093104*T*T-0.0000062*T*T*T;
        S0=S0-floor(S0/86400)*86400;
        double S = S0+d_day*366.25/365.25+7278;
        if (S>86400){S=S-86400;}
        return S;
};

QString getStrFromS(double s, QString spl_symb, int ndig)
{
  int h = (int)(s/3600);
  s = s-h*3600;
  int m = (int)(s/60);
  s = s-m*60;
  QString str = QString( "%1" ).arg( h,2,10,QLatin1Char( '0' ))+spl_symb;
   str = str + QString( "%1" ).arg( m,2,10,QLatin1Char( '0' ))+spl_symb;
   if (ndig!=0) str = str + QString( "%1" ).arg( s,3+ndig,'f',ndig,QLatin1Char( '0' ));
       else     str = str + QString( "%1" ).arg( int(s),2,10,QLatin1Char( '0' ));
  return str;
}

QString getRaDecMag(QString str)
{
  bool fl = false;
  for(int i=0;i<str.length();i++)
  {
        if((str[i]==' ')||(str[i]=='+')||(str[i]=='-'))
        {
           if(fl)
           {
            if ((str[i]!='+')&&(str[i]!='-'))
                {
                        str[i] = '0';
                        fl = false;
                }
           }
           else
           {
            fl = true;
           }
        }
        else fl = false;
  }
  return str;
};

QString getStrFromDE(double de, QString spl_symb)
{
        int d = (int)(de/3600);
        de = de - d*3600;
        int am = (int)(de/60);
        de = de - am*60;
        int as = (int)de;
    QString str = QString( "%1" ).arg( d,2,10,QLatin1Char( '0' ))+spl_symb;
     str = str + QString( "%1" ).arg( am,2,10,QLatin1Char( '0' ))+spl_symb;
     str = str + QString( "%1" ).arg( as,2,10,QLatin1Char( '0' ));
    return str;
};


///////////////////time_a
double UTC2TDB(double jdUTC)
{
    double jdTDT;
    jdTDT = UTC2TDT(jdUTC);
    double T=(jdTDT-2451545)/36525.0;
    double g = (357.528 + 35999.050*T)*2.0*PI/360.0;
    return(jdTDT + 0.001658*sin(g+0.0167*sin(g))/86400.0);
}

double UTC2TDT(double jdUTC)
{
    return(jdUTC + TAImUTC(jd2mjd(jdUTC)) + 32.184/86400.0);
}

double TDT2UTC(double jdTDT)
{
    return(jdTDT - TAImUTC(jd2mjd(jdTDT)) - 32.184/86400.0);
}

double TDB2TDT(double jdTDB)
{
    double T=(jdTDB-2451545)/36525.0;
    double g = (357.528 + 35999.050*T)*2.0*PI/360.0;
    double res = jdTDB - 0.001658*sin(g+0.0167*sin(g))/86400.0;
    return(res);
}

double TDT2TDB(double jdTDT)
{
    double T=(jdTDT-2451545)/36525.0;
    double g = (357.528 + 35999.050*T)*2.0*PI/360.0;
    double res = jdTDT + 0.001658*sin(g+0.0167*sin(g))/86400.0;
    return(res);
}

double TDB2UTC(double jdTDB)
{
    return(TDT2UTC(TDB2TDT(jdTDB)));
}

double TDB2TT(double jdTDB)
{
    double T=(jdTDB-2451545)/36525;
    double g = (357.528 + 35999.050*T)*2.0*PI/360.0;
    return(jdTDB - 0.001658*sin(g+0.0167*sin(g))/86400.0);
}

double dUT1() {return -0.39;}

double TAImUTC(double mjd) // [day]
{
    double res;

    if(mjd>=jd2mjd(2437300.5)&&mjd<jd2mjd(2437512.5)) res = 1.4228180 + (mjd - 37300.0)*0.001296;
    if(mjd>=jd2mjd(2437512.5)&&mjd<jd2mjd(2437665.5)) res = 1.3728180 + (mjd - 37300.0)*0.001296;
    if(mjd>=jd2mjd(2437665.5)&&mjd<jd2mjd(2438334.5)) res = 1.8458580 + (mjd - 37665.0)*0.0011232;
    if(mjd>=jd2mjd(2438334.5)&&mjd<jd2mjd(2438395.5)) res = 1.9458580 + (mjd - 37665.0)*0.0011232;
    if(mjd>=jd2mjd(2438395.5)&&mjd<jd2mjd(2438486.5)) res = 3.2401300 + (mjd - 38761.0)*0.001296;
    if(mjd>=jd2mjd(2438486.5)&&mjd<jd2mjd(2438639.5)) res = 3.3401300 + (mjd - 38761.0)*0.001296;
    if(mjd>=jd2mjd(2438639.5)&&mjd<jd2mjd(2438761.5)) res = 3.4401300 + (mjd - 38761.0)*0.001296;
    if(mjd>=jd2mjd(2438761.5)&&mjd<jd2mjd(2438820.5)) res = 3.5401300 + (mjd - 38761.0)*0.001296;
    if(mjd>=jd2mjd(2438820.5)&&mjd<jd2mjd(2438942.5)) res = 3.6401300 + (mjd - 38761.0)*0.001296;
    if(mjd>=jd2mjd(2438942.5)&&mjd<jd2mjd(2439004.5)) res = 3.7401300 + (mjd - 38761.0)*0.001296;
    if(mjd>=jd2mjd(2439004.5)&&mjd<jd2mjd(2439126.5)) res = 3.8401300 + (mjd - 38761.0)*0.001296;
    if(mjd>=jd2mjd(2439126.5)&&mjd<jd2mjd(2439887.5)) res = 4.3131700 + (mjd - 39126.0)*0.002592;
    if(mjd>=jd2mjd(2439887.5)&&mjd<jd2mjd(2441317.5)) res = 4.2131700 + (mjd - 39126.0)*0.002592;
    if(mjd>=jd2mjd(2441317.5)&&mjd<jd2mjd(2441499.5)) res = 10.0;
    if(mjd>=jd2mjd(2441499.5)&&mjd<jd2mjd(2441683.5)) res = 11.0;
    if(mjd>=jd2mjd(2441683.5)&&mjd<jd2mjd(2442048.5)) res = 12.0;
    if(mjd>=jd2mjd(2442048.5)&&mjd<jd2mjd(2442413.5)) res = 13.0;
    if(mjd>=jd2mjd(2442413.5)&&mjd<jd2mjd(2442778.5)) res = 14.0;
    if(mjd>=jd2mjd(2442778.5)&&mjd<jd2mjd(2443144.5)) res = 15.0;
    if(mjd>=jd2mjd(2443144.5)&&mjd<jd2mjd(2443509.5)) res = 16.0;
    if(mjd>=jd2mjd(2443509.5)&&mjd<jd2mjd(2443874.5)) res = 17.0;
    if(mjd>=jd2mjd(2443874.5)&&mjd<jd2mjd(2444239.5)) res = 18.0;
    if(mjd>=jd2mjd(2444239.5)&&mjd<jd2mjd(2444786.5)) res = 19.0;
    if(mjd>=jd2mjd(2444786.5)&&mjd<jd2mjd(2445151.5)) res = 20.0;
    if(mjd>=jd2mjd(2445151.5)&&mjd<jd2mjd(2445516.5)) res = 21.0;
    if(mjd>=jd2mjd(2445516.5)&&mjd<jd2mjd(2446247.5)) res = 22.0;
    if(mjd>=jd2mjd(2446247.5)&&mjd<jd2mjd(2447161.5)) res = 23.0;
    if(mjd>=jd2mjd(2447161.5)&&mjd<jd2mjd(2447892.5)) res = 24.0;
    if(mjd>=jd2mjd(2447892.5)&&mjd<jd2mjd(2448257.5)) res = 25.0;
    if(mjd>=jd2mjd(2448257.5)&&mjd<jd2mjd(2448804.5)) res = 26.0;
    if(mjd>=jd2mjd(2448804.5)&&mjd<jd2mjd(2449169.5)) res = 27.0;
    if(mjd>=jd2mjd(2449169.5)&&mjd<jd2mjd(2449534.5)) res = 28.0;
    if(mjd>=jd2mjd(2449534.5)&&mjd<jd2mjd(2450083.5)) res = 29.0;
    if(mjd>=jd2mjd(2450083.5)&&mjd<jd2mjd(2450630.5)) res = 30.0;
    if(mjd>=jd2mjd(2450630.5)&&mjd<jd2mjd(2451179.5)) res = 31.0;
    if(mjd>=jd2mjd(2451179.5)&&mjd<jd2mjd(2453736.5)) res = 32.0;
    if(mjd>=jd2mjd(2453736.5)&&mjd<jd2mjd(2454832.5)) res = 33.0;
    if(mjd>=jd2mjd(2454832.5)&&mjd<jd2mjd(2456109.5)) res = 34.0;
    if(mjd>=jd2mjd(2456109.5)) res = 35.0;

    return(res/86400.0);
}
