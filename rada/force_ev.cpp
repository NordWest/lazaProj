#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rada.h"
#include <dele.h>
#include <omp.h>

#define OMPLIB

#define ka 0.017202098955

  extern double a, CC, omega, Ltilde, A;
  //extern int nofzbody;
  extern int SK, CENTER, centr_num;
  extern int useEPM;

  extern ever_params *eparam;
  int iterNum;

  extern int nofzbody;
  extern double *mass;
  extern dele eph;
  extern double *massF;
  extern int *fbodynum;
  extern int nfbody;

double dist3(double X0[], double X1[])
{
    return(sqrt(pow(X1[0] - X0[0], 2) + pow(X1[1] - X0[1], 2) + pow(X1[2] - X0[2], 2)));
}

double norm3(double *v)
{
    double nm = 0.0;
    int len;
    len = 3;

    for(int i = 0; i<len; i++)
    {
            nm += v[i]*v[i];
    }

    return(sqrt(nm));
}


void force_N(double X[], double V[], double F[]);
int force_PPN(double X[], double V[], double F[]);
void force_GN(double X[], double V[], double F[]);
void force_GN_dele(double X[], double V[], double TS, double F[]);


void Everhardt::force(double X[], double V[], double TS, double F[])
{
    iterNum = 0;

    force_GN_dele(X, V, TS, F);
}
//GELIOCENTR
void force_GN(double X[], double V[], double F[])
{
  int iNum = nofzbody;
  int Ni = iNum*3;

  #pragma omp parallel for
  for(int teloi=0; teloi<iNum; teloi++)
  {
      int i=teloi*3;
      double Ri = norm3(&X[i]);

      if(Ri>(eparam->vout))
      {
          printf("WARN!!!! V OUT!!!!\n");
          printf("Ri[%d]: %f > %f\n", teloi, Ri, eparam->vout);
          exit(1);
      }

      double massI = mass[teloi];//0.0;
      if(massI<0)massI=0;
      for(int komp=0; komp<3; komp++)
      {
          double res0, res1;
              res0 = res1 = 0.0;
//                         #pragma omp parallel for reduction(+:res0)
              for(int teloj=0; teloj<iNum; teloj++)
              {
                 int j=teloj*3;
                 if(teloi!=teloj&&mass[teloj]>0)
                 {
                    double Rij = dist3(&X[i], &X[j]);
                    double Rj = norm3(&X[j]);

                    if(Rij<eparam->col)
                    {

                        printf("teloi= %d\tteloj= %d\n", teloi, teloj);
                        printf("Rij= %f\n", Rij);
                        printf("WARN!!!! CRASH!!!!\n");
                        exit(1);
                    }



                    res0 += mass[teloj]*((X[j+komp] - X[i+komp])/(pow(Rij,3)) - X[j+komp]/(pow(Rj, 3)));

                 }

              }

              res1 = -((1.0 + massI)*X[i+komp])/(pow(Ri, 3));


              F[i+komp] = ka*ka*(res0+res1);

          }
  }
}

void force_GN_dele(double X[], double V[], double TS, double F[])
{
  int iNum = nofzbody;
  int Ni = iNum*3;
  int fNum = nfbody;
  int Nf = fNum*3;
  double *XF = new double[Nf];

  for(int i=0; i<fNum; i++)
  {
      eph.detR(&XF[i*3], &XF[i*3+1], &XF[i*3+2], TS, fbodynum[i], 0, CENTER, SK);
  }

  #pragma omp parallel for
  for(int teloi=0; teloi<iNum; teloi++)
  {
      int i=teloi*3;
      double Ri = norm3(&X[i]);

      if(Ri>(eparam->vout))
      {
          printf("WARN!!!! V OUT!!!!\n");
          printf("Ri[%d]: %f > %f\n", teloi, Ri, eparam->vout);
          exit(1);
      }

      double massI = mass[teloi];//0.0;
      if(massI<0)massI=0;
      for(int komp=0; komp<3; komp++)
      {
          double res0, res1;
              res0 = res1 = 0.0;
//                         #pragma omp parallel for reduction(+:res0)
              for(int teloj=0; teloj<iNum; teloj++)
              {
                 int j=teloj*3;
                 if(teloi!=teloj&&mass[teloj]>0)
                 {
                    double Rij = dist3(&X[i], &X[j]);
                    double Rj = norm3(&X[j]);

                    if(Rij<eparam->col)
                    {

                        printf("teloi= %d\tteloj= %d\n", teloi, teloj);
                        printf("Rij= %f\n", Rij);
                        printf("WARN!!!! CRASH!!!!\n");
                        exit(1);
                    }



                    res0 += mass[teloj]*((X[j+komp] - X[i+komp])/(pow(Rij,3)) - X[j+komp]/(pow(Rj, 3)));

                 }

              }

              for(int telof=0; telof<fNum; telof++)
              {
                 int fn=telof*3;
                 if(massF[telof]>0)
                 {
                    double Rij = dist3(&X[i], &XF[fn]);
                    double Rj = norm3(&XF[fn]);

                    if(Rij<eparam->col)
                    {

                        printf("teloi= %d\ttelof= %d\n", teloi, telof);
                        printf("Rij= %f\n", Rij);
                        printf("WARN!!!! CRASH!!!!\n");
                        exit(1);
                    }



                    res0 += massF[telof]*((XF[fn+komp] - X[i+komp])/(pow(Rij,3)) - XF[fn+komp]/(pow(Rj, 3)));

                 }

              }

              res1 = -((1.0 + massI)*X[i+komp])/(pow(Ri, 3));


              F[i+komp] = ka*ka*(res0+res1);

          }
  }
}



/*==============================================================*/
/* EOF 								*/
/*==============================================================*/
