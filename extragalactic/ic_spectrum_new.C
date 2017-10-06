//************************************
// SSC code   Author: Lab Saha
// version 1.0.3
// Date:  15/11/2014
// Last Modified: 31/07/2015 
//***********************************
#include<iostream>
#include<cmath>
#include <fstream>
#include <complex>
//#include "spectrum.h"   
#include <cstdio>
#include <limits>
#include <cstdlib>
//#include "test1.h"

//#ifndef _NR_H_
//#define _NR_H_
//#endif /* _NR_H_ */
//*****************
extern int FLAG; 
extern int FLAG_SSA; 
extern double B;
extern double Gamma;
extern double zshift;
extern double es;
extern double dL;
extern double deltaD;
extern double mec2;
extern double Rb;
extern double gammaMax;
extern double PI;
extern double erg_to_eV;
extern double h;
extern double c;
extern double e;
extern double dm;
extern double dmc2;
extern double dsigt;
extern double dN;
extern double dp;
extern double dq;
extern double GammaMax;
extern double GammaMin;
extern double GammaBreak;
extern double GammaBinSize;
extern int nGammaBin;
extern int n_bin;
extern double nseed[50];
extern double dNdg[50];
extern double gammaN[50];
extern double OmegaM;
extern double Omega_Lambda;
double func_bessel(int,double, double,double );
double FunctTheta(double );
double  func_gamma(double );
double func_simpsint_1F(double , double , int,double );
double func_simpsint_ic(double , double , int,double );
double FuncSimpson1D(double , double , int,double,double );
double func_simpint2D_1D_1F(double , double , double , double,double );
double func_gamma_bessel(double );
double Simpint1D_IC(double , double , double , double, double );
double funct_ic(int , double,int,double);
double func_ssc(double);
double func_Syn_seed(double);
double absorption_coeff(double , double );
double luminosity_dist_cal(double , double );
double etotal_cal(double , int , double, double);
double ndensity_cal(double, int , double , double);
double fn_tau_dx(double, double,double);
double fn_tau_de(double,double,double,double);
double fn_tau(double, double, double,double);
double extragal_absorption(double,double);
//***************************************
int    FLAG=0;
int    FLAG_SSA=0;
double B = 2.2e-2;
double Gamma = 10;
double zshift = 0.25;
double es = 1.0;
double dL = 1.0;
double deltaD=1.0;
double mec2=8.01e-7;
double Rb = 5.9000e+16; // in cm
double gammaMax = 5.9000e+6; //
double PI=3.1415927;
double erg_to_eV=1.0/1.602e-12;
double h=4.14131e-15;
double c = 3.0e10;
double e = 4.803e-10;
double dm = 9.109e-28;
double dmc2=.51e6;
double dsigt=6.65e-25;
double dN =1.0;
double dp =1;
double dq =1;
double GammaMin=1.0;
double GammaMax=1.0;
double GammaBreak =1.0;
double nseed[50]={0};
double dNdg[50]={0};
double gammaN[50]={0};
double GammaBinSize=1.0;
int nGammaBin=50;
int n_bin=50;
double OmegaM=0.3;
double Omega_Lambda=0.7;
//***********************
using namespace std;

template <class T>
class NRVec {
private:
	int nn;	// size of array. upper index is nn-1
	T *v;
public:
	NRVec();
	explicit NRVec(int n);		// Zero-based array
	NRVec(const T &a, int n);	//initialize to constant value
	NRVec(const T *a, int n);	// Initialize to array
	NRVec(const NRVec &rhs);	// Copy constructor
	NRVec & operator=(const NRVec &rhs);	//assignment
	NRVec & operator=(const T &a);	//assign a to every element
	inline T & operator[](const int i);	//i'th element
	inline const T & operator[](const int i) const;
	inline int size() const;
	~NRVec();
};
template <class T>
NRVec<T>::NRVec(const T *a, int n) : nn(n), v(new T[n])
{
	for(int i=0; i<n; i++)
		v[i] = *a++;
}


template <class T>
inline const T & NRVec<T>::operator[](const int i) const	//subscripting
{
	return v[i];
}


template <class T>
NRVec<T>::~NRVec()
{
	if (v != 0)
		delete[] (v);
}


typedef double DP;
typedef NRVec<DP> Vec_DP, Vec_O_DP, Vec_IO_DP;
typedef const NRVec<DP> Vec_I_DP;
namespace NR {
	inline void nrerror(const string error_text)
	// Numerical Recipes standard error handler
	{
		cerr << "Numerical Recipes run-time error..." << endl;
		cerr << error_text << endl;
		cerr << "...now exiting to system..." << endl;
		exit(1);
	}


void beschb(const DP x, DP &gam1, DP &gam2, DP &gampl, DP &gammi);
DP bessi(const int n, const DP x);
DP bessi0(const DP x);
DP bessi1(const DP x);
void bessik(const DP x, const DP xnu, DP &ri, DP &rk, DP &rip, DP &rkp);
DP chebev(const DP a, const DP b, Vec_I_DP &c, const int m, const DP x);
DP gammln(const DP xx);
}

//*****************************

double func_gamma_bessel(double e1);

int main(){

      double signal, sigerr, dsqrt_Sp, dms;
      char character;
      char line[1024];
      double e1max,e1min,binsize,ene_l,ene_h,area,dL;
      double dis,e1,e1p,beta,theta,Fac;
      int i,NumArgc;
      double dL_Mpc, energy_density;
      double E_total,vol;
      double number_density;
      double GammaMinLog,GammaMaxLog,GammaBreakLog;

//      NumArgc = IARGC()
//      IF (NumArgc < 1) THEN
//      write(*,'(a)') "usage:  ./foo outFileName"
//      call exit(0)
//      endif

//      call GETARG(1,InFile)
//      read(*,*) InFile
//      open(unit=12,file=InFile,status='unknown')

//************** Reading input parameters ****************
      FILE *fp=fopen("InputICSSC","r");
      fgets(line,1024,fp);
      fscanf(fp,"%d \n",&FLAG);
      fgets(line,1024,fp);
      fscanf(fp,"%d \n",&FLAG_SSA);
      fgets(line,1024,fp);
      fscanf(fp,"%lf\n",&B);
      fgets(line,1024,fp);
      fscanf(fp,"%lf %lf\n",&dp,&dq);
      fgets(line,1024,fp);
      fscanf(fp,"%lf\n",&energy_density);
      fgets(line,1024,fp);
      fscanf(fp,"%lf %lf\n",&e1min,&e1max);
      fgets(line,1024,fp);
      fscanf(fp,"%lf\n",&Rb);
      fgets(line,1024,fp);
      fscanf(fp,"%lf %lf\n",&GammaMinLog,&GammaMaxLog);
      fgets(line,1024,fp);
      fscanf(fp,"%lf \n",&GammaBreakLog);
      fgets(line,1024,fp);
      fscanf(fp,"%lf \n",&Gamma);
      fgets(line,1024,fp);
      fscanf(fp,"%lf \n",&zshift);
      fgets(line,1024,fp);
      fscanf(fp,"%lf \n",&theta);
      fgets(line,1024,fp);
      fscanf(fp,"%lf \n",&OmegaM);
      fgets(line,1024,fp);
      fscanf(fp,"%lf \n",&Omega_Lambda);
//************************************************
      GammaMin   = exp(GammaMinLog);
      GammaMax   = exp(GammaMaxLog);
      GammaBreak = exp(GammaBreakLog);
//  
//********************************
       binsize = (GammaMaxLog - GammaMinLog)/nGammaBin;
      GammaBinSize = binsize;
      
      for(int i=0;i<nGammaBin;i++){
      gammaN[i] = exp(GammaMinLog+i*binsize);
       if(FLAG == 0){
       if(gammaN[i] >= GammaMin && gammaN[i] <= GammaMax){
       dNdg[i] = pow(gammaN[i],-dp);}
       else {       dNdg[i] = 0.;}
        }
       if(FLAG == 1){
       if(gammaN[i] >= GammaMin && gammaN[i] <= GammaBreak){
       dNdg[i] = pow(gammaN[i],-dp);}
       else if(gammaN[i] > GammaBreak && gammaN[i] <= GammaMax){
        dNdg[i] = pow(GammaBreak,-dp+dq)*pow(gammaN[i],-dq);}
       else {       dNdg[i] = 0.;}
        }
}
//********************************



//********************************
//      cout<<" DEBUG " <<GammaMin<< "  "<<GammaMax<<endl;

      E_total = etotal_cal(1.0,FLAG,GammaMin, GammaMax);
//      cout<< " E_total " <<E_total<<endl;
      dN = energy_density/E_total;
      number_density=ndensity_cal(dN,FLAG,GammaMin, GammaMax);
      cout<<" Number Density(number/cm^3): "<<number_density<<endl;
      dL_Mpc = luminosity_dist_cal(OmegaM,Omega_Lambda);
      dL = dL_Mpc*1e6*3.0860*1e18; // luminosity  distance in cm
      area = 4.0*PI*pow(dL,2);
//    Corresponding  Factor for blazar  //
       beta   = sqrt(pow(Gamma,2)-1)/Gamma;  
       deltaD = 1.0/(Gamma * (1- beta* cos(theta*PI/180.0)));
        
       
       Fac = pow(deltaD,4); 
//************************************** 
//*     
        double fsyn,e2;
        double y0,yn,dy,y,zq,yp;
        y0=1e-8;
        yn=1e8;
        dy= (log(yn)-log(y0))/n_bin;
        for(int i=0;i<n_bin;i++)
	{
          y=exp(log(y0)+i*dy);
          nseed[i]=erg_to_eV*func_gamma_bessel(y)*2.24/(4*PI*pow(Rb,2)*c*h);

        // 2.24 comes due to integration over isotropic emissitivity
        // see Atoyan and Aharonian  1996 MNRAS 278, 525-541 (1996)
        }

//
        double sum0,sum1,sum2,sum3,sum,x,z,f;
        float signal_SynIC,signal_Syn;
        double signalIC_ebl,signal_SynIC_ebl, tau;
//        
	binsize =(log(e1max) - log(e1min))/n_bin;
        for(int i=1;i<n_bin;i++)
	{
         ene_l = log(e1min)+i*binsize;
         ene_h = log(e1min)+(i+1)*binsize;
         e1 = (exp(ene_h)+exp(ene_l))/2.00;
//        
	 e1p = e1*(1+zshift)/deltaD;  
        
         signal_Syn = func_gamma_bessel(e1p);
	
	 signal = 4./3*PI*pow(Rb,3)*func_ssc(e1p);
        
	 tau = extragal_absorption(e1,zshift);
        
	 signalIC_ebl = signal*exp(-tau);

        
	 signal_SynIC = (erg_to_eV*1e-12*Fac*(e1p/h)*signal_Syn/area) + 1e-12*e1p*Fac*signal/area; 
        
	 signal_SynIC_ebl = (erg_to_eV*1e-12*Fac*(e1p/h)*signal_Syn/area) + 1e-12*e1p*Fac*signalIC_ebl/area; 
        
// NOTE: e1 should always be in eV.
// We can divide it by area to get TeV/cm^2 unit.

         cout<< (e1/h) <<" " << e1 <<"  " << erg_to_eV*1e-12*Fac*(e1p/h)*signal_Syn/area <<"   "<< (1.0e-12*e1p*Fac*signal/area)<<"  " << signal_SynIC <<"  "<< 1e-12*e1p*Fac*signalIC_ebl/area << "  " << signal_SynIC_ebl << endl;

        
 //   ***********************
 //   Bolometric luminosity
 //   //
 //
       z = log(e1);
       x = exp(z);
       f = x * signal_SynIC_ebl/pow(e1,2);

       if(i==0){ sum0=sum0+exp(z)*f;}
       if(i%2==1){ sum1=sum1+4*exp(z)*f;}
       if(i%2==0 && i!=n_bin-1 && i!=0){sum2=sum2+2*exp(z)*f;}
       if(i==n_bin-1){sum3=sum3+exp(z)*f;}
        
 
        }

        sum=(sum0+sum1+sum2+sum3)*binsize/3;

        cout<<" Bolometric Luminosity(erg/s) " << sum*area*1.602 << endl;
//************************

}

//***************

double etotal_cal(double norm, int FLAG, double x0, double xn)
{
  double f,sum = 0,sum0=0,sum1=0,sum2=0,sum3=0,nud,dNt;
  double x_min=x0,x_max=xn,gamma,y0,yn;
  double z_min,z_max,dz;
  double kk;
  for(int j=0;j<nGammaBin;j++){
  f=gammaN[j]*norm*dNdg[j]*dmc2/erg_to_eV;
//  cout<<"F "<<f<<"  "<<gammaN[j]<<"  "<<dNdg[j]<<endl;
  kk=j%2;
  if(j==0){ sum0=sum0+gammaN[j]*f;}
  if(j%2==1){ sum1=sum1+4*gammaN[j]*f;}
  if(j%2==0 && j!=nGammaBin-1 && j!=0){sum2=sum2+2*gammaN[j]*f;}
  if(j==nGammaBin-1){sum3=sum3+gammaN[j]*f;}
//  cout<<j<<" " <<func_gamma(x)<<endl;
  }
  sum=(sum0+sum1+sum2+sum3)*GammaBinSize/3;
//  cout<<"sum "<< sum0<<" "<<sum1<<" "<<sum2<<" "<<sum3<<endl;
//  cout<<sum0<<"  "<<sum1<<endl;
  return sum;
}
//************************

double ndensity_cal(double norm, int FLAG, double x0, double xn)
{
  double f,sum = 0,sum0=0,sum1=0,sum2=0,sum3=0,nud,dNt;
  double x_min=x0,x_max=xn,gamma,y0,yn;
  double z_min,z_max,dz;
  double kk;
  for(int j=0;j<nGammaBin;j++){

  f=norm*dNdg[j];
  kk=j%2;
  if(j==0){ sum0=sum0+gammaN[j]*f;}
  if(j%2==1){ sum1=sum1+4*gammaN[j]*f;}
  if(j%2==0 && j!=nGammaBin-1 && j!=0){sum2=sum2+2*gammaN[j]*f;}
  if(j==nGammaBin-1){sum3=sum3+gammaN[j]*f;}
//  cout<<j<<" " <<func_gamma(x)<<endl;
  }
  sum=(sum0+sum1+sum2+sum3)*GammaBinSize/3;
//  cout<<"sum "<< sum0<<" "<<sum1<<" "<<sum2<<" "<<sum3<<endl;

  return sum;
}
//************************


//*******
void NR::beschb(const DP x, DP &gam1, DP &gam2, DP &gampl, DP &gammi)
{
	const int NUSE1=7, NUSE2=8;
	static const DP c1_d[7] = {
		-1.142022680371168e0,6.5165112670737e-3,
		3.087090173086e-4,-3.4706269649e-6,6.9437664e-9,
		3.67795e-11,-1.356e-13};
	static const DP c2_d[8] = {
		1.843740587300905e0,-7.68528408447867e-2,
		1.2719271366546e-3,-4.9717367042e-6,-3.31261198e-8,
		2.423096e-10,-1.702e-13,-1.49e-15};
	DP xx;
	static Vec_DP c1(c1_d,7),c2(c2_d,8);

	xx=8.0*x*x-1.0;
	gam1=chebev(-1.0,1.0,c1,NUSE1,xx);
	gam2=chebev(-1.0,1.0,c2,NUSE2,xx);
	gampl= gam2-x*gam1;
	gammi= gam2+x*gam1;
}
//**********
#include <cmath>
#include <limits>
//#include "nr.h"
using namespace std;

void NR::bessik(const DP x, const DP xnu, DP &ri, DP &rk, DP &rip, DP &rkp)
{
	const int MAXIT=10000;
	const DP EPS=numeric_limits<DP>::epsilon();
	const DP FPMIN=numeric_limits<DP>::min()/EPS;
	const DP XMIN=2.0, PI=3.141592653589793;
	DP a,a1,b,c,d,del,del1,delh,dels,e,f,fact,fact2,ff,gam1,gam2,
		gammi,gampl,h,p,pimu,q,q1,q2,qnew,ril,ril1,rimu,rip1,ripl,
		ritemp,rk1,rkmu,rkmup,rktemp,s,sum,sum1,x2,xi,xi2,xmu,xmu2;
	int i,l,nl;

	if (x <= 0.0 || xnu < 0.0) nrerror("bad arguments in bessik");
	nl=int(xnu+0.5);
	xmu=xnu-nl;
	xmu2=xmu*xmu;
	xi=1.0/x;
	xi2=2.0*xi;
	h=xnu*xi;
	if (h < FPMIN) h=FPMIN;
	b=xi2*xnu;
	d=0.0;
	c=h;
	for (i=0;i<MAXIT;i++) {
		b += xi2;
		d=1.0/(b+d);
		c=b+1.0/c;
		del=c*d;
		h=del*h;
		if (fabs(del-1.0) <= EPS) break;
	}
	if (i >= MAXIT)
		nrerror("x too large in bessik; try asymptotic expansion");
	ril=FPMIN;
	ripl=h*ril;
	ril1=ril;
	rip1=ripl;
	fact=xnu*xi;
	for (l=nl-1;l >= 0;l--) {
		ritemp=fact*ril+ripl;
		fact -= xi;
		ripl=fact*ritemp+ril;
		ril=ritemp;
	}
	f=ripl/ril;
	if (x < XMIN) {
		x2=0.5*x;
		pimu=PI*xmu;
		fact = (fabs(pimu) < EPS ? 1.0 : pimu/sin(pimu));
		d = -log(x2);
		e=xmu*d;
		fact2 = (fabs(e) < EPS ? 1.0 : sinh(e)/e);
		beschb(xmu,gam1,gam2,gampl,gammi);
		ff=fact*(gam1*cosh(e)+gam2*fact2*d);
		sum=ff;
		e=exp(e);
		p=0.5*e/gampl;
		q=0.5/(e*gammi);
		c=1.0;
		d=x2*x2;
		sum1=p;
		for (i=1;i<=MAXIT;i++) {
			ff=(i*ff+p+q)/(i*i-xmu2);
			c *= (d/i);
			p /= (i-xmu);
			q /= (i+xmu);
			del=c*ff;
			sum += del;
			del1=c*(p-i*ff);
			sum1 += del1;
			if (fabs(del) < fabs(sum)*EPS) break;
		}
		if (i > MAXIT) nrerror("bessk series failed to converge");
		rkmu=sum;
		rk1=sum1*xi2;
	} else {
		b=2.0*(1.0+x);
		d=1.0/b;
		h=delh=d;
		q1=0.0;
		q2=1.0;
		a1=0.25-xmu2;
		q=c=a1;
		a = -a1;
		s=1.0+q*delh;
		for (i=1;i<MAXIT;i++) {
			a -= 2*i;
			c = -a*c/(i+1.0);
			qnew=(q1-b*q2)/a;
			q1=q2;
			q2=qnew;
			q += c*qnew;
			b += 2.0;
			d=1.0/(b+a*d);
			delh=(b*d-1.0)*delh;
			h += delh;
			dels=q*delh;
			s += dels;
			if (fabs(dels/s) <= EPS) break;
		}
		if (i >= MAXIT) nrerror("bessik: failure to converge in cf2");
		h=a1*h;
		rkmu=sqrt(PI/(2.0*x))*exp(-x)/s;
		rk1=rkmu*(xmu+x+0.5-h)*xi;
	}
	rkmup=xmu*xi*rkmu-rk1;
	rimu=xi/(f*rkmu-rkmup);
	ri=(rimu*ril1)/ril;
	rip=(rimu*rip1)/ril;
	for (i=1;i <= nl;i++) {
		rktemp=(xmu+i)*xi2*rk1+rkmu;
		rkmu=rk1;
		rk1=rktemp;
	}
	rk=rkmu;
	rkp=xnu*xi*rkmu-rk1;
}
//****************
DP NR::chebev(const DP a, const DP b, Vec_I_DP &c, const int m, const DP x)
{
	DP d=0.0,dd=0.0,sv,y,y2;
	int j;

	if ((x-a)*(x-b) > 0.0)
		nrerror("x not in range in routine chebev");
	y2=2.0*(y=(2.0*x-a-b)/(b-a));
	for (j=m-1;j>0;j--) {
		sv=d;
		d=y2*d-dd+c[j];
		dd=sv;
	}
	return y*d-dd+0.5*c[0];
}
//*******************
// using simpson's 1/3 rule for 2d integration by calling 1D twice

//double func_simpint2D_1D_1F(double x0, double xn, double dx, double xbin, double e1)
//{

//}
//****************
// using simpson's 1/3 rule for 2d integration

double func_simpsint_1F(double y0, double yn, int ll, double e1)
{
  double f, f1, sum = 0,sum0=0,sum1=0,sum2=0,sum3=0;
  double x_min=y0,x_max=yn,x,z0,zn,nud;
  double z_min,z_max,dz;
//  int n_bin=50;
  double z[n_bin],kk;
  z_min=log(x_min);
  z_max=log(x_max);
  dz= (z_max-z_min)/n_bin;

//  for ( int i=0;i<n_bin;i++){
//  z[i]=z_min+i*dz;
//  }
  for(int j=0;j<n_bin;j++){
  z[j]= z_min+j*dz;
  x=exp(z[j]);
  nud=3.0*e*B*sin(x)/(4.0*PI*dm*c);
  z0 = e1/(h*nud*pow(gammaN[ll],2));    // Bessel function variabl
  zn = 300.0;              // for Synchrotron photons
//  cout<<" e1 "<<e1<<" y0 " <<y0<<" e "<<e<<" yn:: "<<yn<<" h:: "<<h<<endl; 
  if (z0 < zn){
//  cout<<" e1 "<<e1<<" y0 " <<z0<<" e "<<e<<" x:: "<<x<<" h:: "<<h<<endl; 
  f=FuncSimpson1D(z0,zn,ll,x,e1);
  kk=j%2;
  if(j==0){ sum0=sum0+exp(z[j])*f;}
  if(j%2==1){ sum1=sum1+4*exp(z[j])*f;}
  if(j%2==0 && j!=n_bin-1 && j!=0){sum2=sum2+2*exp(z[j])*f;}
  if(j==n_bin-1){sum3=sum3+exp(z[j])*f;}
//  cout<<j<<" " <<kk<<endl;
  }}
  sum=(sum0+sum1+sum2+sum3)*dz/3;
//  cout<<"sum2 "<< sum<<endl;
  
  return sum;
}
//********************************
#include<iostream>
#include<cmath>
//#include "spectrum.h"
//#include "nr.h"
double func_bessel(int pp, double theta1, double y,double e1)
{
      double I_bessel;
      double rk,ri,rip,rkp,FIFTHRD;
      double nud,dNt,dconst;
//*     ==========All Constants============================== 
      FIFTHRD=5.0/3.00;//      ! Fractional order of Bessel Function
//*     ============== Some Required function definition =============
//*    ============= Call Bessel function (rk gives the value) =======
      NR::bessik(y,FIFTHRD,ri,rk,rip,rkp);    

//*     =============== Main Function definition =======================
//      nud=3.00*e*B*sqrt(2.00/3.0)/(4.0*PI*dm*c);
      nud=3.00*e*B*sin(theta1)/(4.0*PI*dm*c);
      dconst=sqrt(3.0)*pow(e,3)*B/(dm*c*c);
//      dconst=sqrt(2.0)*pow(e,3)*B/(dm*c*c)*sqrt(3.0/2);
//      cout<<" dN "<<dN<< endl;         

       double kappa_nu,fac,tau;
       kappa_nu=absorption_coeff(gammaN[pp],e1);
       tau = 2 * kappa_nu* Rb;
 
       if (FLAG_SSA ==0){       
       fac=4.0/3*Rb;
       I_bessel=fac*PI*pow(Rb,2)*dN*dNdg[pp]*(e1/(nud*h))*pow(gammaN[pp],-2)*dconst*rk*sin(theta1)*sin(theta1)/2;// ! For Electron Energy Integration
       }
       else{
       if(kappa_nu <=0){ fac=1.0;}
       else {
       if(tau<1E-3) {   fac = 4./3*Rb ;}
       else { 
//       tau = 2 * kappa_nu* Rb; 
       fac = (1.-2./pow(tau,2)*(1.-exp(-tau)*(tau+1.)))/kappa_nu;
//       cout<<kappa_nu<<endl;
            }}
//* * * * * * * * * * * * * * * * * * * * * *
//   
       I_bessel=fac*PI*pow(Rb,2)*dN*dNdg[pp]*(e1/(nud*h))*pow(gammaN[pp],-2)*dconst*rk*sin(theta1)*sin(theta1)/2;// ! For Electron Energy Integration
          }
//*     ================================================================      
      return I_bessel;
      
}
//**********************      
#include<iostream>
#include<math.h>
//#include "spectrum.h"
using namespace std;
//double func_simpint2D_1D_1F(double x0, double xn, double dx, double xbin);
double func_gamma_bessel(double e1)
{
  double y0,yn;
  double f,sum = 0,sum0=0,sum1=0,sum2=0,sum3=0,nud;
  int kk;

  y0=0.001;
  yn=PI;
  for(int j=0;j<nGammaBin;j++){
  //  x=exp(z[j]);
  //  cout<<" e1 "<<e1<<" y0 " <<y0<<" e "<<e<<" x:: "<<x<<" h:: "<<h<<endl; 
  f=func_simpsint_1F(y0,yn,j,e1);
  kk=j%2;
  if(j==0){ sum0=sum0+gammaN[j]*f;}
  if(j%2==1){ sum1=sum1+4*gammaN[j]*f;}
  if(j%2==0 && j!=nGammaBin-1 && j!=0){sum2=sum2+2*gammaN[j]*f;}
  if(j==nGammaBin-1){sum3=sum3+gammaN[j]*f;}
  //  cout<<j<<" " <<func_gamma(x)<<endl;
  }
  sum=(sum0+sum1+sum2+sum3)*GammaBinSize/3;
  //  cout<<"sum "<< sum0<<" "<<sum1<<" "<<sum2<<" "<<sum3<<endl;
  return sum;
}
//*************************
// using simpson's 1/3 rule for 2d integration
double FuncSimpson1D(double z0, double zn, int pp,double theta, double e1)
{
  double f, f1, sum = 0,sum0=0,sum1=0,sum2=0,sum3=0;
  double x_min=z0,x_max=zn,x,x1;
  double z_min,z_max,dz;
//  int n_bin=100;
  double z[n_bin],kk;
  z_min=log(x_min);
  z_max=log(x_max);
  dz= (z_max-z_min)/n_bin;
//  for ( int i=0;i<n_bin;i++){
//  z[i]=z_min+i*dz;
//  }
  for(int j=0;j<n_bin;j++){
  z[j] = z_min+j*dz;
  x=exp(z[j]);
  f=func_bessel(pp,theta,x,e1);
  kk=j%2;
  if(j==0){ sum0=sum0+exp(z[j])*f;}
  if(j%2==1){ sum1=sum1+4*exp(z[j])*f;}
  if(j%2==0 && j!=n_bin-1 && j!=0){sum2=sum2+2*exp(z[j])*f;}
  if(j==n_bin-1){sum3=sum3+exp(z[j])*f;}
//  cout<<j<<" " <<kk<<endl;
  }
  sum=(sum0+sum1+sum2+sum3)*dz/3;
//  cout<<"sum2 "<< sum<<endl;
  
  return sum;
}
//***********************************
#include<iostream>
#include<cmath>
//#include "spectrum.h"
using namespace std;
      
double funct_ic(int mm, double epsilon,int kk, double e1)
{
      double dconst,dist,I_funct_ic;
      double done_ovr_4G,done_min_gmc; 
        
      double Fic;
      double q;
       
      double seed_const,KT; 
      double dN_r1;
//      double ga=gammaN[mm];
      
//      dist=dis*3.086e21;//  ! Distance to the source in cm 
      KT=8.617343*1.e-5 * 2.7;
       
//*    ===== One Part of constants ============================
      dconst=3.0*c*dsigt/4.0;

//*     ===============================================================
      done_min_gmc = (1.0-e1/(gammaN[mm]*dmc2));
//c      print*,done_min_gmc,nseed
      if (done_min_gmc == 0.0  ){
      I_funct_ic=0.0;}
      
      else {
   
//*     ============== Some Required function definition =============
      q=e1/(4.00*epsilon*pow(gammaN[mm],2)*done_min_gmc);
      Fic=2.0*q*log(q)+(1.0+2.0*q)*(1.0-q)+.50*(1.0-q)*pow(4.0*epsilon*gammaN[mm]*q/dmc2,2)/(1.0+(4.0*epsilon*gammaN[mm]*q/dmc2));
      
//      cout<< "Fic= "<<Fic<<" q  " <<q<<endl;
//*    
//*      
//*    ========== Condition Check ====================================
      done_ovr_4G = 1.0/4.0/pow(gammaN[mm],2);
//c      print*,q,done_ovr_4G
      if (q <= done_ovr_4G || q >= 1.0){
      I_funct_ic=0.0;}
       else {
//*    ===============================================================
//*    ============= Call Bessel function (rk gives the value) =======
//c      call bessik(y,FIFTHRD,ri,rk,rip,rkp)   ! y should be dimensionless 
//c    ===============================================================
       
//*     =============== Main Function definition ======================= 
//       I_funct_ic=e1*dN*dNdg[mm]*pow(gammaN[mm],-2)*(nseed[kk]/epsilon)*Fic*dconst;                          
       I_funct_ic=e1*dN*dNdg[mm]*pow(gammaN[mm],-2)*(nseed[kk]/pow(epsilon,2))*Fic*dconst;                          
//       cout<<"funct_ic "<< I_funct_ic<<" dN_r1  "<< dconst<<endl;

       }}
         return I_funct_ic;
      
}
//**************************      
// using simpson's 1/3 rule for 2d integration by calling 1D twice
//double Simpint1D_IC(double x0, double xn, double dx, double xbin,double e1)
//{

//}
//***********************
// using simpson's 1/3 rule for 2d integration
double func_simpsint_ic(double y0, double yn, int mm,double e1)
{
  double f, f1, sum = 0,sum0=0,sum1=0,sum2=0,sum3=0;
  double x_min=y0,x_max=yn,x,z0,zn,nud;
  double z_min,z_max,dz;
//  int n_bin=100;
  double z[n_bin],kk;
  z_min=log(x_min);
  z_max=log(x_max);
  dz= (z_max-z_min)/n_bin;
//  for ( int i=0;i<n_bin;i++){
//  z[i]=z_min+i*dz;
//  }
  for(int j=0;j<n_bin;j++){
  z[j] = z_min+j*dz;
  x=exp(z[j]);
//  cout<<" e1 "<<e1<<" y0 " <<z0<<" e "<<e<<" x:: "<<x<<" h:: "<<h<<endl; 
  f=funct_ic(mm,x,j,e1);
  kk=j%2;
  if(j==0){ sum0=sum0+exp(z[j])*f;}
  if(j%2==1){ sum1=sum1+4*exp(z[j])*f;}
  if(j%2==0 && j!=n_bin-1 && j!=0){sum2=sum2+2*exp(z[j])*f;}
  if(j==n_bin-1){sum3=sum3+exp(z[j])*f;}
//  cout<<j<<" " <<kk<<endl;
  }
  sum=(sum0+sum1+sum2+sum3)*dz/3;
//  cout<<"sum2 "<< sum<<endl;
  
  return sum;
}

//**********************************
#include<iostream>
#include<math.h>
//#include "spectrum.h"
using namespace std;
//double func_simpint2D_1D_1F(double x0, double xn, double dx, double xbin);
double func_ssc(double e1)
{
  double y0,yn;  
  double f,sum = 0,sum0=0,sum1=0,sum2=0,sum3=0,nud;
  double kk;

  y0=1e-8;
  yn=1e8; 
  for(int j=0;j<nGammaBin;j++){
//  cout<<" e1 "<<e1<<" y0 " <<y0<<" e "<<e<<" x:: "<<x<<" h:: "<<h<<endl; 
  f=func_simpsint_ic(y0,yn,j,e1);
  kk=j%2;
  if(j==0){ sum0=sum0+gammaN[j]*f;}
  if(j%2==1){ sum1=sum1+4*gammaN[j]*f;}
  if(j%2==0 && j!=nGammaBin-1 && j!=0){sum2=sum2+2*gammaN[j]*f;}
  if(j==nGammaBin-1){sum3=sum3+gammaN[j]*f;}
//  cout<<j<<" " <<func_gamma(x)<<endl;
  }
  sum=(sum0+sum1+sum2+sum3)*GammaBinSize/3;
//  cout<<"sum "<< sum0<<" "<<sum1<<" "<<sum2<<" "<<sum3<<endl;
  
  return sum;

}
//**********************
double func_Syn_seed(double epsilon)

{
      double signal, sigerr, dsqrt_Sp, dms;
      char character;
      char line[1024];
      double e1max,e1min,binsize,ene_l,ene_h,area,dist;
      double dis;
      int i,n_bin,NumArgc;
    

         signal = erg_to_eV*func_gamma_bessel(epsilon)/(h*epsilon*4*PI*pow(Rb,2)*c);
 
         return signal;


}

//**********************
double absorption_coeff(double ga, double e1)
{
//int main(){
        
       
      double kappa_nu,fac1,fac2;
      double p,dN1=0.0;
//      ga = 100.0;
//      p = 2.05;

       if (FLAG == 0) { 
       p=dp;
       dN1=dN;}
//       cout<<" DEBUG0 "<< ga<<" "<<GammaMin<<endl; 
       if (FLAG == 1) { 
       if(ga >= GammaMin && ga <= GammaBreak){
       p= dp;
       dN1 = dN;}
//       cout <<" TEST " <<dN1<<" "<<" p " <<p<<endl;}
       else if(ga > GammaBreak && ga <=GammaMax){
       p= dq;
       dN1 = dN*pow(GammaBreak,-dp+dq);}}      


      fac1 = dN*pow(e,2)*sqrt(PI)*pow(3,(p+1)/2)*pow(e*B/(2*PI*dm*c),(p+2)/2)/(8*dm*c);


//      fac1 = pow(3*e/(2*PI*pow(dm,3)*pow(c,4)),p/2);
      fac2 = exp(NR::gammln((3*p+22)/12))*exp(NR::gammln((3*p+2)/12))*exp(NR::gammln((p+6)/4))/exp(NR::gammln((p+8)/4));
//      kappa_nu = sqrt(3*PI)*pow(e,3)*dN1*pow(B,(p+2)/2)*c*fac1*fac2*pow((e1/h),-(p+4)/2)/(16*PI*dm); 
      kappa_nu =fac1*fac2*pow((e1/h),-(p+4)/2); 

         return kappa_nu;
//         cout<<"Kappa_nu "<<kappa_nu<<endl;
      
}
//**************
DP NR::gammln(const DP xx)
{
	int j;
	DP x,y,tmp,ser;
	static const DP cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,0.1208650973866179e-2,
		-0.5395239384953e-5};

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<6;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}
//***********************    

double luminosity_dist_cal(double Omega_M,double Omega_Lambda)
{
	// Reference: Wright, E. L., 2006, PASP, 118, 1711.
	int n=1000;	// number of points in integrals
	double c = 299792.458; // velocity of light in km/sec
	
//	double H0 = h0*100;	// Hubble constant
        double H0 = 69.6;	
	// Densities
	double WM = Omega_M;		// Omega(matter)
	double WV = Omega_Lambda;	// Omega(vacuum) or lambda
	double h = H0/100;	// H0/100
	double WR = 4.165E-5/(h*h);	// Omega(radiation), includes 3 massless neutrino species, T0 = 2.72528
	double WK = 1-WM-WR-WV;	// Omega curvaturve = 1-Omega(total)
	
	double a = 1.0;	// the scale factor of the Universe
	double az = 1.0/(1+1.0*zshift);
	
//	double age = 0;
	double adot;

        double Tyr = 977.8; // coefficent for converting 1/H into Gyr
        double  DTT = 0.5;	// time from z to now in units of 1/H0
	double DTT_Gyr = 0.0;	// value of DTT in Gyr
        double age = 0.5;	// age of Universe in units of 1/H0
	double age_Gyr = 0.0;	// value of age in Gyr
        double zage = 0.1;	// age of Universe at redshift z in units of 1/H0
	double zage_Gyr = 0.0;	// value of zage in Gyr
        double DCMR = 0.0;	// comoving radial distance in units of c/H0
	double DCMR_Mpc = 0.0;
	double DCMR_Gyr = 0.0;
        double DA = 0.0;	// angular size distance
	double DA_Mpc = 0.0;
	double DA_Gyr = 0.0;
	double kpc_DA = 0.0;
        double DL = 0.0;	// luminosity distance
	double DL_Mpc = 0.0;
	double DL_Gyr = 0.0;	// DL in units of billions of light years
        double V_Gpc = 0.0;
	for (int i = 0; i < n; i++) 
	{
		a = az*((double)i+0.5)/(double)n;
		adot = sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a));
		age = age + 1/adot;
	}
	zage = az*age/n;
	// correction for annihilations of particles not present now like e+/e-
	// added 13-Aug-03 based on T_vs_t.f
	double lpz = log((1+1.0*zshift))/log(10.0);
	double dzage = 0;
	if (lpz >  7.500) dzage = 0.002 * (lpz -  7.500);
	if (lpz >  8.000) dzage = 0.014 * (lpz -  8.000) +  0.001;
	if (lpz >  8.500) dzage = 0.040 * (lpz -  8.500) +  0.008;
	if (lpz >  9.000) dzage = 0.020 * (lpz -  9.000) +  0.028;
	if (lpz >  9.500) dzage = 0.019 * (lpz -  9.500) +  0.039;
	if (lpz > 10.000) dzage = 0.048;
	if (lpz > 10.775) dzage = 0.035 * (lpz - 10.775) +  0.048;
	if (lpz > 11.851) dzage = 0.069 * (lpz - 11.851) +  0.086;
	if (lpz > 12.258) dzage = 0.461 * (lpz - 12.258) +  0.114;
	if (lpz > 12.382) dzage = 0.024 * (lpz - 12.382) +  0.171;
	if (lpz > 13.055) dzage = 0.013 * (lpz - 13.055) +  0.188;
	if (lpz > 14.081) dzage = 0.013 * (lpz - 14.081) +  0.201;
	if (lpz > 15.107) dzage = 0.214;
	zage = zage*pow(10.0,dzage);
	zage_Gyr = (Tyr/H0)*zage;

        DTT = 0.0;	
	// do integral over a=1/(1+z) from az to 1 in n steps, midpoint rule
	for (int i = 0; i < n; i++) {
		a = az+(1.-az)*((double)i+0.5)/(double)n;
		adot = sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a));
		DTT = DTT + 1./adot;
		DCMR = DCMR + 1./(a*adot);
	}
	
	DTT = (1-az)*DTT/(double)n;
	DCMR = (1-az)*DCMR/(double)n;
	
	age = DTT+zage;
	age_Gyr = age*(Tyr/H0);
	
	
	DTT_Gyr = (Tyr/H0)*DTT;
	
	DCMR_Gyr = (Tyr/H0)*DCMR;
	DCMR_Mpc = (c/H0)*DCMR;
	
	
	// tangential comoving distance
	double ratio = 1.00;
	double x;
	double y;
	x = sqrt(fabs(WK))*DCMR;
	if (x > 0.1) 
	{
		if (WK>0) ratio = 0.5*(exp(x)-exp(-x))/x;
		else ratio = sin(x)/x;
		y = ratio*DCMR;
	}
	else	
	{
		y = x*x;
		// statement below fixed 13-Aug-03 to correct sign error in expansion
		if (WK < 0) y = -y;
		ratio = 1 + y/6 + y*y/120;
		y= ratio*DCMR;
	}
	double DCMT = y;
	
	DA = az*DCMT;
	DA_Mpc = (c/H0)*DA;
	kpc_DA = DA_Mpc/206.264806;
	DA_Gyr = (Tyr/H0)*DA;
	DL = DA/(az*az);
	DL_Mpc = (c/H0)*DL;
	DL_Gyr = (Tyr/H0)*DL;
	
//	printf("Distance parameters: age_Gyr %e zage_Gyr %e DTT_Gyr %e DA_Mpc %e kpc_DA %e DL_Mpc %e \n", 
//		   age_Gyr,zage_Gyr,DTT_Gyr,DA_Mpc,kpc_DA,DL_Mpc);
	
	return DL_Mpc;
}
//*******************************************
double extragal_absorption(double Egamma,double ze)
{
   //Egamma has to be in eV
   
   
   ifstream fp("Gamma-gamma-opacity-z=0-2.dat");

   double ebl_TeV[2000][50],ebl_eV[2000][50],tau[2000][50],exp_tau[2000][50];
   double y1,y2,y3,y4,y_result;
   double u,v,t;
//   std::string data;
   double z[2000],en,dz=0.001;
   int l,m,n,i,j,p,count=0;
  
   if (fp.is_open())
        {
	
       for(int i=0;i<2000;i++)
         {
         z[i]=0+(i+1)*dz;
         std::string data;
             getline(fp, data);
//             cout<<data<<endl;
             getline(fp, data);
//             cout<<data<<endl;
             
             for(int j=0;j<50;j++)
             {
               fp >> ebl_TeV[i][j] >> ebl_eV[i][j] >> tau[i][j] >>exp_tau[i][j];
//               cout<<j<<"  "<<z[i]<<"  "<<ebl_TeV[count][j]<<endl;
             }
           getline(fp, data);
         }

//   double ze=z;
      en = Egamma;
      if( en < ebl_eV[0][0] || en > ebl_eV[0][49] ) return 0.0;
 
        int zbin = (int)(ze/dz);
        if (zbin>2000-2) zbin=2000-2;

        double z1    = dz*(double)zbin;
        double z2    = z1+dz;
        t = (ze-z1)/(z2-z1);

        for(l=0;l<2000-1;l++)
        if(ze == z[l]){
          for(m=0;m<50-2;m++)
          if(en > ebl_eV[l][m] && en <ebl_eV[l][m+1]) break;
          
          y_result = tau[l][m] + (tau[l][m+1]-tau[l][m])*(en-ebl_eV[l][m])/(ebl_eV[l][m+1]-ebl_eV[l][m]);
//          cout<<ebl_eV[l][m]<<"  "<<(en-ebl_eV[l][m])<<"  "<<y_result<<endl;
//          exit(1);
          return y_result;
                      }  
        
//      for(l=0;l<2000-1;l++)
//      if(ze >= z[l] && ze < z[l+1]) break;

//   cout<<l <<" "<< z[l]<< "  "<<en<<endl;
   for(m=0;m<50-2;m++)
   if(en < ebl_eV[zbin][m]) break;

//   for(n=0;n<31-2;n++)
//   if(en <e_ebl[l+1][n]) break;

//   cout<<"l "<<l<< " m "<<m<<" n "<<n<<"  "<<e_ebl[l+1][n]<<"  "<<e_ebl[l+1][n-1]<<endl;
   y1 = tau[zbin][m];
   y2 = tau[zbin+1][m];
   y3 = tau[zbin+1][m+1];
   y4 = tau[zbin][m+1];
//
//   t = (ze-z[l])/(z[l+1]-z[l]);
   u = (en - ebl_eV[zbin][m])/(ebl_eV[zbin][m+1] - ebl_eV[zbin][m]);
//   v = (en-e_ebl[l+1][n])/(e_ebl[l+1][n+1]-e_ebl[l+1][n]);
//   cout<< t << "  "<<u<<"  "<<v<<endl;

   y_result = (1 - t) * (1 - u) * y1 + t * (1 - u) * y2 + t * u * y3 + (1 - t) * u * y4;
 
   return y_result;
   }
   {
    cout << "Unable to open file :: 'Gamma-gamma-opacity-z=0-2.dat'"<<endl;
    cout << "EBL correction won't be there "<<endl;
    return 0.0;
   } 


}

//*******************************************
/* // This part is not required since tau value is calculated already at different z (0-2range) and different energies. It is taken from "http://www.astro.unipd.it/background/"
 
//for integaring over z and getting tau value
double extragal_absorption(double Egamma,double ze, double ebl)
{
//integration over  z
double z_min=0,dz,z,f;
double sum0=0.0,sum1=0.0,sum2=0.0,sum3=0.0,sum;
  dz = (ze-z_min)/n_bin;
 
  for(int j=0;j<n_bin;j++){
  z = z_min+j*dz;
  f=fn_tau_dx(Egamma,z,ebl);

  if(j==0){ sum0=sum0+f;}
  if(j%2==1){ sum1=sum1+4*f;}
  if(j%2==0 && j!=n_bin-1 && j!=0){sum2=sum2+2*f;}
  if(j==n_bin-1){sum3=sum3+f;}
//  cout<<j<<" " <<kk<<endl;
  }
  sum=(sum0+sum1+sum2+sum3)*dz/3;
//  cout<<"sum2 "<< sum<<endl;
 return sum;
}
// *******************************************
//for integration over x
double fn_tau_dx(double Egamma,double z,double ebl)
{
double sum0=0.0,sum1=0.0,sum2=0.0,sum3=0.0,sum;
double x_min=0,x_max=2,dx,x,f;
dx = (x_max-x_min)/n_bin;

  for(int j=0;j<n_bin;j++){
  x = x_min+j*dx;
  f=fn_tau_de(Egamma,z,x,ebl);

  if(j==0){ sum0=sum0+f;}
  if(j%2==1){ sum1=sum1+4*f;}
  if(j%2==0 && j!=n_bin-1 && j!=0){sum2=sum2+2*f;}
  if(j==n_bin-1){sum3=sum3+f;}
//  cout<<j<<" " <<kk<<endl;
  }
  sum=(sum0+sum1+sum2+sum3)*dx/3;
//  cout<<"sum2 "<< sum<<endl;

  return sum;
}
// ****************************
// for integration over energy 
double fn_tau_de(double Egamma,double z, double x, double ebl)
{
double z1_min,z1_max,dz1,z1[n_bin],e,f;
double sum0=0.0,sum1=0.0,sum2=0.0,sum3=0.0,sum; 
double x_max=1.27,x_min;
x_min=2*dmc2*dmc2/(2*Egamma*ebl*x*(1+z));
z1_min=log(x_min);
z1_max=log(x_max);
dz1=(z1_max-z1_min)/n_bin;

  for(int j=0;j<n_bin;j++){
  z1[j] = z1_min+j*dz1; //z is for mking it to logarithmic bin
  e=exp(z1[j]);
  f=fn_tau(Egamma,e,z,x,ebl);

  if(j==0){ sum0=sum0+exp(z1[j])*f;}
  if(j%2==1){ sum1=sum1+4*exp(z1[j])*f;}
  if(j%2==0 && j!=n_bin-1 && j!=0){sum2=sum2+2*exp(z1[j])*f;}
  if(j==n_bin-1){sum3=sum3+exp(z1[j])*f;}
//  cout<<j<<" " <<kk<<endl;
  }
  sum=(sum0+sum1+sum2+sum3)*dz1/3;
//  cout<<"sum2 "<< sum<<endl;
  return sum;
}
// ******************************

double fn_tau(double Egamma, double e, double ze, double xe,double ebl)
{
double dtdz,sigma_gg,s,beta;
double tau;
double H0 = 69.6;
s=2*Egamma*ebl*xe*(1+ze);
beta = sqrt(1-4*dmc2*dmc2/s);
sigma_gg= (3.*dsigt/16) *(1-pow(beta,2))*(2*beta*(pow(beta,2)-2)+(3-pow(beta,4))*log((1+beta)/(1-beta)));
dtdz=(1.0/(H0*(1+ze)))*pow(pow(1+ze,2)*(1+OmegaM*ze)-ze*(ze+2)*Omega_Lambda,-0.5);
//dndgebl();
//tau = c*dtdz*0.5*xe*dngdebl*sigma_gg;
return tau;
}

//
double ebl_photon_density(double ze,double en)
//
{
FILE *fp=fopen("franceschini_new.dat","r");
   float z[11] ={0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0};
   float e_ebl[11][31];
   float den_ebl[11][31];
   float ya[11][31];
   double y1,y2,y3,y4,y_result;
   double u,v,t;
   double ze,en;
   int l,m,n,i,j,p;
   for(int j=0;j<31;j++){
   for(int i=0;i<11;i++){
   fscanf(fp,"%e %e",&e_ebl[i][j],&den_ebl[i][j]);
   ya[i][j]=den_ebl[i][j];
   }}
//   ze=1.1;
   en=log10(en);
//   for(j=0;j<31;j++) 
//    en[j]=e_ebl[9][j];

      for(l=0;l<11-1;l++)
      if(ze > z[l] && ze < z[l+1]) break;

//   cout<<l <<" "<< z[5]<< "  "<<en<<endl;
   for(m=0;m<31-2;m++)
   if(en < e_ebl[l][m]) break;
   for(n=0;n<31-2;n++)
   if(en <e_ebl[l+1][n]) break;

//   cout<<"l "<<l<< " m "<<m<<" n "<<n<<"  "<<e_ebl[l+1][n]<<"  "<<e_ebl[l+1][n-1]<<endl;
   y1 = ya[l][m];
   y2 = ya[l][m+1];
   y3 = ya[l+1][n];
   y4 = ya[l+1][n+1];
//
   t = (ze-z[l])/(z[l+1]-z[l]);
   u = (en-e_ebl[l][m])/(e_ebl[l][m+1]-e_ebl[l][m]);
   v = (en-e_ebl[l+1][n])/(e_ebl[l+1][n+1]-e_ebl[l+1][n]);
//   cout<< t << "  "<<u<<"  "<<v<<endl;

   y_result = (1-t)*(1-u)*y1+t*(1-v)*y3+t*v*y4+(1-t)*u*y2;

   return pow(10,y_result)/pow(10,en);

}
*/ // This part is not required since tau value is calculated already at different z (0-2range) and different energies. It is taken from "http://www.astro.unipd.it/background/"

 
