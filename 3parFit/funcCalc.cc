#include<iostream> 
#include<iomanip> 
#include<fstream> 
#include<stdlib.h> 
#include<stdio.h> 
#include<math.h> 
#include<fcntl.h> 
#include <cstdlib>

#include "funcCalc.h"



double sv123(const double &t, const double &t01, const double &tb1, const double &t02, 
		const double &tb2, const double &td1, const double &ts1){
	
    double sv123 = 0.;
	double  dks0, dks1, dksm,
			dw0, dw1, dwp, dwm, das1, dac1, das0, dac0, dzna, dksm2, ds, dd,
			dcs0, dsn0, dzn0, td, ts, dr,
			dcs0s, dsn0s, dcs0d, dsn0d, dcs1s, dsn1s, dcs1d, dsn1d;


	if (t < 0.) return 0.;

	dr = (ts1 - td1) / td1;
	if (fabs(dr) >= 0.00001) {
		td = td1;
		ts = ts1;
	} else {
		td = td1;
		if (ts1 > td1) {
		ts = td1 * 1.00001;
		} else {
		ts = td1 * 0.99999;
		}
	}



	dr = ((t01 - t02) * (t01 - t02) + (tb1 - tb2) * (tb1 - tb2)) / ((t01) * (t01) + (tb1) * (tb1));
	dks0 = 1.0 / t01;
	dks1 = 1.0 / t02;

	if (dr < 0.0000000001) {

		if (dks0 > dks1) {
		dks0 = dks1 * 1.00001;
		} else {
		dks0 = dks1 * 0.99999;
		}
	}

	//  printf(" ts=%e dr=%e dks0=%e dks1=%e \n",ts,dr,dks0,dks1);

	if (t < 0.) return 0;



	dksm = dks1 - dks0;

	ds = 1. / ts;
	dd = 1. / td;

	dw0 = 1. / tb1;
	dw1 = 1. / tb2;
	dwp = dw0 + dw1;
	dwm = dw1 - dw0;

	dksm2 = dksm * dksm;

	dzna = (dksm2 + dwm * dwm) * (dksm2 + dwp * dwp);


	das0 = dw1 * (dksm2 + dwp * dwm);
	dac0 = -2 * dksm * dw0 * dw1;
	das1 = dw0 * (dksm2 - dwp * dwm);
	dac1 = -dac0;





	dsn0 = (ds - dks0);
	dcs0 = -dw0;
	dzn0 = dcs0 * dcs0 + dsn0 * dsn0;

	dsn0s = (dsn0 * das0 - dcs0 * dac0) / dzn0;
	dcs0s = (dcs0 * das0 + dsn0 * dac0) / dzn0;

	dsn0 = (ds - dks1);
	dcs0 = -dw1;
	dzn0 = dcs0 * dcs0 + dsn0 * dsn0;

	dsn1s = (dsn0 * das1 - dcs0 * dac1) / dzn0;
	dcs1s = (dcs0 * das1 + dsn0 * dac1) / dzn0;


	dsn0 = (dd - dks0);
	dcs0 = -dw0;
	dzn0 = dcs0 * dcs0 + dsn0 * dsn0;

	dsn0d = (dsn0 * das0 - dcs0 * dac0) / dzn0;
	dcs0d = (dcs0 * das0 + dsn0 * dac0) / dzn0;

	dsn0 = (dd - dks1);
	dcs0 = -dw1;
	dzn0 = dcs0 * dcs0 + dsn0 * dsn0;

	dsn1d = (dsn0 * das1 - dcs0 * dac1) / dzn0;
	dcs1d = (dcs0 * das1 + dsn0 * dac1) / dzn0;

	sv123 = ((((dsn0s - dsn0d) * sin(dw0 * t)
				+ (dcs0s - dcs0d) * cos(dw0 * t)) * exp(-t * dks0)
				- (dcs0s + dcs1s) * exp(-t * ds) + (dcs0d + dcs1d) * exp(-t * dd)
				+ ((dsn1s - dsn1d) * sin(dw1 * t)
				+ (dcs1s - dcs1d) * cos(dw1 * t)) * exp(-t * dks1)) / dzna / (ts - td));

			
	sv123=sv123/(-.109+.919*t01-.261*t01*t01)
		/(-.109+.919*t02-.261*t02*t02)
		/(.262+.174*tb1-.208*tb1*tb1)
		/(.262+.174*tb2-.208*tb2*tb2)
		/(4.56-1.58*td1)/(1.391-0.434*ts1)
		/(1.06-0.578*(t01-tb1)*(t01-tb1))
		/(1.06-0.578*(t02-tb2)*(t02-tb2))
		/(1.2140-0.79645*t01+0.63440*t01*t01)
		/(1.2140-0.79645*t02+0.63440*t02*t02);

	return sv123;

}

double ShaperDSP_F(double *xx, double *ss){

	double tr1 = xx[0];
	//double precision FITPAD(12)

	double FITFUN;

	//      common/norm/ped,amp,ts0,td,t0,b1,t1,amm,tmm,t01,b2,a
	double ped,amp,ts0,td,t0,b1,t1,amm,tmm,t01,b2,a;
	double tr,x;
	double tr2, tr3;

	
	ped = 0.0;
	amp = * (ss + 10);
	ts0 = * (ss + 0);
	td  = * (ss + 1);
	t0  = * (ss + 2);
	b1  = * (ss + 3);
	t1  = * (ss + 4);

	amm = * (ss + 5);
	tmm = * (ss + 6);

	t01 = * (ss + 7);
	b2  = * (ss + 8);
	a   = * (ss + 9);


	tr=tr1-ts0;
	tr2=tr+0.2;
	//**      
	tr3=tr-0.2;
	//**

	if(tr2 <= 0) return (ped ); 

	FITFUN=( sv123(tr,t0,b1,t01,b2,td,t1)*(1-a) + 
			a*0.5*(sv123(tr2,t0,b1,t01,b2,td,t1) + sv123(tr3,t0,b1,t01,b2,td,t1)) );
	x=tr/t0;

	FITFUN=amp*(FITFUN-amm*(exp(-tr/tmm)*(1-exp(-x)*(1+x+x*x/2+x*x*x/6+x*x*x*x/24+x*x*x*x*x/120))));


	FITFUN = FITFUN + ped;

	return FITFUN;
}
