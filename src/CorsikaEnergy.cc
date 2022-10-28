#include <iostream>
#include <cmath>
#include "CorsikaEnergy.h"

char _ParticleName[75][30] = {   "",
                             "Photon",                    /*  1 */
                             "e+",                        /*  2 */
			     "e-",                        /*  3 */
			     "Not Defined",               /*  4 */
			     "mu+",                       /*  5 */
			     "mu-",                       /*  6 */
			     "pi0",                       /*  7 */
			     "pi+",                       /*  8 */
			     "pi-",                       /*  9 */
                             "K0l",                       /* 10 */
			     "K+",                        /* 11 */
			     "K-",                        /* 12 */
			     "n",                         /* 13 */
			     "p",                         /* 14 */
			     "Anti-p",                    /* 15 */
			     "K0S",                       /* 16 */
			     "Eta",                       /* 17 */
			     "Lambda",                    /* 18 */
			     "Sigma+",                    /* 19 */
			     "Sigma0",                    /* 20 */
			     "Sigma-",                    /* 21 */
                             "Xi0",                       /* 22 */
			     "Xi-",                       /* 23 */
			     "Omega-",                    /* 24 */
			     "Anti-n",                    /* 25 */
			     "Anti Lambda",               /* 26 */
			     "Anti-Sigma-",               /* 27 */
			     "Anti-Sigma0",               /* 28 */
			     "Anti-Sigma+",               /* 29 */
                             "Anti Xi 0",                 /* 30 */
			     "Anti Xi +",                 /* 31 */
			     "Anti-Omega +",              /* 32 */
			     "Not Defined",               /* 33 */
			     "Not Defined",               /* 34 */
			     "Not Defined",               /* 35 */
			     "Not Defined",               /* 36 */
			     "Not Defined",               /* 37 */
			     "Not Defined",               /* 38 */
			     "Not Defined",               /* 39 */
			     "Not Defined",               /* 40 */
			     "Not Defined",               /* 41 */
			     "Not Defined",               /* 42 */
			     "Not Defined",               /* 43 */
			     "Not Defined",               /* 44 */
			     "Not Defined",               /* 45 */
			     "Not Defined",               /* 46 */
			     "Not Defined",               /* 47 */
			     "Not Defined",               /* 48 */
			     "Not Defined",               /* 49 */
			     "Omega",                     /* 50 */
			     "Rho 0",                     /* 51 */
			     "Rho +",                     /* 52 */
			     "Rho -",                     /* 53 */
			     "Delta ++",                  /* 54 */
			     "Delta +",                   /* 55 */
		             "Delta 0",                   /* 56 */
     			     "Delta -",                   /* 57 */
			     "Anti-Delta--",              /* 58 */
			     "Anti-Delta-",               /* 59 */
			     "Anti-Delta0",               /* 60 */
			     "Anti-Delta+",               /* 61 */
		             "K*0",                       /* 62 */ 
		             "K*+",                       /* 63 */
			     "K*-",                       /* 64 */
                             "Anti K*0",                  /* 65 */
			     "Nu e",                      /* 66 */
			     "Anti-Nu e",                 /* 67 */
			     "Nu mu",                     /* 68 */
		             "Anti-Nu mu",                /* 69 */
		             "Not defined",               /* 70 */
		             "eta -> 2 gamma",            /* 71 */
                             "Eta -> 3 pi0",              /* 72 */
		             "Eta -> pi+ pi- p0",         /* 73 */
		             "Eta -> pi+ pi- photon"  };  /* 74 */
  
int _ParticleCharge[75] = {   	 0,	/* "",      0 */
                         0,     /* "Photon",   	    1 */
                         1,     /* "e+",	    2 */
			-1,     /* "e-", 	    3 */
			 0,    	/* "Not Defined",   4 */
			 1,    	/* "mu+",	    5 */
			-1,    	/* "mu-",	    6 */
			 0,    	/* "pi0",	    7 */
			 1,    	/* "pi+",	    8 */
			-1,    	/* "pi-",	    9 */
                         0,    	/* "K0l",	   10 */
			 1,    	/* "K+",	   11 */
			-1,     /* "K-",	   12 */
			 0,    	/* "n",		   13 */
			 1,     /* "p",		   14 */
			-1,     /* "Anti-p",	   15 */
			 0,    	/* "K0S",	   16 */
			 0,    	/* "Eta",	   17 */
			 0,    	/* "Lambda",	   18 */
			 1,    	/* "Sigma+",	   19 */
			 0,    	/* "Sigma0",	   20 */
			-1,    	/* "Sigma-",	   21 */
                         0,    	/* "Xi 0",	   22 */
			-1,     /* "X1 -",	   23 */
			-1,     /* "Omega-",	   24 */
			 0,    	/* "Anti-n",	   25 */
			 0,    	/* "Anti Lambda",  26 */
			-1,    	/* "Anti-Sigma-",  27 */
			 0,    	/* "Anti-Sigma0",  28 */
			 1,    	/* "Anti-Sigma+",  29 */
                         0,    	/* "Anti Xi 0",    30 */
			 1,    	/* "Anti Xi +",    31 */
			 1,    	/* "Anti-Omega +", 32 */
			 0,     /* "Not defined",  33 */
			 0,     /* "Not defined",  34 */
			 0,     /* "Not defined",  35 */
			 0,     /* "Not defined",  36 */
			 0,     /* "Not defined",  37 */
			 0,     /* "Not defined",  38 */
			 0,     /* "Not defined",  39 */
			 0,     /* "Not defined",  40 */
			 0,     /* "Not defined",  41 */
			 0,     /* "Not defined",  42 */
			 0,     /* "Not defined",  43 */
			 0,     /* "Not defined",  44 */
			 0,     /* "Not defined",  45 */
			 0,     /* "Not defined",  46 */
			 0,     /* "Not defined",  47 */
			 0,     /* "Not defined",  48 */
			 0,     /* "Not defined",  49 */
			 0,     /* "omega",        50 */
			 0,    	/* "Rho 0",	   51 */
			 1,    	/* "Rho +",	   52 */
			-1,     /* "Rho -",	   53 */
			 2,    	/* "Delta ++",	   54 */
			 1,    	/* "Delta +",	   55 */
		         0,    	/* "Delta 0",	   56 */
     			-1,    	/* "Delta -",	   57 */
			-2,     /* "Anti-Delta--", 58 */
			-1,     /* "Anti-Delta-",  59 */
			 0,    	/* "Anti-Delta0",  60 */
			+1,    	/* "Anti-Delta+",  61 */
		         0,    	/* "K*0",	   62 */ 
		         1,    	/* "K*+",	   63 */
			-1,     /* "K*-",	   64 */
                         0,    	/* "Anti K*0",	   65 */
			 0,    	/* "Nu e",	   66 */
			 0,    	/* "Anti-Nu e",	   67 */
		         0,    	/* "Nu mu"	   68 */
			 0,    	/* "Anti Nu mu",   69 */
		         0,    	/* "Not defined",  70 */
		         0,    	/* "eta -> 2 gamma", 71 */
                         0,    	/* "Eta -> 3 pi0",   72 */
		         0,    	/* "Eta -> pi+ pi- p0", 73 */
		         0    	/* "Eta -> pi+ pi- photon" 74*/
};


int _ParticleBaryonNumber[75] = {   	 0,	/* "", 		*/
                         0,     /* "Photon", 	*/
                         0,     /* "e+",	*/
			 0,     /* "e-", 	*/
			 0,    	/* "Not Defined",*/
			 0,    	/* "mu+",	*/
			 0,    	/* "mu-",	*/
			 0,    	/* "pi0",	*/
			 0,    	/* "pi+",	*/
			 0,    	/* "pi-",	*/
                         0,    	/* "K0l",	*/
			 0,    	/* "K+",	*/
			 0,     /* "K-",	*/
			 1,    	/* "n",		*/
			 1,     /* "p",		*/
			-1,     /* "Anti-p",	*/
			 0,    	/* "K0S",	*/
			 0,    	/* "Eta",	*/
			 1,    	/* "Lambda",	*/
			 1,    	/* "Sigma+",	*/
			 1,    	/* "Sigma0",	*/
			 1,    	/* "Sigma-",	*/
                         0,    	/* "Ksi0",	*/
			 0,     /* "Ksi-",	*/
			 0,     /* "Omega-",	*/
			-1,    	/* "Anti-n",	*/
			-1,    	/* "Anti Lambda",*/
			-1,    	/* "Anti-Sigma-",*/
			-1,    	/* "Anti-Sigma0",*/
			-1,    	/* "Anti-Sigma+",*/
                        -1,    	/* "Anti Ksi 0",*/
			-1,    	/* "Anti Ksi +",*/
			 0,    	/* "Anti-Omega +",*/
			 0,    	/* "Rho 0",	*/
			 0,    	/* "Rho +",	*/
			 0,     /* "Rho -",	*/
			 1,    	/* "Delta ++",	*/
			 1,    	/* "Delta +",	*/
		         1,    	/* "Delta 0",	*/
     			 1,    	/* "Delta -",	*/
			-1,     /* "Anti-Delta--",*/
			-1,     /* "Anti-Delta-",*/
			-1,    	/* "Anti-Delta0",*/
			-1,    	/* "Anti-Delta+",*/
		         0,    	/* "K*0",	*/ 
		         0,    	/* "K*+",	*/
			 0,     /* "K*-",	*/
                         0,    	/* "Anti K*0",	*/
			 0,    	/* "nu e",	*/
			 0,    	/* "Anti nu e",	*/
		         0,    	/* "Anti-nu mu"	*/
			 0,    	/* "nu mu",	*/
		         0,    	/* "Anti-nu mu",*/
		         0,    	/* "Not defined",*/
		         0,    	/* "eta -> 2 gamma",*/
                         0,    	/* "Eta -> 3 pi0",*/
		         0,    	/* "Eta -> pi+ pi- p0",*/
		         0    	/* "Eta -> pi+ pi- photon" */
};

//Unit of particle mass is GeV
float _ParticleMass[75] = {0,            /* "",      0 */
                         0,            /* "Photon",        1 */
                         0.511e-3,     /* "e+",            2 */
                         0.511e-3,     /* "e-",            3 */
                         0,            /* "Not Defined",   4 */
                         105.7e-3,     /* "mu+",           5 */
                         105.7e-3,     /* "mu-",           6 */
                         135.0e-3,     /* "pi0",           7 */
                         139.6e-3,     /* "pi+",           8 */
                         139.6e-3,     /* "pi-",           9 */
                         497.6e-3,     /* "K0l",          10 */
                         493.7e-3,     /* "K+",           11 */
                         493.7e-3,     /* "K-",           12 */
                         939.6e-3,     /* "n",            13 */
                         938.3e-3,     /* "p",            14 */
                         938.3e-3,     /* "Anti-p",       15 */
                         497.6e-3,     /* "K0S",          16 */
                         547.9e-3,     /* "Eta",          17 */
                         1116.e-3,     /* "Lambda",       18 */
                         1189.e-3,     /* "Sigma+",       19 */
                         1193.e-3,     /* "Sigma0",       20 */
                         1197.e-3,     /* "Sigma-",       21 */
                         1315.e-3,     /* "Xi 0",         22 */
                         1322.e-3,     /* "X1 -",         23 */
                         1672.e-3,     /* "Omega-",       24 */
                         938.3e-3,     /* "Anti-n",       25 */
                         1116.e-3,     /* "Anti Lambda",  26 */
                         1189.e-3,     /* "Anti-Sigma-",  27 */
                         1193.e-3,     /* "Anti-Sigma0",  28 */
                         1197.e-3,     /* "Anti-Sigma+",  29 */
                         1315.e-3,     /* "Anti Xi 0",    30 */
                         1322.e-3,     /* "Anti Xi +",    31 */
                         1672.e-3,     /* "Anti-Omega +", 32 */
                         0,            /* "Not defined",  33 */
                         0,            /* "Not defined",  34 */
                         0,            /* "Not defined",  35 */
                         0,            /* "Not defined",  36 */
                         0,            /* "Not defined",  37 */
                         0,            /* "Not defined",  38 */
                         0,            /* "Not defined",  39 */
                         0,            /* "Not defined",  40 */
                         0,            /* "Not defined",  41 */
                         0,            /* "Not defined",  42 */
                         0,            /* "Not defined",  43 */
                         0,            /* "Not defined",  44 */
                         0,            /* "Not defined",  45 */
                         0,            /* "Not defined",  46 */
                         0,            /* "Not defined",  47 */
                         0,            /* "Not defined",  48 */
                         0,            /* "Not defined",  49 */
                         782.7e-3,     /* "omega",        50 */
                         775.5e-3,     /* "Rho 0",        51 */
                         775.5e-3,     /* "Rho +",        52 */
                         775.5e-3,     /* "Rho -",        53 */
                         1232.e-3,     /* "Delta ++",     54 */
                         1232.e-3,     /* "Delta +",      55 */
                         1232.e-3,     /* "Delta 0",      56 */
                         1232.e-3,     /* "Delta -",      57 */
                         1232.e-3,     /* "Anti-Delta--", 58 */
                         1232.e-3,     /* "Anti-Delta-",  59 */
                         1232.e-3,     /* "Anti-Delta0",  60 */
                         1232.e-3,     /* "Anti-Delta+",  61 */
                         895.9e-3,     /* "K*0",          62 */
                         891.7e-3,     /* "K*+",          63 */
                         891.7e-3,     /* "K*-",          64 */
                         895.9e-3,     /* "Anti K*0",     65 */
                         0,     /* "Nu e",         66 */
                         0,     /* "Anti-Nu e",    67 */
                         0,     /* "Nu mu"         68 */
                         0,     /* "Anti Nu mu",   69 */
                         0,     /* "Not defined",  70 */
                         0,     /* "eta -> 2 gamma", 71 */
                         0,     /* "Eta -> 3 pi0",   72 */
                         0,     /* "Eta -> pi+ pi- p0", 73 */
                         0      /* "Eta -> pi+ pi- photon" 74*/
};

float GetPartcileMass(int id)
{
  float mass=0;
  if(id<75)
    mass=_ParticleMass[id]; 
  if(id/100>200)
    mass=id/100*931.5;
  return mass;
}
float GetParticleEnergy(int id,float px, float py,float pz)
{
  float energy;
  float mass=GetPartcileMass(id);
  energy=sqrt(px*px+py*py+pz*pz+mass*mass);
  return energy;
}
