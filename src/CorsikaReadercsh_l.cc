//*****************************************************************
// minicorread c++ version
// Function: read out the corsika output data file
// Author: Zhiguo Yao, 2008/10/13
// Written on base of a fortran version of Min Zha & Zhiguo Yao
//*****************************************************************

//#define _THIN_
//#define _CERENKOV_

#ifdef _CERENKOV_
#define _DEFAULT_FILE_ "CER00000X"
#else
#define _DEFAULT_FILE_ "DAT00000X"
#endif

//The length of a particle record (in unit of "word" = 4 bytes)
#ifdef _THIN_
#define _PARTICLE_LENGTH_ 8
#else
#define _PARTICLE_LENGTH_ 7
#endif

//Number of particles per subblock
#define _NPARTICLE_ 39

//Number of subblocks per record
#define _NSUBBLOCK_ 21

//The length of a subblock
#define _SUBBLOCK_LENGTH_   _PARTICLE_LENGTH_*_NPARTICLE_

//The length of a record
#define _RECORD_LENGTH_     _SUBBLOCK_LENGTH_*_NSUBBLOCK_

#define PI 3.14159265358979312
#define RADDEG 180/PI

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

#include <TMath.h>
#include "CorsikaEnergy.h"

using namespace std;

void Getmean_sigam(double offset_x,vector<double> &xs_,vector<double> &es_,double *mean,double *sigma){
  double sumx=0,sume=0;
  for(int j=0;j<xs_.size();j++){
    sumx+=(xs_.at(j)-offset_x)*es_.at(j);
    sume+=es_.at(j);
  }
  double meanx=sumx/sume;
  sumx=0;sume=0;
  for(int j=0;j<xs_.size();j++){
    sumx+=pow((xs_.at(j)-offset_x-meanx)*es_.at(j),2);
    sume+=pow(es_.at(j),2);
  }
  double sigmax=sqrt(sumx/sume);
  *mean=meanx;
  *sigma=sigmax;
  //printf(" ---- meanx %f\n",meanx,sigmax);
}

namespace cordata {
	//RUNH
	int irun;
	int idate;
	float version;
	int nlevel;

	//EVTH
	int ievent, priPID;
	int ipartp;
	float ep;
	float pxp, pyp, pzp;
	float thetap, phip;

	//EVTE
	int nphoton, nelectron, nhadron, nmuon, nparticle;

	//RUNE
	int nevent;

	//Particle bank
#ifdef _CERENKOV_
	float nclight;
	float x, y, u, v, t, height, weight;
#else
	int ipart;
	int igen;
	int ilevel, secPID;
	float px, py, pz, x, y, t, weight;
#endif
};


//int minicooread_rms(char* inputfile = _DEFAULT_FILE_, char* outputfile= _DEFAULT_FILE_, char* outputfile2= _DEFAULT_FILE_){
int minicooread_rms(char* inputfile = _DEFAULT_FILE_, char* outputfile= _DEFAULT_FILE_){
	double zenith, azimuth;
	int i, j,  k, temp, temp0, Nevent,ninbunch,istep, priPID, secPID;
	double x0, y0, z0,x1, y1, z1, xp, yp, zp,xc, yc, radius;
	double m, n, l, theta, n1, m1, l1,t,h,Xmax,Nmax,nclight;
	double time,the,phy,r,dens;
	double totalclight,partclight[20000],x0_sh, px, py, pz ;
	//  float Nclight, thetap, phip, ep, xs, ys, thetas, phis, es, hs, pxs, pys, pzs, xs2, ys2; 
	double Nclight, thetap, phip, ep, xs, ys, thetas, phis, es, hs, pxs, pys, pzs, xs2, ys2,weight; 
	//double roteang=1.04837;
	double roteang=0.0;
int A,q;
	int nsec=0;
	totalclight=0.;q=0.;dens=0.;
        int ievt;
	double corex[20]={0};double corey[20]={0};
	double sumex = 0,sumey = 0,sume = 0;
	double offset_x,offset_y,offset_x1,offset_y1,offset_x2,offset_y2;
	char Name0[500]="root://eos01.ihep.ac.cn/";
	char Name1[500];
//	strcpy(Name1,Name0);
//	strcat(Name1,outputfile);
	//TFile *newfile = TFile::Open(Name1,"recreate");
int iter=0;
vector<double>secPID_;
vector<double>xs_;
vector<double>ys_;
vector<double>es_;
vector<double>hs_;
vector<double>thetas_;
vector<double>phis_;
vector<double>pxs_;
vector<double>pys_;
vector<double>pzs_;

	TFile *file= TFile::Open(outputfile,"recreate");
	TTree *events = new TTree("events","RayTrace");
	events->Branch("priPID",&priPID,"priPID/I");
	events->Branch("ep",&ep,"ep/D");
	events->Branch("thetap",&thetap,"thetap/D");
	events->Branch("phip",&phip,"phip/D");
	events->Branch("corex",&corex,"corex[20]/D");
	events->Branch("corey",&corey,"corey[20]/D");
	events->Branch("offset_x",&offset_x,"offset_x/D");
	events->Branch("offset_y",&offset_y,"offset_y/D");
	events->Branch("xs",&xs_);
	events->Branch("ys",&ys_);
	events->Branch("pxs",&pxs_);
	events->Branch("pys",&pys_);
	events->Branch("pzs",&pzs_);
	events->Branch("thetas",&thetas_);
	events->Branch("phis",&phis_);
	events->Branch("es",&es_);
	//events->Branch("hs",&hs_);
	events->Branch("secPID",&secPID_);
	events->Branch("iter",&iter,"iter/I");

//double sigmax,sigmay,
double offsetx[10]={0},offsety[10]={0},sigmax[10]={0},sigmay[10]={0.};
double sigmax1,sigmay1;
TTree *events1 = new TTree("events1","RayTrace");
        events1->Branch("priPID",&priPID,"priPID/I");
        events1->Branch("ep",&ep,"ep/D");
        events1->Branch("thetap",&thetap,"thetap/D");
        events1->Branch("phip",&phip,"phip/D");
        events1->Branch("corex",&corex,"corex[20]/D");
        events1->Branch("corey",&corey,"corey[20]/D");
	events1->Branch("iter",&iter,"iter/I");
        /*events1->Branch("es",&es_);
        events1->Branch("xs",&xs_);
        events1->Branch("ys",&ys_);
        
	events1->Branch("offset_x",&offset_x,"offset_x/D");
        events1->Branch("offset_y",&offset_y,"offset_y/D");
	events1->Branch("offset_x1",&offset_x1,"offset_x1/D");
        events1->Branch("offset_y1",&offset_y1,"offset_y1/D");
	events1->Branch("offset_x2",&offset_x2,"offset_x2/D");
        events1->Branch("offset_y2",&offset_y2,"offset_y2/D");

        events1->Branch("sigmax",&sigmax,"sigmax/D");
        events1->Branch("sigmay",&sigmay,"sigmay/D");
        events1->Branch("sigmax1",&sigmax1,"sigmax1/D");
        events1->Branch("sigmay1",&sigmay1,"sigmay1/D");
*/
        events1->Branch("sigmax",&sigmax,"sigmax[10]/D");
        events1->Branch("sigmay",&sigmay,"sigmay[10]/D");
        events1->Branch("offsetx",&offsetx,"offsetx[10]/D");
        events1->Branch("offsety",&offsety,"offsety[10]/D");
        events1->Branch("offset_x",&offset_x,"offset_x/D");
        events1->Branch("offset_y",&offset_y,"offset_y/D");
        events1->Branch("offset_x1",&offset_x1,"offset_x1/D");
        events1->Branch("offset_y1",&offset_y1,"offset_y1/D");

	ifstream fin;
	fin.open(inputfile);
	if (!fin.is_open()) {
		cerr << "Error in opening file " << inputfile << endl;
		exit(1);
	}

	fin.seekg(0,ios::end);
	size_t filesize = fin.tellg();
	//if (jdebug>0) cout << "The file size is " << filesize << " bytes" << endl;
	if (!filesize) {
		cerr << "Warning: The file is empty!" << endl;
		fin.close();
		return(0);
	}
	fin.seekg(0,ios::beg);

	union {
		float f[_RECORD_LENGTH_];
		char  c[4*_RECORD_LENGTH_];
	} rec;


	union {
		int i;
		char c[4];
	} padding;

	char sss[5];

	int iblock = 0;
	int ievent = 0;
	int jrune = 0;
	int jrunh = 0;
	int jevte = 0;
	int jevth = 0;


	while (fin.good()) {
		//There are padding words before and after a record in the fortran output
		fin.read(padding.c,4);
		if (!fin.good()) {
			if (jrunh!=1||jrune!=1||iblock==0) {
				cerr << "Error in reading corsika file \"" << inputfile << "\"" << endl;
				cerr << "Is the file truncated?" << endl;
				if (iblock==0) cerr << "Is the file empty?" << endl;
				if (jrunh!=1) cerr << "There is no RUNH sub-block!" << endl;
				if (jrune!=1) cerr << "There is no RUNE sub-block!" << endl;
				exit(1);
			}
			//A successful reading would end here
			// if (jdebug>0)
			cerr << "Succeed in reading corsika file \"" << inputfile << "\"" << endl;
			break;
		}

		int iflag = 0;
		if (padding.i!=4*_RECORD_LENGTH_) iflag = 1;
		fin.read(rec.c,4*_RECORD_LENGTH_);
		if (!fin.good()) iflag = 1;
		fin.read(padding.c,4);
		if (padding.i!=4*_RECORD_LENGTH_) iflag = 1;

		if (iflag) {
			cerr << "Error in reading corsika file \"" << inputfile << "\"" << endl;
			cerr << "Is the file truncated?" << endl;
			exit(1);
		}

		iblock++;
int iNNN=0;
		for (int isubblock = 1; isubblock <= _NSUBBLOCK_; isubblock++) {

			int iptr = _SUBBLOCK_LENGTH_*(isubblock-1) - 1;

			memcpy(sss,&rec.f[iptr+1],4);
			sss[4] = '\0';

			if (strcmp(sss,"RUNH")==0) {
				jrunh++;
			}
			else if (strcmp(sss,"EVTH")==0) {
				thetap = rec.f[iptr+11];
				phip = rec.f[iptr+12];
				Nevent = rec.f[iptr+2];
				priPID = rec.f[iptr+3];
				ep = rec.f[iptr+4];

				for(int iuse=1; iuse<21; iuse++)
				{ 
					corex[iuse-1]=rec.f[iptr+98+iuse];
					corey[iuse-1]=rec.f[iptr+118+iuse];
				}

				sumex = 0,sumey = 0,sume = 0;
				secPID_.clear();	xs_.clear();	ys_.clear();	es_.clear();	thetas_.clear();	 phis_.clear();	 pxs_.clear();	 pys_.clear(); pzs_.clear(); 

			}

			else if (strcmp(sss,"EVTE")==0) {
				Nmax=rec.f[iptr+255+1];
				x0_sh=rec.f[iptr+255+2];
				Xmax=rec.f[iptr+255+3];

				offset_x=sumex/sume;
				offset_y=sumey/sume;
			        events->Fill();
			}

			else if (strcmp(sss,"RUNE")==0) {
				jrune++;
				ievt = rec.f[iptr+3];
			}

			else {
				for (int iptcl=1; iptcl<=39; iptcl++) {
					int iptrnow = (iptcl-1)*_PARTICLE_LENGTH_ + iptr;
#ifdef _THIN_
					cordata::weight = rec.f[iptrnow+8];
#else
					cordata::weight = 1;
#endif

#ifdef _CERENKOV_
					cordata::nclight = rec.f[iptrnow+1];
					nclight=cordata::nclight;
#else
					int idd = rec.f[iptrnow+1];
					if (idd) {
						px=rec.f[iptrnow+2]; 
						py=rec.f[iptrnow+3];
						pz=rec.f[iptrnow+4];
						pxs=px; pys=py; pzs=pz;
						xs=rec.f[iptrnow+5];   //cm
						ys=rec.f[iptrnow+6];   //cm
						secPID = int(idd/1000);
						thetas=acos(pz/sqrt(px*px+py*py+pz*pz));
						phis=atan2(py,px);
						xs2 = xs*cos(roteang) - ys*sin(roteang);
						ys2 = ys*cos(roteang) + xs*sin(roteang);
						es = GetParticleEnergy(secPID,px,py,pz);  //GeV
						es=es*1000.0;   //MeV
						if(secPID==1||secPID==2||secPID==3||secPID==5||secPID==6){
							if(es<=0) continue;
							sumex+=es*xs;
							sumey+=es*ys;
							sume+=es;
							secPID_.push_back(secPID);	xs_.push_back(xs);	ys_.push_back(ys);	es_.push_back(es);	hs_.push_back(hs);	thetas_.push_back(thetas);	phis_.push_back(phis);

						}
						hs=rec.f[iptrnow+7];
						nsec++;
					}
					else {
						cordata::weight = 0;
					}
#endif

				}
			} 
		}

	}

//	cout<<"jrunh="<<jrunh<<",jevth="<<jevth<<endl;

fin.close();
file->Write();

int Entr=events->GetEntries();
cout<<"jrunh="<<jrunh<<",jevth="<<jevth<<" ent "<<Entr<<endl; 

double meanx,meany;
vector<double> vx;
vector<double> vex;
vector<double> vy;
vector<double> vey;
double sumx,sumy;

for(int i=0;i<Entr;i++){
  events->GetEntry(i);
  for(int ii=0;ii<10;ii++){offsetx[ii]=0; offsety[ii]=0; sigmax[ii]=0;sigmay[ii]=0; }
  iter=0;
  offsetx[iter]=offset_x;
  offsety[iter]=offset_y;
  Getmean_sigam(offsetx[iter],xs_,es_,&meanx,&sigmax[iter]);
  Getmean_sigam(offsety[iter],ys_,es_,&meany,&sigmay[iter]);
  printf("i %d ep %f iter %d x %f y %f sx %f sy %f \n",i,ep,iter,offsetx[iter],offsety[iter],sigmax[iter],sigmay[iter]);
  sumx=0;sumex=0;sumey=0;sumy=0;sume=0;
  vx.clear();vy.clear();vex.clear();vey.clear();

  for(int j=0;j<xs_.size();j++){
    if(fabs((xs_.at(j)-offsetx[iter])*es_.at(j))<3*sigmax[iter]) { sumx+=xs_.at(j)*es_.at(j); sumex+=es_.at(j); vx.push_back(xs_.at(j)); vex.push_back(es_.at(j)); }
    if(fabs((ys_.at(j)-offsety[iter])*es_.at(j))<3*sigmay[iter]) { sumy+=ys_.at(j)*es_.at(j); sumey+=es_.at(j); vy.push_back(ys_.at(j)); vey.push_back(es_.at(j)); }
  }
  iter=1;
  offsetx[iter]=sumx/sumex;
  offsety[iter]=sumy/sumey;

  Getmean_sigam(offsetx[iter],vx,vex,&meanx,&sigmax[iter]);
  Getmean_sigam(offsety[iter],vy,vey,&meany,&sigmay[iter]);

  printf("i %d ep %f iter %d x %f y %f sx %f sy %f \n",i,ep,iter,offsetx[iter],offsety[iter],sigmax[iter],sigmay[iter]);

  while(iter<9&&sqrt(pow(offsetx[iter]-offsetx[iter-1],2)+pow(offsety[iter]-offsety[iter-1],2))>50)
  {
  sumx=0;sumex=0;sumey=0;sumy=0;sume=0;
  vx.clear();vy.clear();vex.clear();vey.clear();
  for(int j=0;j<xs_.size();j++){
    if(fabs((xs_.at(j)-offsetx[iter])*es_.at(j))<3*sigmax[iter]) { sumx+=xs_.at(j)*es_.at(j); sumex+=es_.at(j); vx.push_back(xs_.at(j)); vex.push_back(es_.at(j)); }
    if(fabs((ys_.at(j)-offsety[iter])*es_.at(j))<3*sigmay[iter]) { sumy+=ys_.at(j)*es_.at(j); sumey+=es_.at(j); vy.push_back(ys_.at(j)); vey.push_back(es_.at(j)); }
  }
  iter++;
  offsetx[iter]=sumx/sumex;
  offsety[iter]=sumy/sumey;
  Getmean_sigam(offsetx[iter],vx,vex,&meanx,&sigmax[iter]);
  Getmean_sigam(offsety[iter],vy,vey,&meany,&sigmay[iter]);
  printf("i %d ep %f iter %d x %f y %f sx %f sy %f \n",i,ep,iter,offsetx[iter],offsety[iter],sigmax[iter],sigmay[iter]);
  }
  offset_x=offsetx[0];offset_y=offsety[0];
  offset_x1=offsetx[iter];offset_y1=offsety[iter];

  events1->Fill();
}
file->Write();
file->Close();

}



#ifndef __CINT__
int main(int argc, char** argv) {
	/*
	   if (argc>3) {
	   minicooread_rms(argv[1],argv[2],strtol(argv[3],NULL,0));
	//minicorread(argv[1],argv[2],strtol(argv[3],NULL,0));
	}

	if (argc==4) {
	minicooread_rms(argv[1],argv[2],argv[3]);
	//minicorread(argv[1],argv[2]);
	}
	*/ 
	// else{
	if (argc==3) {
		minicooread_rms(argv[1],argv[2]);
	}
	//}
	else {
		fprintf(stderr,"Usae %s inputfile outputfile\n",argv[0]);
		exit(0);
	}

	exit(0);
}
#endif







