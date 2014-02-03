#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <complex>
#include <cmath>
#include <sys/time.h>
#include "simt.h"
#include "time.h"

//#define serial
#define parallel

//#define OUTPUTDATA

#define x_dim 640
#define y_dim 480

#define k1 -0.303621f	
#define k2 0.0939467f

#define Xc 304.514f
#define Yc 264.369f

#define Focal 452.16f

typedef char uint8;
typedef unsigned short uint16;
typedef unsigned int uint32;
typedef unsigned long long uint64;

int min(int a,int b){return a<b ? a:b;}
int max(int a,int b){return a>b ? a:b;}

void serial_calc_table(float *pXtable,float *pYtable){
  
  float xd,yd,r2,r4,factor;

  for(int x=0;x<x_dim;x++){
    for(int y=0;y<y_dim;y++){
      xd = ((float)(x-Xc))/Focal;      
      yd = ((float)(y-Yc))/Focal;
      r2 = xd*xd + yd*yd;
      r4 = r2*r2;
      factor = 1.0f + k1*r2 + k2*r4;
      pXtable[x*y_dim+y] = xd*factor*Focal + Xc;
      pYtable[x*y_dim+y] = yd*factor*Focal + Yc;
    }
  }

}

void serial_bilinear_interp(uint8 **pInp,uint8 **pOut,float *pXtable,float *pYtable){

int x1,x2,y1,y2;

  for(int x=0;x<x_dim;x++){
    for(int y=0;y<y_dim;y++){
      x1 = (int)floor((pXtable[x*y_dim+y]));           
      x1 = max(0,x1);
      x2 = (int)ceil((pXtable[x*y_dim+y]));
      x2 = min(x2,(x_dim-1));

      y1 = (int)floor((pYtable[x*y_dim+y]));
      y1 = max(0,y1);
      y2 = (int)ceil((pYtable[x*y_dim+y]));
      y2 = min(y2,(y_dim-1));


    }
  }

}

int main(int argc, char **argv){

uint8 **pInp=new uint8*[x_dim];
//uint8 **pOut=new uint8*[x_dim];
for(int x=0;x<x_dim;x++){
  *(pInp+x)=new uint8[y_dim];
//  *(pOut+x)=new uint8[y_dim];
}

float *pXtable=new float[x_dim*y_dim];
float *pYtable=new float[x_dim*y_dim];

srand(time(NULL));
/////////generate data/////////
for(int x=0;x<x_dim;x++){
  for(int y=0;y<y_dim;y++){
    pInp[x][y] = (uint16)(rand()%256);
  }
}
fprintf(stderr,"input image size %d x %d\n",x_dim,y_dim);

#ifdef OUTPUTDATA
fprintf(stderr,"input image\n");
for(int xidx=0;xidx<x_dim;xidx++){
   for(int yidx=0;yidx<y_dim;yidx++){
    fprintf(stderr,"%d\t",pInp[xidx][yidx]);
   }
}
fprintf(stderr,"\n");
#endif

struct timeval tim1,tim2;
#ifdef serial
gettimeofday(&tim1, NULL); 
serial_calc_table(pXtable,pYtable);
gettimeofday(&tim2, NULL); 
fprintf(stderr,"serial version\n");
#endif

/*
#ifdef OUTPUTDATA
fprintf(stderr,"output image\n");
for(int xidx=0;xidx<x_dim;xidx++){
   for(int yidx=0;yidx<y_dim;yidx++){
    fprintf(stderr,"%d\t",pOut[xidx][yidx]);
   }
}
fprintf(stderr,"\n");
#endif
*/
for(int x=0;x<x_dim;x++){
  delete[] *(pInp+x);
//  delete[] *(pOut+x);
}
delete[] pInp;
//delete[] pOut;

delete[] pXtable;
delete[] pYtable;

return 0;

}
