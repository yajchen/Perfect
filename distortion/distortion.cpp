#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <complex>
#include <cmath>
#include <sys/time.h>
#include "simt.h"
#include "gem5.h"
#include "time.h"

//#define serial
#define parallel

//#define OUTPUTDATA

#define x_dim 640
#define y_dim 480

#define stride_enable
#define stride 2

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
float a,b;
uint16 tmp;

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

      a = (float)(y-y1);
      b = (float)(x-x1);
      
      tmp = (uint16)((1-a)*(1-b)*((float)pInp[x1][y1]) + (1-a)*b*((float)pInp[x2][y1]) + a*(1-b)*((float)pInp[x1][y2]) + a*b*((float)pInp[x2][y2]));
      
      pOut[x][y] = (uint8)(tmp<=255 ? tmp:255);

    }
  }

}

void parallel_calc_table(size_t idx, float *pXtable,float *pYtable){
size_t x = idx/y_dim;
size_t y = idx%y_dim;

float xd = ((float)(x-Xc))/Focal;      
float yd = ((float)(y-Yc))/Focal;
float r2 = xd*xd + yd*yd;
float r4 = r2*r2;
float factor = 1.0f + k1*r2 + k2*r4;

pXtable[x*y_dim+y] = xd*factor*Focal + Xc;
pYtable[x*y_dim+y] = yd*factor*Focal + Yc;
}


void column_bilinear_interp(uint8 **pInp,uint8 **pOut,size_t col, float *pXtable, float *pYtable){
for(int x=0;x<x_dim;x++){
  int x1 = (int)floor((pXtable[x*y_dim+col]));           
      x1 = max(0,x1);
  int x2 = (int)ceil((pXtable[x*y_dim+col]));
      x2 = min(x2,(x_dim-1));

  int y1 = (int)floor((pYtable[x*y_dim+col]));
      y1 = max(0,y1);
  int y2 = (int)ceil((pYtable[x*y_dim+col]));
      y2 = min(y2,(y_dim-1));

  float a = (float)(col-y1);
  float b = (float)(x-x1);
  uint16 tmp = (uint16)((1-a)*(1-b)*((float)pInp[x1][y1]) + (1-a)*b*((float)pInp[x2][y1]) + a*(1-b)*((float)pInp[x1][y2]) + a*b*((float)pInp[x2][y2]));
      
  pOut[x][col] = (uint8)(tmp<=255 ? tmp:255);

}

}

#ifdef stride_enable
void column_bilinear_interp_stride(uint8 **pInp,uint8 **pOut,size_t tid, float *pXtable, float *pYtable){
for(int x=0;x<x_dim;x++){
  for(int s=0;s<stride;s++){
    size_t col = tid*stride+s;
    int x1 = (int)floor((pXtable[x*y_dim+col]));           
        x1 = max(0,x1);
    int x2 = (int)ceil((pXtable[x*y_dim+col]));
        x2 = min(x2,(x_dim-1));

    int y1 = (int)floor((pYtable[x*y_dim+col]));
        y1 = max(0,y1);
    int y2 = (int)ceil((pYtable[x*y_dim+col]));
        y2 = min(y2,(y_dim-1));

    float a = (float)(col-y1);
    float b = (float)(x-x1);
    uint16 tmp = (uint16)((1-a)*(1-b)*((float)pInp[x1][y1]) + (1-a)*b*((float)pInp[x2][y1]) + a*(1-b)*((float)pInp[x1][y2]) + a*b*((float)pInp[x2][y2]));
      
    pOut[x][col] = (uint8)(tmp<=255 ? tmp:255);

  }

}

}
#endif

int main(int argc, char **argv){

uint8 **pInp=new uint8*[x_dim];
uint8 **pOut=new uint8*[x_dim];
for(int x=0;x<x_dim;x++){
  *(pInp+x)=new uint8[y_dim];
  *(pOut+x)=new uint8[y_dim];
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

struct timeval tim1,tim2,tim3,tim4;
#ifdef serial
gem5::ResetStats();
//above is initial stat.txt

//start collect stat.txt for distortion-build-table
gem5::DumpStats();
gettimeofday(&tim1, NULL); 
serial_calc_table(pXtable,pYtable);
gettimeofday(&tim2, NULL); 
gem5::ResetStats();
//end collect stat.txt for distortion-build-table

//start collect stat.txt for distortion-correction
gem5::DumpStats();
gettimeofday(&tim3, NULL); 
serial_bilinear_interp(pInp,pOut,pXtable,pYtable);
gettimeofday(&tim4, NULL);
gem5::ResetStats();
//end collect stat.txt for distortion-correction

gem5::DumpStats();
//below is tail stat.txt 
fprintf(stderr,"serial version\n");
#endif


#ifdef parallel
gem5::ResetStats();
//above is initial stat.txt

//start collecting stat.txt for parallel distortion-build-table
gem5::DumpStats();
gettimeofday(&tim1, NULL); 
size_t total_pixel = x_dim*y_dim;
ar::simt_tau::par_for(total_pixel, [&](size_t idx) {
    parallel_calc_table(idx,pXtable,pYtable);
});
gettimeofday(&tim2, NULL); 
gem5::ResetStats();
//end collecting stat.txt for parallel distortion-build-table

//start collecting stat.txt for parallel distortion-correction
gem5::DumpStats();
gettimeofday(&tim3, NULL); 
//parallel_bilinear_interp(pInp,pOut,pXtable,pYtable);
#ifdef stride_enable
size_t num_threads=y_dim/stride;
ar::simt_tau::par_for(num_threads, [&](size_t tid) {
    column_bilinear_interp_stride(pInp,pOut,tid,pXtable,pYtable);
});
#else
ar::simt_tau::par_for(y_dim, [&](size_t col) {
    column_bilinear_interp(pInp,pOut,col,pXtable,pYtable);
});
#endif
gettimeofday(&tim4, NULL);
gem5::ResetStats();
//end collecting stat.txt for parallel distortion-correction

gem5::DumpStats();
//below is tail stat.txt 
#ifdef stride_enable
fprintf(stderr,"parallel-stride-2 version\n");
#else
fprintf(stderr,"parallel version\n");
#endif

#endif

double t1=tim1.tv_sec*1000000.0+tim1.tv_usec;
double t2=tim2.tv_sec*1000000.0+tim2.tv_usec; 
double t3=tim3.tv_sec*1000000.0+tim3.tv_usec;  
double t4=tim4.tv_sec*1000000.0+tim4.tv_usec;  
fprintf(stderr,"calculate_table@@@ %f microseconds elapsed\n", t2-t1);
fprintf(stderr,"bilinear_interp@@@ %f microseconds elapsed\n", t4-t3);

//fprintf(stderr,"%d\n",pOut[88][97]);

#ifdef OUTPUTDATA
fprintf(stderr,"output image\n");
for(int xidx=0;xidx<x_dim;xidx++){
   for(int yidx=0;yidx<y_dim;yidx++){
    fprintf(stderr,"%d\t",pOut[xidx][yidx]);
   }
}
fprintf(stderr,"\n");
#endif

for(int x=0;x<x_dim;x++){
  delete[] *(pInp+x);
  delete[] *(pOut+x);
}
delete[] pInp;
delete[] pOut;

delete[] pXtable;
delete[] pYtable;

return 0;

}
