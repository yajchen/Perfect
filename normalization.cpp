#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <complex>
#include <cmath>
#include <sys/time.h>
#include "simt.h"
#include "time.h"

#define serial
//#define parallel

//#define OUTPUTDATA

#define x_dim 512
#define y_dim 512

typedef char uint8;
typedef unsigned short uint16;
typedef unsigned int uint32;
typedef unsigned long long uint64;

void serial_stat(uint8 **pInp,uint16 *pmean,uint16 *pstd){
  uint64 Sum=0;
  uint64 SumCube=0;
  for(int x=0;x<x_dim;x++){
    for(int y=0;y<y_dim;y++){
      Sum+=pInp[x][y];
      SumCube+=(pInp[x][y]*pInp[x][y]);
    }
  }
  
  uint32 total_pixel=x_dim*y_dim;

  uint16 tmp = (uint16)(Sum/total_pixel);
  *pmean = tmp;
  *pstd  = (uint16)sqrt((SumCube/total_pixel) - (tmp*tmp));

}

void column_stat(uint8 **pInp,uint32 col,uint32 *pSumPixelPerCol,uint32 *pSumCubePixPerCol){
  uint32 Sum=0;
  uint32 SumCube=0; 
  for(int x=0;x<x_dim;x++){
     Sum+=pInp[x][col];
     SumCube+=(pInp[x][col]*pInp[x][col]);
  }
  *(pSumPixelPerCol+col)  =Sum;
  *(pSumCubePixPerCol+col)=SumCube;
}

void parallel_stat(uint8 **pInp,uint16 *pmean,uint16 *pstd){
uint32 *pSumPixelPerCol=new uint32[y_dim];
uint32 *pSumCubePixPerCol=new uint32[y_dim];

ar::simt_tau::par_for(y_dim, [&](size_t col) {
     column_stat(pInp,col,pSumPixelPerCol,pSumCubePixPerCol);
});

uint64 Sum=0;
uint64 SumCube=0;
for(int y=0;y<y_dim;y++){
  Sum     += pSumPixelPerCol[y];
  SumCube += pSumCubePixPerCol[y];
}
 
uint32 total_pixel=x_dim*y_dim;

uint16 tmp = (uint16)(Sum/total_pixel);
*pmean = tmp;
*pstd  = (uint16)sqrt((SumCube/total_pixel) - (tmp*tmp));

delete[] pSumPixelPerCol;
delete[] pSumCubePixPerCol;
}

int main(int argc, char **argv){

uint8 **pInp=new uint8*[x_dim];
uint8 **pOut=new uint8*[x_dim];
for(int x=0;x<x_dim;x++){
  *(pInp+x)=new uint8[y_dim];
  *(pOut+x)=new uint8[y_dim];
}

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

uint16 mean,std;

struct timeval tim1,tim2,tim3;
#ifdef serial
gettimeofday(&tim1, NULL); 
serial_stat(pInp,&mean,&std);
gettimeofday(&tim2, NULL); 
fprintf(stderr,"serial version\n");
#endif


#ifdef parallel
gettimeofday(&tim1, NULL);
parallel_stat(pInp,&mean,&std);
gettimeofday(&tim2, NULL);
fprintf(stderr,"parallel version\n");  
#endif


double t1=tim1.tv_sec+(tim1.tv_usec/1000000.0);
double t2=tim2.tv_sec+(tim2.tv_usec/1000000.0); 
fprintf(stderr,"@@@ %f seconds elapsed\n", t2-t1);

fprintf(stderr,"mean=%d\tstd=%d\n",mean,std);



for(int x=0;x<x_dim;x++){
  delete[] *(pInp+x);
  delete[] *(pOut+x);
}
delete[] pInp;
delete[] pOut;

return 0;
}
