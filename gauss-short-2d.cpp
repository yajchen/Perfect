#include <stdio.h>
#include <stdlib.h>
//#include <time.h>
#include <cmath>
#include <sys/time.h>
#include "simt.h"
#define x_dim 512
#define y_dim 512
#define fir_val 6 
#define edge 2
#define NRUN 1
//#define serial
#define parallel

typedef unsigned short u16;
typedef unsigned int u32;

void conv2d_pixel(u16 **pInp, u16 **pOut, u16 pV[fir_val], u32 idx, u32 out_y_dim){
  
  u32 v0,v1,v2,v3,v4,v5;

  u32 x=idx/out_y_dim;
  u32 y=idx%out_y_dim;

  v0=pInp[x+2][y+2];
  v1=pInp[x+1][y+2]+pInp[x+2][y+1]+pInp[x+2][y+3]+pInp[x+3][y+2];
  v2=pInp[x+1][y+1]+pInp[x+1][y+3]+pInp[x+3][y+1]+pInp[x+3][y+3];
  v3=pInp[x][y+2]+pInp[x+2][y+0]+pInp[x+2][y+4]+pInp[x+4][y+2];
  v4=pInp[x][y+1]+pInp[x][y+3]+pInp[x+1][y]+pInp[x+1][y+4]+pInp[x+3][y]+pInp[x+3][y+4]+pInp[x+4][y+1]+pInp[x+4][y+3];
  v5=pInp[x][y]+pInp[x][y+4]+pInp[x+4][y]+pInp[x+4][y+4];

  pOut[x][y] = (u16)((v0*pV[0]+v1*pV[1]+v2*pV[2]+v3*pV[3]+v4*pV[4]+v5*pV[5])/289); 
}

int main(int argc, char **argv){


unsigned short pV[fir_val];
pV[0]=49; pV[1]=28; pV[2]=16; pV[3]=7; pV[4]=4; pV[5]=1;


unsigned short **pInp=new unsigned short*[x_dim];

for(int i=0;i<x_dim;i++)
   *(pInp+i)=new unsigned short[y_dim];


u32 out_x_dim = x_dim-2*edge;
u32 out_y_dim = y_dim-2*edge;

unsigned short **pOut=new unsigned short*[out_x_dim];
for(unsigned int i=0;i<out_x_dim;i++)
   *(pOut+i)=new unsigned short[out_y_dim];

for(unsigned int xidx=0;xidx<x_dim;xidx++)
   for(unsigned int yidx=0;yidx<y_dim;yidx++)
     pInp[xidx][yidx]=((xidx*31)+(yidx*13)+3319)%256;

 fprintf(stderr,"read in image size=%dx%d\n",x_dim,y_dim);


  fprintf(stderr,"partial input image\n");
  int ix=0;
  int iy=0;
  fprintf(stderr,"Inp[%d][%d]=%d ",ix,iy,pInp[ix][iy]);

  ix=(ix+1129)%x_dim; iy=(iy+2333)%y_dim;
  fprintf(stderr,"Inp[%d][%d]=%d ",ix,iy,pInp[ix][iy]);

  ix=(ix+1129)%x_dim; iy=(iy+2333)%y_dim;
  fprintf(stderr,"Inp[%d][%d]=%d ",ix,iy,pInp[ix][iy]);

  ix=(ix+1129)%x_dim; iy=(iy+2333)%y_dim;
  fprintf(stderr,"Inp[%d][%d]=%d ",ix,iy,pInp[ix][iy]);

  ix=(ix+1129)%x_dim; iy=(iy+2333)%y_dim;
  fprintf(stderr,"Inp[%d][%d]=%d ",ix,iy,pInp[ix][iy]);


/*
 fprintf(stderr,"input image\n");
for(unsigned int xidx=0;xidx<x_dim;xidx++){
   for(unsigned int yidx=0;yidx<y_dim;yidx++){
    fprintf(stderr,"%d\t",pInp[xidx][yidx]);
   }
}
*/
struct timeval tim1,tim2;
#ifdef serial
gettimeofday(&tim1, NULL);  
  
for(unsigned int idx=0;idx<out_x_dim*out_y_dim;idx++)
   conv2d_pixel(pInp,pOut,pV,idx,out_y_dim);

gettimeofday(&tim2, NULL);  
fprintf(stderr,"serial version\n");
#endif

#ifdef parallel
u32 out_dim=out_x_dim*out_y_dim;
gettimeofday(&tim1, NULL);  
    ar::simt_tau::par_for(out_dim, [&](size_t idx) {
      conv2d_pixel(pInp,pOut,pV,idx,out_y_dim);
    });
gettimeofday(&tim2, NULL);  
fprintf(stderr,"parallel version\n"); 
#endif
double t1=tim1.tv_sec+(tim1.tv_usec/1000000.0);
double t2=tim2.tv_sec+(tim2.tv_usec/1000000.0);  

fprintf(stderr,"@@@2dconv %f seconds elapsed\n", t2-t1);

/*
fprintf(stderr,"\nOut\n");
for(unsigned int x=0;x<out_x_dim;x++){
  for(unsigned int y=0;y<out_y_dim;y++){
    fprintf(stderr,"%d\t",pOut[x][y]);
  }
  fprintf(stderr,"\n");
}

*/
fprintf(stderr,"\nOut\n");



  ix=0;
  iy=0;
  fprintf(stderr,"Partial Out[%d][%d]=%d ",ix,iy,pOut[ix][iy]);

  ix=(ix+3631)%out_x_dim; iy=(iy+4297)%out_y_dim;
  fprintf(stderr,"Out[%d][%d]=%d ",ix,iy,pOut[ix][iy]);

  ix=(ix+3631)%out_x_dim; iy=(iy+4297)%out_y_dim;
  fprintf(stderr,"Out[%d][%d]=%d ",ix,iy,pOut[ix][iy]);

  ix=(ix+3631)%out_x_dim; iy=(iy+4297)%out_y_dim;
  fprintf(stderr,"Out[%d][%d]=%d ",ix,iy,pOut[ix][iy]);

  ix=(ix+3631)%out_x_dim; iy=(iy+4297)%out_y_dim;
  fprintf(stderr,"Out[%d][%d]=%d ",ix,iy,pOut[ix][iy]);


for(unsigned int i=0;i<x_dim;i++)
  delete[] *(pInp+i);

delete[] pInp;

for(unsigned int i=0;i<out_x_dim;i++)
  delete[] *(pOut+i);

delete[] pOut;

return 0;
}
