#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <sys/time.h>
#include "simt.h"
#include "gem5.h"
#include "time.h"

#define x_dim 2048
#define y_dim 2048
#define gf_dim 5
#define edge 2
#define win_dim 7

typedef unsigned short uint16;
typedef unsigned int uint32;

#define serial
//#define parallel

void gauss_row(uint16 **pInp,uint16 **pMid,uint16 pV[gf_dim], uint32 x, uint32 win_floor){

  uint16 r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10;
  uint32 y=0;

  r7=pInp[x][y];
  r8=pInp[x][y+1];
  r9=pInp[x][y+2];
  r10=pInp[x][y+3];

  uint32 yidx;

  for(uint32 win_idx=0;win_idx<win_floor;win_idx++){
    r0=r7;
    r1=r8;
    r2=r9;
    r3=r10;

    y=win_idx*win_dim+2*edge;
    r4=pInp[x][y];
    r5=pInp[x][y+1];
    r6=pInp[x][y+2];
    r7=pInp[x][y+3];
    r8=pInp[x][y+4]; 
    r9=pInp[x][y+5]; 
    r10=pInp[x][y+6]; 
   
    yidx= win_idx*win_dim+edge;  
    pMid[x][yidx]=  ((r0+r4)*pV[0]+(r1+r3)*pV[1]+r2*pV[2])/17;
    pMid[x][yidx+1]=((r1+r5)*pV[0]+(r2+r4)*pV[1]+r3*pV[2])/17;
    pMid[x][yidx+2]=((r2+r6)*pV[0]+(r3+r5)*pV[1]+r4*pV[2])/17;
    pMid[x][yidx+3]=((r3+r7)*pV[0]+(r4+r6)*pV[1]+r5*pV[2])/17;
    pMid[x][yidx+4]=((r4+r8)*pV[0]+(r5+r7)*pV[1]+r6*pV[2])/17;
    pMid[x][yidx+5]=((r5+r9)*pV[0]+(r6+r8)*pV[1]+r7*pV[2])/17;
    pMid[x][yidx+6]=((r6+r10)*pV[0]+(r7+r9)*pV[1]+r8*pV[2])/17;

  }

  for(uint32 y=win_floor*win_dim+edge;y<y_dim-edge;y++){
    pMid[x][y]=(((pInp[x][y+2]+pInp[x][y-2])*pV[0])+((pInp[x][y+1]+pInp[x][y-1])*pV[1])+((pInp[x][y])*pV[2]))/17;
  } 
}


void gauss_col(uint16 **pMid,uint16 **pOut,uint16 pV[gf_dim], uint32 y, uint32 win_floor){

  uint16 r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10; //win_dim + gf_dim - 1 =7+5-1=11 regs
  uint32 x=0;

  r7=pMid[x][y];
  r8=pMid[x+1][y];
  r9=pMid[x+2][y];
  r10=pMid[x+3][y];

  uint32 xidx;

  for(uint32 win_idx=0;win_idx<win_floor;win_idx++){
    r0=r7;
    r1=r8;
    r2=r9;
    r3=r10;
  
    x=win_idx*win_dim+2*edge;
    r4 =pMid[x][y];
    r5 =pMid[x+1][y];
    r6 =pMid[x+2][y];
    r7 =pMid[x+3][y];
    r8 =pMid[x+4][y]; 
    r9 =pMid[x+5][y];
    r10=pMid[x+6][y];

    xidx= win_idx*win_dim+edge;
    pOut[xidx][y] =  ((r0+r4)*pV[0]+(r1+r3)*pV[1]+r2*pV[2])/17;
    pOut[xidx+1][y]= ((r1+r5)*pV[0]+(r2+r4)*pV[1]+r3*pV[2])/17;
    pOut[xidx+2][y]= ((r2+r6)*pV[0]+(r3+r5)*pV[1]+r4*pV[2])/17;
    pOut[xidx+3][y]= ((r3+r7)*pV[0]+(r4+r6)*pV[1]+r5*pV[2])/17;
    pOut[xidx+4][y]= ((r4+r8)*pV[0]+(r5+r7)*pV[1]+r6*pV[2])/17;
    pOut[xidx+5][y]= ((r5+r9)*pV[0]+(r6+r8)*pV[1]+r7*pV[2])/17;
    pOut[xidx+6][y]=((r6+r10)*pV[0]+(r7+r9)*pV[1]+r8*pV[2])/17;

  }
  for(uint32 x=win_floor*win_dim+edge;x<x_dim-edge;x++){
    pOut[x][y]=(((pMid[x+2][y]+pMid[x-2][y])*pV[0])+((pMid[x+1][y]+pMid[x-1][y])*pV[1])+((pMid[x][y])*pV[2]))/17;
  } 

}


int main(int argc, char **argv){
uint16 pV[gf_dim];
pV[0]=1; pV[1]=4; pV[2]=7; pV[3]=4; pV[4]=1;
uint16 **pInp=new uint16*[x_dim];
uint16 **pMid=new uint16*[x_dim];
uint16 **pOut=new uint16*[x_dim];

for(int i=0;i<x_dim;i++){
   *(pInp+i)=new uint16[y_dim];
   *(pMid+i)=new uint16[y_dim];
   *(pOut+i)=new uint16[y_dim];
}

//FILE *file;
for(int xidx=0;xidx<x_dim;xidx++){
   for(int yidx=0;yidx<y_dim;yidx++){
     pInp[xidx][yidx]=(uint16)(rand()%256);
     pMid[xidx][yidx]=0;
     pOut[xidx][yidx]=0;
   }
}
 fprintf(stderr,"read in image size=%dx%d\n",x_dim,y_dim);
/*
 fprintf(stderr,"input image\n");
for(int xidx=0;xidx<x_dim;xidx++){
   for(int yidx=0;yidx<y_dim;yidx++){
    fprintf(stderr,"%d\t",pInp[xidx][yidx]);
   }
}
 fprintf(stderr,"\n");
*/
/*
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
*/


uint32 ywin_floor=(uint32)floor((double)(x_dim-2*edge)/win_dim);
uint32 xwin_floor=(uint32)floor((double)(y_dim-2*edge)/win_dim);


struct timeval tim1,tim2,tim3,tim4;
gem5::DumpStats();

#ifdef serial
gem5::ResetStats();
gettimeofday(&tim1, NULL); 
/////////along row vector///////////

for(int x=0;x<x_dim;x++){
  gauss_row(pInp,pMid,pV,x,xwin_floor);
}

gettimeofday(&tim2, NULL); 
gem5::DumpStats();

/////////along col vector/////////////
gem5::ResetStats();
gettimeofday(&tim3, NULL);  
for(int y=0;y<y_dim;y++){
  gauss_col(pMid,pOut,pV,y,ywin_floor);
}
gettimeofday(&tim4, NULL);  
gem5::DumpStats();

gem5::ResetStats();
fprintf(stderr,"serial version\n");
#endif

#ifdef parallel
gem5::ResetStats();
gettimeofday(&tim1, NULL);  
/////////along row vector///////////

    ar::simt_tau::par_for(x_dim, [&](size_t x) {
      gauss_row(pInp,pMid,pV,x,xwin_floor); 
    });  

gettimeofday(&tim2, NULL);
gem5::DumpStats();

/////////along col vector/////////////
gem5::ResetStats();
gettimeofday(&tim3, NULL);   
    ar::simt_tau::par_for(y_dim, [&](size_t y) {
      gauss_col(pMid,pOut,pV,y,ywin_floor);
    });
gettimeofday(&tim4, NULL);   
gem5::DumpStats();

gem5::ResetStats();
fprintf(stderr,"parallel version\n"); 
#endif


double t1=tim1.tv_sec*1000000.0+tim1.tv_usec;
double t2=tim2.tv_sec*1000000.0+tim2.tv_usec; 
double t3=tim3.tv_sec*1000000.0+tim3.tv_usec;  
double t4=tim4.tv_sec*1000000.0+tim4.tv_usec;  
 
fprintf(stderr,"row@@@ %f seconds elapsed\n", t2-t1);
fprintf(stderr,"col@@@ %f seconds elapsed\n", t4-t3);

fprintf(stderr,"%d\n",pOut[88][97]);

/*
  ix=0;
  iy=0;
  fprintf(stderr,"Partial Out[%d][%d]=%d ",ix,iy,pOut[ix][iy]);

  ix=(ix+3631)%(x_dim-2*edge); iy=(iy+4297)%(y_dim-2*edge);
  fprintf(stderr,"Out[%d][%d]=%d ",ix,iy,pOut[ix][iy]);

  ix=(ix+3631)%(x_dim-2*edge); iy=(iy+4297)%(y_dim-2*edge);
  fprintf(stderr,"Out[%d][%d]=%d ",ix,iy,pOut[ix][iy]);

  ix=(ix+3631)%(x_dim-2*edge); iy=(iy+4297)%(y_dim-2*edge);
  fprintf(stderr,"Out[%d][%d]=%d ",ix,iy,pOut[ix][iy]);

  ix=(ix+3631)%(x_dim-2*edge); iy=(iy+4297)%(y_dim-2*edge);
  fprintf(stderr,"Out[%d][%d]=%d ",ix,iy,pOut[ix][iy]);
*/
/*
fprintf(stderr,"\nOut\n");
for(int x=edge;x<x_dim-edge;x++){
  for(int y=edge;y<y_dim-edge;y++){
    fprintf(stderr,"%d\t",pOut[x][y]);
  }
  fprintf(stderr,"\n");
}

fprintf(stderr,"\nOut\n");
*/

for(int i=0;i<x_dim;i++){
  delete[] *(pMid+i);
  delete[] *(pOut+i);
  delete[] *(pInp+i);
}
delete[] pInp;
delete[] pMid;
delete[] pOut;

return 0;

}
