#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <sys/time.h>
#include "simt.h"
#define x_dim 512
#define y_dim 512
#define gf_dim 5
#define edge 2

#define uint16 unsigned short
#define uint32 unsigned int
typedef unsigned short uint16
typedef unsigned int uint32

#define serial_pixel
#define parallel_pixel



void gauss_row(uint16 **pInp,uint16 **pMid,uint16 pV, uint32 x, uint32 win_floor){

uint16 r0,r1,r2,r3,r4,r5,r6,r7,f8;
uint32 y=0;
r0=pInp[x][y];
r1=pInp[x][y+1];
r2=pInp[x][y+2];
r3=pInp[x][y+3];
r4=pInp[x][y+4];
r5=pInp[x][y+5];
r6=pInp[x][y+6];
r7=pInp[x][y+7];
r8=pInp[x][y+8];
//uint32 win_floor=(uint32)floor((double)(y_dim-2*edge)/gf_dim);
//uint32 win_ceil = (uint32)ceil((double)(y_dim-2*edge)/gf_dim);

for(uint32 win_idx=0;win_idx<win_floor;){
  pMid[x][win_idx*gf_dim+edge]=  ((r0+r4)*pV[0]+(r1+r3)*pV[1]+r2*pV[2])/17;
  pMid[x][win_idx*gf_dim+edge+1]=((r1+r5)*pV[0]+(r2+r4)*pV[1]+r3*pV[2])/17;
  pMid[x][win_idx*gf_dim+edge+2]=((r2+r6)*pV[0]+(r3+r5)*pV[1]+r4*pV[2])/17;
  pMid[x][win_idx*gf_dim+edge+3]=((r3+r7)*pV[0]+(r4+r6)*pV[1]+r5*pV[2])/17;
  pMid[x][win_idx*gf_dim+edge+4]=((r4+r8)*pV[0]+(r5+r7)*pV[1]+r6*pV[2])/17;

  r0=r5;
  r1=r6;
  r2=r7;
  r3=r8;
  win_idx++;
  y=win_idx*gf_dim+2*edge
  r4=pInp[x][y];
  r5=pInp[x][y+1];
  r6=pInp[x][y+2];
  r7=pInp[x][y+3];
  r8=pInp[x][y+4]; 
}

}

void gauss_row(uint16 **pMid,uint16 **pOut,uint16 pV, uint32 y, uint32 win_floor){

uint16 r0,r1,r2,r3,r4,r5,r6,r7,f8;
uint32 x=0;
r0=pMid[x][y];
r1=pMid[x+1][y];
r2=pMid[x+2][y];
r3=pMid[x+3][y];
r4=pMid[x+4][y];
r5=pMid[x+5][y];
r6=pMid[x+6][y];
r7=pMid[x+7][y];
r8=pMid[x+8][y];
//uint32 win_floor=(uint32)floor((double)(y_dim-2*edge)/gf_dim);
//uint32 win_ceil = (uint32)ceil((double)(y_dim-2*edge)/gf_dim);

for(uint32 win_idx=0;win_idx<win_floor;){
  pOut[win_idx*gf_dim+edge][y]=  ((r0+r4)*pV[0]+(r1+r3)*pV[1]+r2*pV[2])/17;
  pOut[win_idx*gf_dim+edge+1][y]=((r1+r5)*pV[0]+(r2+r4)*pV[1]+r3*pV[2])/17;
  pOut[win_idx*gf_dim+edge+2][y]=((r2+r6)*pV[0]+(r3+r5)*pV[1]+r4*pV[2])/17;
  pOut[win_idx*gf_dim+edge+3][y]=((r3+r7)*pV[0]+(r4+r6)*pV[1]+r5*pV[2])/17;
  pOut[win_idx*gf_dim+edge+4][y]=((r4+r8)*pV[0]+(r5+r7)*pV[1]+r6*pV[2])/17;

  r0=r5;
  r1=r6;
  r2=r7;
  r3=r8;
  win_idx++;
  x=win_idx*gf_dim+2*edge
  r4=pInp[x][y];
  r5=pInp[x+1][y];
  r6=pInp[x+2][y];
  r7=pInp[x+3][y];
  r8=pInp[x+4][y]; 
}

}
