#include <stdio.h>
#include <stdlib.h>
//#include <time.h>
#include <cmath>
#include <sys/time.h>
#include "simt.h"
#define x_dim 512
#define y_dim 512
#define NRUN 1
//#define serial
#define parallel

void conv_row_pixel(unsigned short **pInp,unsigned short pV[5],unsigned short **pMid,int edge,unsigned int idx){
    unsigned int x=idx/((unsigned int)y_dim);
    unsigned int y_t=idx%((unsigned int)y_dim);
    if(y_t<y_dim-2*edge){
      pMid[x][y_t+edge]=(((((unsigned short)pInp[x][edge+y_t+edge]+(unsigned short)pInp[x][edge+y_t+edge-4])*pV[0])+(((unsigned short)pInp[x][edge+y_t+edge-1]+(unsigned short)pInp[x][edge+y_t+edge-3])*pV[1])+(((unsigned short)pInp[x][edge+y_t+edge-2])*pV[2]))/17);
/*
     if(tmp<256){pMid[x][y_t+edge]=tmp;}
     else{pMid[x][y_t+edge]=255;}  
    }
    else{}
*/
    }
}
void conv_col_pixel(unsigned short **pMid,unsigned short pV[5],unsigned short **pOut,int edge,unsigned int idx){
    unsigned int y=idx/((unsigned int)x_dim);
    unsigned int x_t=idx%((unsigned int)x_dim);
    if(x_t<x_dim-2*edge){
      pOut[x_t+edge][y]=(((((unsigned short)pMid[edge+x_t+edge][y]+(unsigned short)pMid[edge+x_t+edge-4][y])*pV[0])+(((unsigned short)pMid[edge+x_t+edge-1][y]+(unsigned short)pMid[edge+x_t+edge-3][y])*pV[1])+((unsigned short)pMid[edge+x_t+edge-2][y]*pV[2]))/17);
/*
      if(tmp<256){pOut[x_t+edge][y]=(unsigned char)tmp;}
      else{pOut[x_t+edge][y]=255;}
    }
    else{}
*/
    }    
}

int main(int argc, char **argv){

int nV=5;
unsigned short pV[5];
pV[0]=1; pV[1]=4; pV[2]=7; pV[3]=4; pV[4]=1;

int edge=(nV-1)/2;

unsigned short **pInp=new unsigned short*[x_dim];
unsigned short **pMid=new unsigned short*[x_dim];
unsigned short **pOut=new unsigned short*[x_dim];

for(int i=0;i<x_dim;i++){
   *(pInp+i)=new unsigned short[y_dim];
   *(pMid+i)=new unsigned short[y_dim];
   *(pOut+i)=new unsigned short[y_dim];
}


//FILE *file;
for(int xidx=0;xidx<x_dim;xidx++){
   for(int yidx=0;yidx<y_dim;yidx++){
     pInp[xidx][yidx]=((xidx*31)+(yidx*13)+3319)%256;
     pMid[xidx][yidx]=0;
     pOut[xidx][yidx]=0;
   }
}
 fprintf(stderr,"read in image size=%d\n",x_dim);

 fprintf(stderr,"END read in image\n");



for(int nrun=0;nrun<NRUN;nrun++){

 fprintf(stderr,"nrun=%d\n",nrun);

  fprintf(stderr,"partial input image\n");
  int ix=0;
  int iy=0;
  fprintf(stderr,"Inp[%d][%d]=%d ",ix,iy,pInp[ix][iy]);

  ix+=nrun; iy+=nrun; ix=(ix+1129)%x_dim; iy=(iy+2333)%y_dim;
  fprintf(stderr,"Inp[%d][%d]=%d ",ix,iy,pInp[ix][iy]);

  ix+=nrun; iy+=nrun; ix=(ix+1129)%x_dim; iy=(iy+2333)%y_dim;
  fprintf(stderr,"Inp[%d][%d]=%d ",ix,iy,pInp[ix][iy]);

  ix+=nrun; iy+=nrun; ix=(ix+1129)%x_dim; iy=(iy+2333)%y_dim;
  fprintf(stderr,"Inp[%d][%d]=%d ",ix,iy,pInp[ix][iy]);

  ix+=nrun; iy+=nrun; ix=(ix+1129)%x_dim; iy=(iy+2333)%y_dim;
  fprintf(stderr,"Inp[%d][%d]=%d ",ix,iy,pInp[ix][iy]);

/*
 fprintf(stderr,"input image\n");
for(int xidx=0;xidx<x_dim;xidx++){
   for(int yidx=0;yidx<y_dim;yidx++){
    fprintf(stderr,"%d\t",pInp[xidx][yidx]);
   }
}
 fprintf(stderr,"\n");
*/

struct timeval tim1,tim2,tim3;
/////////along row vector///////////
#ifdef serial
gettimeofday(&tim1, NULL);  
for(unsigned int idx=0;idx<x_dim*y_dim;idx++){
   conv_row_pixel(pInp,pV,pMid,edge,idx);
   
}
/////////along col vector/////////////
for(unsigned int idx=0;idx<x_dim*y_dim;idx++){
   conv_col_pixel(pMid,pV,pOut,edge,idx);
}
gettimeofday(&tim2, NULL);  
fprintf(stderr,"serial version\n");
#endif

#ifdef parallel
    unsigned int dim=x_dim*y_dim;

    gettimeofday(&tim1, NULL);  

    ar::simt_tau::par_for(dim, [&](size_t idx) {
      conv_row_pixel(pInp,pV,pMid,edge,idx);     
    });  

    gettimeofday(&tim2, NULL);  

    ar::simt_tau::par_for(dim, [&](size_t idx) {
      conv_col_pixel(pMid,pV,pOut,edge,idx);
    });

    gettimeofday(&tim3, NULL); 
    fprintf(stderr,"parallel version\n"); 
#endif

double t1=tim1.tv_sec+(tim1.tv_usec/1000000.0);
double t2=tim2.tv_sec+(tim2.tv_usec/1000000.0);  
double t3=tim3.tv_sec+(tim3.tv_usec/1000000.0);  

fprintf(stderr,"@@@col %f seconds elapsed\n", t2-t1);
fprintf(stderr,"@@@col %f seconds elapsed\n", t3-t2);


/*
fprintf(stderr,"\nOut\n");
for(int x=0;x<x_dim;x++){
  for(int y=0;y<y_dim;y++){
    fprintf(stderr,"%d\t",pOut[x][y]);
  }
  fprintf(stderr,"\n");
}

fprintf(stderr,"\nOut\n");
*/


  ix=0;
  iy=0;
  fprintf(stderr,"Partial Out[%d][%d]=%d ",ix,iy,pOut[ix][iy]);

  ix+=nrun; iy+=nrun; ix=(ix+3631)%x_dim; iy=(iy+4297)%y_dim;
  fprintf(stderr,"Out[%d][%d]=%d ",ix,iy,pOut[ix][iy]);

  ix+=nrun; iy+=nrun; ix=(ix+3631)%x_dim; iy=(iy+4297)%y_dim;
  fprintf(stderr,"Out[%d][%d]=%d ",ix,iy,pOut[ix][iy]);

  ix+=nrun; iy+=nrun; ix=(ix+3631)%x_dim; iy=(iy+4297)%y_dim;
  fprintf(stderr,"Out[%d][%d]=%d ",ix,iy,pOut[ix][iy]);

  ix+=nrun; iy+=nrun; ix=(ix+3631)%x_dim; iy=(iy+4297)%y_dim;
  fprintf(stderr,"Out[%d][%d]=%d ",ix,iy,pOut[ix][iy]);


 fprintf(stderr,"END nrun=%d\n",nrun);
}

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
