////////////////////////restrictions//////////////////
// V: nV-by-1 vector   Gaussian G=V*V': nV-nV matrix//
// V is strictly symmetric                          //
//Sliding window on row vector for the entire image //
//Then                                              //
//Sliding window on col vector for the entire image //
//////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <complex>
#include <iostream>
#include <string>
#include <fstream>
#include <stdexcept>
#include <complex>
#include "gauss.h"
using namespace std;

int main(){

int nV=5;
int *pV=new int[5];
pV[0]=1; pV[1]=4; pV[2]=7; pV[3]=4; pV[4]=1;

int MA=20;
int NA=25;
int seed=-11;
int sigma=10;
cout<<"A is"<<endl;
int **pA=new int*[MA];
for(int ma=0;ma<MA;ma++){
  *(pA+ma)=new int[NA];
  for(int na=0;na<NA;na++){
    pA[ma][na]=(int)(sigma*(gauss1(&seed)));
    cout<<pA[ma][na]<<"\t";
  }
  cout<<endl;
}
cout<<endl;


int edge=(nV-1)/2;

int MC=MA;
int NC=NA;
int **pC=new int*[MC];
int **pD=new int*[MC];
for(int mc=0;mc<MC;mc++){
  *(pC+mc)=new int[NC];
  *(pD+mc)=new int[NC];
}


/////////along row vector///////////
for(int mc=0;mc<MC;mc++){
  for(int nc=edge;nc<NA-edge;nc++){
//     pC[mc][nc]=0;
//     for(int j=0;j<nV;j++)
//       pC[mc][nc]+=((pA[mc][edge+nc-j])*(pV[j]));
     pC[mc][nc]=((pA[mc][edge+nc-0]+pA[mc][edge+nc-4])*(pV[0])) + ((pA[mc][edge+nc-1]+pA[mc][edge+nc-3])*(pV[1])) + ((pA[mc][edge+nc-2])*(pV[2]));
  
  }
}
cout<<"C is"<<endl;
for(int mc=0;mc<MC;mc++){
  for(int nc=0;nc<NC;nc++){
     cout<<pC[mc][nc]<<"\t";
  }
  cout<<endl;
}


/////////along col vector/////////////
for(int nc=0;nc<NA;nc++){
  for(int mc=edge;mc<MC-edge;mc++){
//    pD[mc][nc]=0;    
//    for(int j=0;j<nV;j++)
//       pD[mc][nc]+=((pC[edge+mc-j][nc])*(pV[j]));
    pD[mc][nc]=((pC[edge+mc-0][nc]+pC[edge+mc-4][nc])*(pV[0])) + ((pC[edge+mc-1][nc]+pC[edge+mc-3][nc])*(pV[1]))  + ((pC[edge+mc-2][nc])*(pV[2]));
  }
}

cout<<"D is"<<endl;
for(int mc=0;mc<MC;mc++){
  for(int nc=0;nc<NC;nc++){
     cout<<pD[mc][nc]<<"\t";
  }
  cout<<endl;
}

delete[] pV;
for(int ma=0;ma<MA;ma++){
  delete[] *(pA+ma);
  delete[] *(pC+ma);
  delete[] *(pD+ma);
}
delete[] pA;
delete[] pC;
delete[] pD;

}
