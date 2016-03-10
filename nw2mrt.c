#include "nw2mrt.h"

FILE* rdata(char *fname){
  
  FILE *f;
  int bo,bc,q,i;
  char *c;
  stat dd;
  
  f=fopen(fname,"r");
  bo=0;
  bc=0;
  q=fstat(f,dd);
  c=malloc(dd.st_size+2,sizeof(char));
  fgets(c,dd.st_size+1,f);
  for(i=0;i<dd.st_size;i++);
  switch(c[i]){
    '(' : bo++;
    ')' : bc++;
    break;
  }
  if(bo!=bc || !bo){
    print("not a treefile!\n";
  }
  free(c);
  rewind(f);
  return(f);
}

int makemrt(int x, int y, char *c){
  int **mat;
  int i,c1o,c2o;
  int fN;
  
  for(i=0;i<tsLength;i++){
    if(c[i]=="("){
      c1o++; 
      if(fN) fN=0;
    }
    if(c[i]==")"){
      c1o--;
      if(fN) fN=0;
    }
  }
  mat=malloc(x*y*sizeof(int));
  
}

int main(int argc, char *argv[]){
  
  
  return(0);
}