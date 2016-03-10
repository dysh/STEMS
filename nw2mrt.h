#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/stat.h>

FILE *fin, *fout;

typedef struct {
  int notu;
  int nranks;
  int **m;
} mrt;

m *mrt;
int co,cc,i,j,tsLength;
