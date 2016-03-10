#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/stat.h>
#include <getopt.h>
#include <fcntl.h>
#include <time.h>
#include <regex.h>

#define VERSION "0.9.0"
#define PATTERN1 "(\[\{\<)"
#define NUCL "ACGT"

typedef struct {
  char nme[100];
  char *seq;
} otu;

typedef struct {
  int x;
  int y;
} point;

typedef struct {
  float x;
  float y;
} fpoint;

typedef struct {
  int snum;
  int slen;
  otu *ali;
  int *brackets;
  int *ranks;
  int *loop;
} data;

struct sta {
  int loops;
  int stems;
  float stemsR;
  float loopsR;
};

struct {
  int start;
  int end;
  int num;
  float len;
  float ts[2];
  float tv[2];
} esta;

typedef struct {
    int start;
    int end;
    float prs[15];
}dbl;

data *dta, *ndta;
dbl dblts;
int i,j,k,l,ofsize,bo,bc,nd,flag,count,nspec,slen,nloops;
int *loo;
int *ranks,r;
fpoint *pairs;
int loopflag,lloop, npairs;
char c,c1,c2;
char line[1024],name[100],*seq,*ifile,*cfile,*sfile,t[100],sq[1000];
char tfname[100];

/*for statistics*/
int l1[6],s1[20];
double l2[6],s2[20];
  
FILE *in, *in1, *out, *tmp, *rep, *svg;
struct stat fbuf;
struct info {
  char alignment[100];
  char structure[100];
  char outfile[100];
} info;

int iflag, nref;
int lflag, minloop;
int uflag;
int eflag;
int sflag;
int rflag; /*RNA?*/
int xflag; /*это на случай если надо оставить только уникальные или достаточно далекие гаплотипы*/

const char *short_options = "hvta:s:o:x:i:l:e:u";
const struct option long_options[] = {
  {"help",no_argument, NULL,'h'},
  {"version",no_argument,NULL,'v'},
  {"clustal_file",required_argument,NULL,'a'},
  {"structure_file",required_argument,NULL,'s'},
  {"output_file",required_argument,NULL,'o'},
  {"single_structure",required_argument,NULL,'i'},
  {"separate_loops",required_argument,NULL,'l'},
  {"extract_fasta",required_argument,NULL,'e'},
  {"otu_list",no_argument,NULL,'u'},
  {"statistics",no_argument, NULL,'t'},
  {"uninque_OTUs",required_argument,NULL,'x'},
  {NULL,0,NULL,0}
};

int rez;
int option_index;
time_t lt;
  
  
  
