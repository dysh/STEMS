#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/stat.h>
#include <getopt.h>
#include <fcntl.h>
#include <time.h>
#include <regex.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define VERSION "0.1.3"
#define NUCL  "ACGT"
#define COMPL "TGCA"

const char *short_options = "hvs:l:n:d:e:o:z:";

const struct option long_options[] = {
  {"help",no_argument, NULL,'h'},
  {"version",no_argument,NULL,'v'},
  {"OTU_num",required_argument,NULL,'n'},
  {"sequence_length",required_argument,NULL,'l'},
  {"loop_mdepth",required_argument,NULL,'d'},
  {"stem_mdepth",required_argument,NULL,'s'},
  {"loop_length",required_argument,NULL,'e'},
  {"outfile",required_argument,NULL,'o'},
  {"noize",required_argument,NULL,'z'},
  {NULL,0,NULL,0}
};

int i,j,k, rez;
int option_index;
int snum, slen, msteps1,msteps2,depth,mcount,hndl;
char *lstr,tmp[11]; 
char *starter;
char **ali,**ttm;
int **tree;
char nw[2048], from[256],to[256],buf[2048],buf1[2048];
int nodeCounter,events;
int c1,c2,S;
/*S - число вариабельных сайтов*/
int scount,count1,count2,cnt,was;
double ind1,ind2,noize,pi;
FILE *fas, *phy, *nex,*struc,*tre, *stt;
char nfas[20],nphy[20],nnex[20],nstruct[20],root[16],tref[20],sttf[20];
int nflag;
char newick[2024];

const gsl_rng_type *T;
gsl_rng *r;
time_t tt;
struct tm *loctime;
