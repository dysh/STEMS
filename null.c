#include "null.h"

void randsetup(){
	tt=time(0);
	gsl_rng_env_setup();
	T=gsl_rng_default;
	r=gsl_rng_alloc(T);
	gsl_rng_set(r,(long)tt);
	
	return;
}

int findSubstr(char *d, char *t){
  int res,l1,l2,l,i,k;
  
  res = 0;
  j=0;
  l1=strlen(d);
  l2=strlen(t);
  
  if(l1<l2) return(-1);
  l=l1-l2;
  for(i=0;i<l;i++){
    j=0;
    for(k=0;k<l2;k++){
      if(d[i+k]!=t[k]){
	j=1;
	break;
      }
    }
    if(j==0) break;
  }
  res=i;
  //printf("%s\t%s\tres=%i\n",d,t,res);
  return(res);  
}
  

int cpseq(int trg, int sr){
  int c;
  for(c=0;c<slen;c++) ali[trg][c]=ali[sr][c];
  return(0);
}

int divpos(char **a, int n, int nseq, int nmin, int nmax){
  int i,j,res;
  
  res=0;
  if(nmax>n) nmax=n;
  for(i=nmin;i<nmax;i++){
    for(j=1;j<nseq;j++){
      if(a[j][i]!=a[0][i]){
	res +=1;
	break;
      }
    }
  }
  return(res);
}

double piCalc(char **a, int n, int lmin, int lmax){
  double res;
  int i,j,k;
  double sum1,sum2;
  
  sum2=0.0;
  for(i=1;i<n;i++){
    for(j=0; j<i;j++){
      sum1=0.0;
      for(k=lmin;k<lmax;k++){
	if(a[i][k]!=a[j][k]) sum1 +=1.0;
      }
      sum1 = sum1/(double)(lmax-lmin);
    }
    sum2 += sum1;
  }
  res=sum2*2.0;
  
  return(res);
}

double tajCalc(double S,double pi, double n){
    
    int j;
    double a1,a2,b1,b2,c1,c2,e1,e2,D;
    
    a1 =0;
    a2 =0;
    for(j=1;j<n;j++){
	a1 += 1.0/(double)j;
	a2 += 1.0/((double)j*(double)j);
    }
    b1=(n+1.0)/(n-1.0)/3.0;
    b2=(n*n+n+3.0)*2.0/(double)(9*n*(n-1));
    c1=b1-1/a1;
    c2=b2-(n+2)/(a1*n)+a2/a1/a1;
    e1=c1/a1;
    e2=c2/(a1*a1+a2);

    D=(pi-S/a1)/(sqrt(e1*S+e2*S*(S-1.0)));
    return(D);
}

int main(int argc, char *argv[]){
  nflag=0;
  while ((rez=getopt_long(argc,argv,short_options,long_options,&option_index))!=-1){
    switch(rez){
      case 'h':{
	printf("\n\tBRIEF HELP\n======================================================================\n\n");
	printf("\t<l>\t(obligatory)\tLength of a sequence (INT)\n");
	printf("\t<n>\t(obligatory)\tnumber of sequences (INT)\n");
	printf("\t<e>\t(obligatory)\tLength of a loop (e<l/2) (INT)\n");
	printf("\t<d>\t(obligatory)\tnumber of mutational steps in loop (unpaired) region (INT)\n");
	printf("\t<s>\t(obligatory)\tnumber of mutational steps in stem (paired) region (INT)\n");
	printf("\t<z>\t(optionally)\tratio of improper mutations in stem (paired) region (DOUBLE<1.0)\n");
	printf("\t<o>\t(optionally)\tcommon root for all output filenames (default:\"outfile\" (STRING)\n");
	
	
	return(0);
	break;
      }
      case 'v':{
	printf("Version %s\n", VERSION);
	return(0);
	break;
      }
      case 'n':{
	if(optarg!=NULL) snum=atoi(optarg);
	else{
	  printf("\nno argument to \"n\"!\n Quitting\n");
	  return(1);
	}
	if(!snum){
	  printf("\nnumber of sequences = 0?\nQuitting\n");
	  return(2);
	}
	break;
      }
      case 'l':{
	if(optarg!=NULL) slen=atoi(optarg);
	else{
	  printf("\nno argument to \"l\"!\n Quitting\n");
	  return(1);
	}
	if(!slen){
	  printf("\nlength of a sequence = 0?\nQuitting\n");
	  return(2);
	}
	break;
      }
      case 'd':{
	if(optarg!=NULL) msteps1=atoi(optarg);
	else{
	  printf("\nno argument to \"d\"!\n Quitting\n");
	  return(1);
	}
	if(!msteps1){
	  printf("\ntreelength for loops = 0?\nQuitting\n");
	  return(2);
	}
	break;
      }
      case 's':{
	if(optarg!=NULL) msteps2=atoi(optarg);
	else{
	  printf("\nno argument to \"s\"!\n Quitting\n");
	  return(1);
	}
	if(!msteps2){
	  printf("\ntreelength for stems = 0?\nQuitting\n");
	  return(2);
	}
	break;
      }
      case 'e':{
	if(optarg!=NULL) depth=atoi(optarg);
	else{
	  printf("\nno argument to \"e\"!\n Quitting\n");
	  return(1);
	}
	if(!depth){
	  printf("\nlength of the loop = 0?\nQuitting\n");
	  return(2);
	}
	break;
      }
      case 'z':{
	if(optarg!=NULL) noize=atof(optarg);
	else{
	  printf("\nno argument to \"z\"!\n Quitting\n");
	  return(1);
	}
	if(!noize){
	  printf("\nnoize to doublet model = 0.0!\n");
	}
	break;
      }
      case 'o':{
	if(strlen(optarg)>16){
	  printf("filename \"%s\" is way too long. Switching to deault \"outfile\"\n",optarg);
	  sprintf(nfas,"outfile.fas");
	  sprintf(nphy,"outfile.phy");
	  sprintf(nnex,"outfile.nex");
	  sprintf(nstruct,"outfile.STR");
	  sprintf(tref,"outfile.tre");
	  sprintf(root,"outfile");
	  sprintf(sttf,"outfile.dat");
	} else {
	  sprintf(nfas,"%s.fas",optarg);
	  sprintf(nphy,"%s.phy",optarg);
	  sprintf(nnex,"%s.nex",optarg);
	  sprintf(nstruct,"%s.STR",optarg);
	  sprintf(root,"%s",optarg);
	  sprintf(tref,"%s.tre", optarg);
	  sprintf(sttf,"%s.dat",optarg);
	}
	fas=fopen(nfas,"w");
	phy=fopen(nphy,"w");
	nex=fopen(nnex,"w");
	struc=fopen(nstruct,"w");
	tre=fopen(tref,"w");
	stt=fopen(sttf,"w");
	nflag=1;
	break;
      }
      default: break;
    }
  }
  if(snum<3 || slen<10 || msteps1<1 || msteps2<1 || depth<1){
    printf("wrong argument(s)!\nsnum=\t%i\nslen=\t%i\nmsteps1=\t%i\nmsteps2=\t%i\ndepth=\t%i\n",snum,slen,msteps1,msteps2,depth);
    return(-1);
  }
  if(noize>1) noize=0.0;
		 
  lstr=malloc(slen*sizeof(int));  
  starter=malloc(slen*sizeof(int));
  if(!nflag){
     sprintf(nfas,"outfile.fas");
     sprintf(nphy,"outfile.phy");
     sprintf(nnex,"outfile.nex");
     sprintf(nstruct,"outfile.STR");
     sprintf(root,"outfile");
     sprintf(tref,"outfile.tre");
     sprintf(sttf,"outfile.dat");
     fas=fopen(nfas,"w");
     phy=fopen(nphy,"w");
     nex=fopen(nnex,"w");
     struc=fopen(nstruct,"w");
     tre=fopen(tref,"w");
     stt=fopen(sttf,"w");
  }
  randsetup();
  for(i=0;i<depth;i++){
    lstr[i]='(';
    lstr[i+depth]=')';
    j=gsl_rng_uniform_int(r,4);
    starter[i]=NUCL[j];
    
    k=2*depth-i-1;
    starter[k]=COMPL[j];
    
  }
  for(i=2*depth;i<slen;i++){
    lstr[i]='.';
    j=gsl_rng_uniform_int(r,4);
    starter[i]=NUCL[j];
  }
  fprintf(struc,">Starting_sequence\n%s\n%s\n",starter,lstr);
  
  ali=malloc(snum*sizeof(char *));
  ttm=malloc(snum*sizeof(char *));
  tree=malloc(snum*sizeof(int *));
  for(i=0;i<snum;i++){
    ali[i]=malloc(slen*sizeof(char));
    ttm[i]=malloc(slen*sizeof(char));
    tree[i]=malloc(snum*sizeof(int));
  }
  for(i=0;i<slen;i++) ali[0][i]=starter[i];
  cpseq(1,0);
  sprintf(nw,"(item0,item1);");
  tree[0][1]=1;
  scount=2;
  ind1=(double)(snum)/(double)(snum+msteps1+msteps2);
//  ind2=(double)msteps1/(double)(msteps1+msteps2);
  ind2=2.0/(double)(msteps1+msteps2);
  count1=0;
  count2=0;
  mcount=msteps1+msteps2;
  events=snum+msteps1+msteps2;
  while(events>0){
    ind2=(double)scount/(double)(msteps1+msteps2);
    if(gsl_rng_uniform(r)<ind1){
      if(scount<snum){
	j=gsl_rng_uniform_int(r,scount);
	cpseq(scount,j);
	tree[j][scount]=scount;
	tree[scount][scount]=j;
	sprintf(from,"item%i",j);
	sprintf(to,"(%s,item%i)",from,scount);
	k=findSubstr(nw,from);
      //sprintf(buf,"%c",'\0');
	buf[0]='\0';
	strncat(buf,nw,k);
	strcat(buf,to);
	k+=strlen(from);
	strcat(buf,nw+k);      
	sprintf(nw,"%s",buf);
	scount++;
	printf("species: %i\n",scount);
      }
    } 
    if(gsl_rng_uniform(r)>ind1){ 
      if(mcount>0 && scount<=snum){
//	if(gsl_rng_uniform(r)>ind2 && count2<=msteps2){
	  k=gsl_rng_uniform_int(r,scount);
	  j=gsl_rng_uniform_int(r,slen);
	  was=j;
	  c2=ali[k][j];
	  c1=ali[k][j];
	  while(c1==c2){ 
	    cnt=gsl_rng_uniform_int(r,4);
	    c2=NUCL[cnt];
	  }
	  if(j<depth && count2<msteps2){
	    ali[k][j]=c2;
	  // printf("mutation in stem pos:\t%i\tof\t%i\n",j,count2);
	  /*next line makes doublet a bit noizy*/
	    if(gsl_rng_uniform(r)>noize) ali[k][2*depth-j-1]=COMPL[cnt];
	    count2++;
	    mcount -=1;
	  } else {
	    if(count1<=msteps1 && j>=depth*2){
	   
	    c2=ali[k][j];
	    c1=ali[k][j];
	    while(c1==c2){ 
	      cnt=gsl_rng_uniform_int(r,4);
	      c2=NUCL[cnt];
	    }
	    ali[k][j]=c2;
	    count1++;
	    mcount -=1;
	    //printf("mutation in loop pos:\t%i\tof\t%i\n",j,count1);
	  }
	}
      }
    } 
    
    events=snum-scount+mcount;
 //   printf("%i\t%i\t%i\t%i\t\n",events,snum,scount,mcount);
  }
  fprintf(tre,"%s\n",nw);
  /*-------------------------*/
  fprintf(nex,"#NEXUS\n\n");
  fprintf(nex,"[created with NULL simulator, where a starting sequence evolved along random tree\n");
  fprintf(nex,"so that certain part of the sequence retained it's secondary structure due to the\n");
  fprintf(nex,"doublet model of substitutions\n");
  fprintf(nex,"=================================PARAMETERS:======================================\n");
  fprintf(nex,"\tTotal length:\t%i\n",slen);
  fprintf(nex,"\tOTU   number:\t%i\n",snum);
  fprintf(nex,"\tstem depth:  \t%i\n",depth);
  fprintf(nex,"max. mutations in stem\t%i\n",msteps2);
  fprintf(nex,"max. mutations in loop\t%i\n",msteps1);
  fprintf(nex,"]\n\nbegin data;\n");
  fprintf(nex,"\tdimensions ntax=%i nchar=%i;\n",snum,slen);
  fprintf(nex,"\tformat datatype=dna gap=- missing=?;\n");
  fprintf(nex,"\tmatrix\n");
  /*-------------------------*/
  for(i=0;i<snum;i++){
    fprintf(fas,">item%i\n",i);
    for(j=0;j<slen;j++) fprintf(fas,"%c",ali[i][j]);
    fprintf(fas,"\n");
  }
  fprintf(phy,"%i  %i\n",snum,slen);
  for(i=0;i<snum;i++){
    sprintf(tmp,"item%i",i);
    while(strlen(tmp)<10){
      strcat(tmp," ");
    }
    fprintf(phy,"%s",tmp);
    fprintf(nex,"\t%s",tmp);
    for(j=0;j<slen;j++){
      fprintf(phy,"%c",ali[i][j]);
      fprintf(nex,"%c",ali[i][j]);
    }
    fprintf(phy,"\n");
    fprintf(nex,"\n");
  }
  fprintf(nex,"\t;\nend;\n\n");
  fprintf(nex,"begin mrbayes;\n");
  fprintf(nex,"\t[Define pairs for the doublet model]\n");
  fprintf(nex,"\tPairs   ");
  for(i=0;i<depth;i++){
    if(i>0) fprintf(nex,", ");
    if(!(i % 12)) fprintf(nex,"\n");
    fprintf(nex,"%i : %i",i+1,2*depth-i);
  }
  fprintf(nex,";\n\n");
  fprintf(nex,"\t[Define character sets]\n");
  fprintf(nex,"\tcharset loops         = ");
  for(i=2*depth;i<slen;i++){
    fprintf(nex,"%i ",i+1);
    if(!(i % 12)) fprintf(nex,"\n");
  }
  fprintf(nex,";\n");
  fprintf(nex,"\tcharset stems                     = ");
  for(i=0;i<2*depth;i++){
    fprintf(nex,"%i ",i+1);
    if(!(i % 12)) fprintf(nex,"\n");
  }
  fprintf(nex,";\n");
  fprintf(nex,"\n\n\tset nowarn=yes;\n\n");
  fprintf(nex,"\t[Define partitions]\n");
  fprintf(nex,"\tpartition smart = 2:stems,loops;\n");
  fprintf(nex,"\tset partition = smart;\n");
  fprintf(nex,"[!   Model taking into account possible co-evolution in stems --------]\n");
  fprintf(nex,"\tlset applyto=(1) nucmodel=doublet;\n");
  fprintf(nex,"\tlset applyto=(2) nucmodel=4by4;\n");
  fprintf(nex,"\tlset nst=1;\n");
  fprintf(nex,"\tmcmc nruns=1 nchains=1 ngen=300000 file=%slps printfreq=1000;\n\tsumt burnin=2000;\n\tsump burnin=2000;\n",root);
  fprintf(nex,"[!   SIMPLE Model  ---------------------------------------------------]\n");
  fprintf(nex,"\tlset applyto=(1) nucmodel=4by4;\n");
  fprintf(nex,"\tlset applyto=(2) nucmodel=4by4;\n");
  fprintf(nex,"\tlset nst=1;\n\n");
  fprintf(nex,"\tmcmc nruns=1 nchains=1 ngen=300000 file=%snols printfreq=1000;\n\tsumt burnin=2000;\n\tsump burnin=2000;",root);
  fprintf(nex,"\nend;\n\n");
  
  fprintf(stt,"Date of analysis:\t%s\n\n",ctime(&tt));
  fprintf(stt,"Length:\t\t\t%i\n",slen);
  fprintf(stt,"Number of sequences:\t%i\n",snum);
  fprintf(stt,"Halflength of loop:\t%i\n",depth);
  fprintf(stt,"Number steps in loop:\t%i\n",msteps1);
  fprintf(stt,"Number steps in stem:\t%i\n",msteps2);

  fprintf(stt,"Total S:\t\t%i\n",divpos(ali,slen,snum,0,slen));
  fprintf(stt,"Stem S:\t\t\t%i\n",divpos(ali,slen,snum,0,depth*2));
  fprintf(stt,"Loop S:\t\t\t%i\n",divpos(ali,slen,snum,depth*2+1,slen));
  fprintf(stt,"Total pi:\t\t%6.5f\n",piCalc(ali,snum,0,slen));
  fprintf(stt,"Stem pi:\t\t%6.5f\n",piCalc(ali,snum,0,depth*2));
  fprintf(stt,"Loop pi:\t\t%6.5f\n",piCalc(ali,snum,depth*2+1,slen));
  S=divpos(ali,slen,snum,0,slen);
  pi=piCalc(ali,snum,0,slen);
  fprintf(stt,"Tajima\'s D total:\t%6.5f\n",tajCalc((double)S,pi,(double)snum));
  S=divpos(ali,slen,snum,0,depth*2);
  pi=piCalc(ali,snum,0,depth*2);
  fprintf(stt,"Tajima\'s D Stem:\t%6.5f\n",tajCalc((double)S,pi,(double)snum));
  S=divpos(ali,slen,snum,depth*2+1,slen);
  pi=piCalc(ali,snum,depth*2+1,slen);
  fprintf(stt,"Tajima\'s D total:\t%6.5f\n",tajCalc((double)S,pi,(double)snum));
  
  /*WAY OUT===================================================================================================*/
  for(i=0;i<snum;i++){ 
    free(ali[i]);
    free(tree[i]);
//    free(ttm[i]);
  }
  /*
  for(i=0;i<nodeCounter;i++) free(nw[i]);
  free(nw);
  */
  free(ali);
  free(tree);
//  free(ttm);
  free(lstr);
  free(starter);
  fclose(fas);
  fclose(phy);
  fclose(nex);
  fclose(struc);
  fclose(tre);
  fclose(stt);
  return(0);
}