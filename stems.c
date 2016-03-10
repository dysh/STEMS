#include "stems.h"

int whichDoublet(int a, int b){
  int i,j;
  
  i=0;j=0;
  while(a != NUCL[i] && i<4) i++;
  while(b != NUCL[j] && j<4) j++;
  return(i*4+j);
}
  
void statistics(otu *src, int start, int end, int sn){
     //system("PAUSE");
  int i,j,k,x1,x2,y1,y2;
  char col[10];
  int count,cc,ccx,ccy,cv,cvx,cvy,vf,vx,vy,a,c,g,t,countx,county;
  point structure[6];
  fpoint fs[6];
  float z1,z2;
  char c1,c2;
  int cnst;
  
  for(i=0;i<6;i++){
    l1[i]=0;
    l2[i]=0;
    structure[i].x=0;
    structure[i].y=0;
    fs[i].x=0;
    fs[i].y=0;
  }
  for(i=0;i<20;i++){
    s1[i]=0;
    s2[i]=0;
  }
  if(end<=start){
    printf("\n\tERROR: the beginning %i is behind of the end %i!\n\n",start,end);
    return;
  }
    
  fprintf(rep,"\n\n<p>----------------------TOTAL-------------------\n");
  count=0;
  countx=0;
  county=0;
  cc=0;
  ccx=0;
  ccy=0;
  cvx=0;
  cvy=0;
  cv=0;
  a=0;c=0;g=0;t=0;
  for(i=start;i<end;i++){
    vf=0;
    vx=0;
    vy=0;
    for(j=1;j<sn;j++){
      c1=src[j].seq[i];
      c1=toupper(c1);
      if(c1=='U') c1='T';
      for(k=0;k<j;k++){
	//printf("%i\t%i\t%i\n",i,j,k);
	c2=src[k].seq[i];
	c2=toupper(c2);
	if(c2=='U'){
	  c2='T';
	  if(!rflag) rflag=1;
	}
	count+=1;
	if(dta->loop[i]!=0) countx +=1; else county += 1;
	if(c1!=c2){
	  vf=1;
	  if(dta->loop[i]!=0) vx=1; else vy=1;
	  if(c1=='A' && c2=='C'){
	    l1[0] +=1;
	    if(dta->loop[i]!=0){
	      structure[0].x +=1;
	    } else {
	      structure[0].y +=1;
	    }
	    break;
	  }
	  if(c1=='A' && c2=='G'){
	    l1[1] +=1;
	    if(dta->loop[i]!=0){
	      structure[1].x +=1;
	    } else {
	      structure[1].y +=1;
	    }
	    break;
	  }
	  if(c1=='A' && c2=='T'){
	    l1[2] +=1;
	    if(dta->loop[i]!=0){
	      structure[2].x +=1;
	    } else {
	      structure[2].y +=1;
	    }
	    break;
	  }
	  if(c1=='C' && c2=='A'){
	    l1[0] +=1;
	    if(dta->loop[i]!=0){
	      structure[0].x +=1;
	    } else {
	      structure[0].y +=1;
	    }
	    break;
	  }
	  if(c1=='C' && c2=='G'){
	    l1[3] +=1;
	    if(dta->loop[i]!=0){
	      structure[3].x +=1;
	    } else {
	      structure[3].y +=1;
	    }
	    break;
	  }
	  if(c1=='C' && c2=='T'){
	    l1[4] +=1;
	    if(dta->loop[i]!=0){
	      structure[4].x +=1;
	    } else {
	      structure[4].y +=1;
	    }
	    break;
	  }
	  if(c1=='G' && c2=='A'){
	    l1[1] +=1;
	    if(dta->loop[i]!=0){
	      structure[1].x +=1;
	    } else {
	      structure[1].y +=1;
	    }
	    break;
	  }
	  if(c1=='G' && c2=='C'){
	    l1[3] +=1;
	    if(dta->loop[i]!=0){
	      structure[3].x +=1;
	    } else {
	      structure[3].y +=1;
	    }
	    break;
	  }
	  if(c1=='G' && c2=='T'){
	    l1[5] +=1;
	    if(dta->loop[i]!=0){
	      structure[5].x +=1;
	    } else {
	      structure[5].y +=1;
	    }
	    break;
	  }
	  if(c1=='T' && c2=='A'){
	    l1[2] +=1;
	    if(dta->loop[i]!=0){
	      structure[2].x +=1;
	    } else {
	      structure[2].y +=1;
	    }
	    break;
	  }
	  if(c1=='T' && c2=='C'){
	    l1[4] +=1;
	    if(dta->loop[i]!=0){
	      structure[4].x +=1;
	    } else {
	      structure[4].y +=1;
	    }
	    break;
	  }
	  if(c1=='T' && c2=='G'){
	    l1[5] +=1;
	    if(dta->loop[i]!=0){
	      structure[5].x +=1;
	    } else {
	      structure[5].y +=1;
	    }
	    break;
	  }
	}
      }
    }
    cc +=1-vf;
    cv +=vf;
    ccx += 1-vx;
    cvx += vx;
    ccy += 1-vy;
    cvy += vy;
  }
  l2[0]=(double)l1[0]/(double)count;
  l2[1]=(double)l1[1]/(double)count;
  l2[2]=(double)l1[2]/(double)count;
  l2[3]=(double)l1[3]/(double)count;
  l2[4]=(double)l1[4]/(double)count;
  l2[5]=(double)l1[5]/(double)count;
  fs[0].x=(float)structure[0].x/(float)countx;
  fs[1].x=(float)structure[1].x/(float)countx;
  fs[2].x=(float)structure[2].x/(float)countx;
  fs[3].x=(float)structure[3].x/(float)countx;
  fs[4].x=(float)structure[4].x/(float)countx;
  fs[5].x=(float)structure[5].x/(float)countx;
  fs[0].y=(float)structure[0].y/(float)county;
  fs[1].y=(float)structure[1].y/(float)county;
  fs[2].y=(float)structure[2].y/(float)county;
  fs[3].y=(float)structure[3].y/(float)county;
  fs[4].y=(float)structure[4].y/(float)county;
  fs[5].y=(float)structure[5].y/(float)county;
  /*TABLE HERE*/
  fprintf(rep,"<table><tr><td>"); 
  fprintf(rep,"<center>\n");
  fprintf(rep,"<table border=\"1\" align=\"left\" cellspacing=\"0\" cellpadding=\"1\" width=\"50\% \">\n");
  fprintf(rep,"<tr bgcolor=\"#aaFFFF\">\n");
  fprintf(rep,"<td> </td><td>A</td><td>C</td><td>G</td><td>T</td></tr>\n");
  fprintf(rep,"<tr><td>A</td><td>0</td><td>%4.4f</td><td>%4.4f</td><td>%4.4f</td></tr>\n",l2[0],l2[1],l2[2]);
  fprintf(rep,"<tr><td>C</td><td></td><td>0</td><td>%4.4f</td><td>%4.4f</td></tr>\n",l2[3],l2[4]);
  fprintf(rep,"<tr><td>G</td><td></td><td></td><td>0</td><td>%4.4f</td></tr>\n",l2[5]);
  fprintf(rep,"<tr><td>T</td><td></td><td></td><td></td><td>0</td></tr>\n");
  fprintf(rep,"</table></center>\n");
  fprintf(rep,"</td></tr>\n<tr><td>\n");
  fprintf(rep,"<h3>Statistics for total sequences</h3>\n<table>\n");
  fprintf(rep,"<tr>\n"); 
  fprintf(rep,"<td bgcolor=\"#ccc\">Transitions:</td><td align=\"right\" bgcolor=\"#008\"><font color=\"lime\">%6.5f</font></td></tr>\n",l2[1]+l2[4]);
  fprintf(rep,"<tr><td bgcolor=\"#ccc\">Transversions:</td><td align=\"right\" bgcolor=\"#008\"><font color=\"lime\">%6.5f</font></td></tr>\n",l2[0]+l2[2]+l2[3]+l2[5]);
  fprintf(rep,"<tr><td bgcolor=\"#ccc\">Ts/Tv ratio:</td><td align=\"right\" bgcolor=\"#008\"><font color=\"lime\">%6.5f</font></td></tr>\n",(l2[1]+l2[4])/(l2[0]+l2[2]+l2[3]+l2[5]));
  fprintf(rep,"<tr><td bgcolor=\"#ccc\">Variable sites:</td><td align=\"right\" bgcolor=\"#008\"><font color=\"lime\">%i</font></td></tr\n",cv);
  fprintf(rep,"<tr><td bgcolor=\"#ccc\">Constant sites:</td><td align=\"right\" bgcolor=\"#008\"><font color=\"lime\">%i</font></td></tr>\n",cc);
  fprintf(rep,"<tr><td bgcolor=\"#ccc\">Ratio of var. sites:</td><td align=\"right\" bgcolor=\"#008\"><font color=\"lime\">%6.5f</font></td></tr>\n",(float)cv/(float)(cc+cv));
  fprintf(rep,"</table>");
  fprintf(rep,"</td></tr>\n<tr><td>\n");
  /*====================================================================================================================*/
  fprintf(rep,"<p>------------------IN-STEMS--------------------</p>\n");
  fprintf(rep,"</td></tr>\n<tr><td>\n");
  fprintf(rep,"<center>\n<table border=\"1\" align=\"left\" cellspacing=\"0\" cellpadding=\"1\" width=\"50\%\">\n");
  fprintf(rep,"<tr bgcolor=\"#FFaaFF\">\n");
  fprintf(rep,"<td> </td><td>A</td><td>C</td><td>G</td><td>T</td></tr>\n");
  fprintf(rep,"<tr></td><td>A</td><td>0</td><td>%4.4f</td><td>%4.4f</td><td>%4.4f</td></tr>\n",fs[0].x,fs[1].x,fs[2].x);
  fprintf(rep,"<tr><td>C</td><td></td><td>0</td><td>%4.4f</td><td>%4.4f</td></tr>\n",fs[3].x,fs[4].x);
  fprintf(rep,"<tr><td>G</td><td></td><td></td><td>0</td><td>%4.4f</td></tr>\n",fs[5].x);
  fprintf(rep,"<tr><td>T</td><td></td><td></td><td></td><td>0</td></tr>\n");
  fprintf(rep,"</table></center>\n");
  fprintf(rep,"</td></tr>\n<tr><td>\n");
  fprintf(rep,"<p>------------------IN-LOOPS--------------------</p>");
  fprintf(rep,"</td></tr>\n<tr><td>\n");
  fprintf(rep,"<center>\n<table border=\"1\" align=\"left\" cellspacing=\"0\" cellpadding=\"1\" width=\"50\%\">\n");
  fprintf(rep,"<tr bgcolor=\"#FFaaFF\">\n");
  fprintf(rep,"<td> </td><td>A</td><td>C</td><td>G</td><td>T</td></tr>\n");
  fprintf(rep,"<tr></td><td>A</td><td>0</td><td>%4.4f</td><td>%4.4f</td><td>%4.4f</td></tr>\n",fs[0].y,fs[1].y,fs[2].y);
  fprintf(rep,"<tr><td>C</td><td></td><td>0</td><td>%4.4f</td><td>%4.4f</td></tr>\n",fs[3].y,fs[4].y);
  fprintf(rep,"<tr><td>G</td><td></td><td></td><td>0</td><td>%4.4f</td></tr>\n",fs[5].y);
  fprintf(rep,"<tr><td>T</td><td></td><td></td><td></td><td>0</td></tr>\n");
  fprintf(rep,"</table></center>");
  fprintf(rep,"</td></tr></table>\n");
  z1=(float)(structure[1].x+structure[4].x)/(float)countx;
  printf("\tTransitions:\t\t%6.5f\n",z1);
  z2=(float)(structure[0].x+structure[2].x+structure[3].x+structure[5].x)/(float)countx;
  printf("\tTransversions:\t\t%6.5f\n",z2);
  printf("\tTs/Tv ratio:\t\t%6.5f\n",z1/z2);
  printf("\tVariable sites:\t\t%i\n",cvx);
  printf("\tConstant sites:\t\t%i\n",ccx);
  printf("\tRatio of var. sites:\t%6.5f\n\n",(float)cvx/(float)(ccx+cvx));
  printf("\n\n<------------------IN-STEMS-------------------->\n");
  printf("\t\tA\tC\tG\tT\n");
  z1=(float)(structure[1].y+structure[4].y)/(float)county;
  printf("\tTransitions:\t\t%6.5f\n",z1);
  z2=(float)(structure[0].y+structure[2].y+structure[3].y+structure[5].y)/(float)county;
  printf("\tTransversions:\t\t%6.5f\n",z2);
  printf("\tTs/Tv ratio:\t\t%6.5f\n",z1/z2);
  printf("\tVariable sites:\t\t%i\n",cvy);
  printf("\tConstant sites:\t\t%i\n",ccy);
  printf("\tRatio of var. sites:\t%6.5f\n\n",(float)cvy/(float)(ccy+cvy));
  for(i=start;i<end;i++){
    for(j=0;j<sn;j++){
      c1=src[j].seq[i];
      c1=toupper(c1);
      switch(c1){
	case 'A' : a +=1;
	case 'C' : c +=1;
	case 'G' : g +=1;
	case 'T' : t +=1;
      }
    }
  }
  printf("\ta:\t%6.5f\n",((float)a/(float)(a+c+g+t)));
  printf("\tc:\t%6.5f\n",((float)c/(float)(a+c+g+t)));
  printf("\tg:\t%6.5f\n",((float)g/(float)(a+c+g+t)));
  printf("\tt:\t%6.5f\n\n",((float)t/(float)(a+c+g+t)));
  
  fprintf(rep,"<H2>SVG scheme of the sequence</h2>\n");
  fprintf(rep,"<p>Each bar is a nucleotide. Red bars designate ascending arm of a loop, \
  the blue ones designate descending arm. Upper line shows polymorphic positions as open circles, \
  and positions containing indels as filled circles</p>\n");
  fprintf(rep,"<p>More graphic information can be found in file <b>results/info.svg</b></p>\n");
  svg=fopen("results/info.svg","w");
  fprintf(svg,"<?xml version=\"1.0\" standalone=\"no\"?>\n <!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n");
  y1=(end/100+1)*50;
  fprintf(rep,"<svg id=\"svgelem\" height=\"%i\" xmlns=\"http://www.w3.org/2000/svg\">\n", y1);
  fprintf(svg,"<svg width = \"12cm\" height=\"20cm\" viewBox= \"0 0 ");
  fprintf(rep,"\n<g style=\"font-family: sans-serif;font-size: 10pt;\"\n");
  for(i=start;i<end;i++){
    x1=50+ 4*i % 400;
    y1= 50+ 40* (int)(i/100);
    //printf("%i\t%i\t%i\n",i,x1,y1);
    x2=x1;
    if(dta->loop[i]>0){
      y2=y1-5;
      sprintf(col,"%s","red");
    }
    if(dta->loop[i]<0){
      y2=y1+5;
      sprintf(col,"%s","blue");
    }
    if(dta->loop[i]==0){
      y2=y1+1;
      y1--;
      sprintf(col,"%s","black");
    }
    if (i%100==1){
      fprintf(rep,"<text x=\"10\" y=\"%i\" style=\"font-family: sans-serif;font-size: 10pt;\"> %i</text>\n",y1,i);
     //printf("<text x=\"0\" y=\"%i\"> %i</text>\n",y1,i);
    }
    fprintf(rep,"  <line x1=\"%i\" y1=\"%i\" x2=\"%i\" y2=\"%i\" style=\"stroke:%s;stroke-width:3\"/>\n",x1,y1,x2,y2,col);
    //полиморфные положения надо отметить, а также положения с дырками
    c1=src[0].seq[i];
    c1=toupper(c1);
    cnst=0;
    for(j=1;j<sn;j++){
      c2=src[j].seq[i];
      c2=toupper(c2);
      if(c1!=c2 && cnst==0) cnst=1;
      if(c1=='-' || c2=='-') cnst=2;
    }
    if(cnst==0) fprintf(rep,"  <line x1=\"%i\" y1=\"%i\" x2=\"%i\" y2=\"%i\" style=\"stroke:green;stroke-width:1\"/>\n",x1,y1-10,x2,y1-7);
    if(cnst==1) fprintf(rep,"  <circle cx=\"%i\" cy=\"%i\" r=\"3\" style=\"stroke:magenta;stroke-width:1;fill:none\"/>\n",x1,y1-10);  
    if(cnst==2) fprintf(rep,"  <circle cx=\"%i\" cy=\"%i\" r=\"3\" style=\"stroke:cyan;stroke-width:1;fill:yellow\"/>\n",x1,y1-12);
  }
  for(i=0;i<10;i++) fprintf(rep,"<line x1=\"%i\" y1=\"%i\" x2=\"%i\" y2=\"%i\" style=\"stroke:black;stroke-width:0.5; \
			    stroke-opacity:0.6; stroke-dasharray:2 2;\"/>\n",50+40*i,10,50+40*i,y2+10);
  //fprintf(rep,"</g>\n");
  fprintf(rep,"</svg>");
  return;
}

void estat(otu *src, int start, int end, int sn, int nl){
  int i,j,k;
  int count,cc,ccx,ccy,cv,cvx,cvy,vf,vx,vy,a,c,g,t,countx,county;
  point structure[6];
  //fpoint fs[6];
  //float z1,z2;
  char c1,c2;
  
  count=0;
  countx=0;
  county=0;
  cc=0;
  ccx=0;
  ccy=0;
  cvx=0;
  cvy=0;
  cv=0;
  a=0;c=0;g=0;t=0;

  for(i=start;i<end;i++){
    vf=0;
    vx=0;
    vy=0;
    for(j=1;j<sn;j++){
      c1=src[j].seq[i];
      c1=toupper(c1);
      for(k=0;k<j;k++){
        //printf("%i\t%i\t%i\n",i,j,k);
        c2=src[k].seq[i];
        c2=toupper(c2);
        count+=1;
        if(dta->loop[i]!=0) countx +=1; else county += 1;
        if(c1!=c2){
          vf=1;
          if(dta->loop[i]!=0) vx=1; else vy=1;
          if(c1=='A' && c2=='C'){
            l1[0] +=1;
            if(dta->loop[i]!=0){
              structure[0].x +=1;
            } else {
              structure[0].y +=1;
            }
            break;
          }
          if(c1=='A' && c2=='G'){
            l1[1] +=1;
            if(dta->loop[i]!=0){
              structure[1].x +=1;
            } else {
              structure[1].y +=1;
            }
            break;
          }
          if(c1=='A' && c2=='T'){
            l1[2] +=1;
            if(dta->loop[i]!=0){
              structure[2].x +=1;
            } else {
              structure[2].y +=1;
            }
            break;
          }
          if(c1=='C' && c2=='A'){
            l1[0] +=1;
            if(dta->loop[i]!=0){
              structure[0].x +=1;
            } else {
              structure[0].y +=1;
            }
            break;
          }
          if(c1=='C' && c2=='G'){
            l1[3] +=1;
            if(dta->loop[i]!=0){
              structure[3].x +=1;
            } else {
              structure[3].y +=1;
            }
            break;
          }
          if(c1=='C' && c2=='T'){
            l1[4] +=1;
            if(dta->loop[i]!=0){
              structure[4].x +=1;
            } else {
              structure[4].y +=1;
            }
            break;
          }
          if(c1=='G' && c2=='A'){
            l1[1] +=1;
            if(dta->loop[i]!=0){
              structure[1].x +=1;
            } else {
              structure[1].y +=1;
            }
            break;
          }
          if(c1=='G' && c2=='C'){
            l1[3] +=1;
            if(dta->loop[i]!=0){
              structure[3].x +=1;
            } else {
              structure[3].y +=1;
            }
            break;
          }
          if(c1=='G' && c2=='T'){
            l1[5] +=1;
            if(dta->loop[i]!=0){
              structure[5].x +=1;
            } else {
              structure[5].y +=1;
            }
            break;
          }
          if(c1=='T' && c2=='A'){
            l1[2] +=1;
            if(dta->loop[i]!=0){
              structure[2].x +=1;
            } else {
              structure[2].y +=1;
            }
            break;
          }
          if(c1=='T' && c2=='C'){
            l1[4] +=1;
            if(dta->loop[i]!=0){
              structure[4].x +=1;
            } else {
              structure[4].y +=1;
            }
            break;
          }
          if(c1=='T' && c2=='G'){
            l1[5] +=1;
            if(dta->loop[i]!=0){
              structure[5].x +=1;
            } else {
              structure[5].y +=1;
            }
            break;
          }
        }
      }
    }
    cc +=1-vf;
    cv +=vf;
    ccx += 1-vx;
    cvx += vx;
    ccy += 1-vy;
    cvy += vy;
  }
  return;
}

int insdot(char *t, int pos){
  int l;
  char dot[] = ".";
  l=strlen(t);
  memmove(t+pos+1,t+pos,l-pos+1);
  memcpy(t+pos,dot,1);
  return(0);
}

int ins0(int *t, size_t l, size_t pos){
  int i;
  
  if(pos>l) return(-1);
  for(i=l;l>pos;i--){
    t[i]=t[i-1];
  }
  t[pos]=1;
  return 0;
}

int main(int argc, char *argv[]){
  char qz,qzz;
  
  iflag=0;
  lflag=0;
  uflag=0;
  eflag=0;
  sflag=0;
  rflag=0;
  
  /*в виндах директорию не получается сделать! Надо что-то придумать...*/
  if(stat("results",&fbuf)<0){
    printf("\n\n\tDirectory \"results\" does not exist, creating\n");
    i=  mkdir("./results",S_IRWXU | S_IRWXG | S_IRWXO);
  }
  rep=fopen("report.html","w");
  printf("\tDumping results to \"report.html\"\n");
  fprintf(rep,"<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n<html>\n<head>\n");
  fprintf(rep,"<title>Results of the last run. Next run will overwrite this!</title>\n");
  fprintf(rep,"<meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\">\n</head>\n\n");
  fprintf(rep,"<body bgcolor=\"#FFFFFF\">\n");
  fprintf(rep,"<STYLE TYPE=\"text/css\">\n<!--\n \
    .hh {\n \
    background-color:#888844;\n \
    color:white;\n \
    }\n \
    rcell {\n \
    background-color:#0000ff;\n \
    color:white;\n \
    align:right;\n \
    }\n \
    .wh {\n  \
    background-color=white;\n \
    color=#ffff00;\n }\n \
    -->\n</STYLE>\n");
    
  fclose(rep);
  rep=fopen("report.html","a");
  time(&lt);
  fprintf(rep,"<table bgcolor=\"#FFFFFF\">\n<tr><td>");
  
  fprintf(rep,"<h2>Date and time: %s</h2>\n",ctime(&lt)); 
  strcpy(info.outfile,"results/stems.nex");
  printf("\treading arguments...\n");
  while ((rez=getopt_long(argc,argv,short_options,long_options,&option_index))!=-1){
    strcpy(info.alignment,"");
    strcpy(info.structure,"");
    //strcpy(info.outfile,"");
    switch(rez){
      case 'h' : {
	  printf ("\n\n\tProgram writes the \"stem block\" for MrBayes\n");
	  printf ("\tthere is special block called \"Stems\"\n");
	  printf ("\tand table of doublets, which are the bases\n");
	  printf ("\topposing each other in loops as described\n");
	  printf ("\tin \"dot-bracket\" file produced for example\n");
	  printf ("\tby Vienna RNA package\n\n");
	  printf ("OPTIONS:\n");
	  printf ("\t -h or --help\t\t\t\tthis message\n");
	  printf ("\t -v or --version\t\t\tprogram version, if anybody cares\n");
	  printf ("========================================================================\n\n");
	  printf ("\t -a or --clustal_file <filename.aln>\tclustal file with alignment, REQUIRED\n");
	  printf ("\t -s or --structure_file <filename>\ttwo line output like is produced with \n\t\t\t\t\t\tRNAfold from Vienna RNA package. REQUIRED\n");
	  printf ("\t -o or --output_file <filename.nex>\tNEXUS file to be made. OPTIONAL\n\n");
	  printf ("\t -i or --single_structure <number>\tloop structure for a lingle OTU is\n\t\t\t\t\t\tused. Its number is given as a parameter. Gaps will be\n\t\t\t\t\t\tintroduced into the template before writing the NEXUS file\n");
	  printf ("\t -l or --separate_loops <number>\tspearate stretches bearing loops\n\t\t\t\t\t\tlonger then the parameter are extracted\n\t\t\t\t\t\tfrom the original alignment and put into separate\n\t\t\t\t\t\tNEXUS files thus allowing one to\n\t\t\t\t\t\texamine them separately\n");
	  printf ("\t -x or --unique_OTUs <number>   \tleave only unique OTUs. In case of argument > 0\n\t\t\t\t\t \
	the haplotypes differing by less then N  mutational\n \
	  \t\t\t\t\tsteps will be considered equal\n");
	  printf ("========================================================================\n\n");
	  printf ("\t -u or --otu_list\t\t\tprints out numbered otu list and quits\n");
	  printf ("\t -e or --extract_fasta <number>\t\textracts a single sequence defined by it\'s number\n");
	  printf ("\t -t or --statistics\t\t\tdisplays substitution statistics on the screen\n\n");
	  fprintf(rep,"Nothing was done\n</body>\n</html>");
	  fclose(rep);
	  return 0;
	  break;
      }
      case 'v' : {
	  printf("\nVersion %s\n\n", VERSION);
	  fprintf(rep,"Vesrion %s\n</body>\n</html>", VERSION);
          fprintf(rep,"\n</body>\n</html>\n");
	  fclose(rep);

	  return 0;
	  break;
      }
      case 'a' : {
	  if(optarg != NULL) {
	    strcpy(info.alignment,optarg);
	    fprintf(rep,"<p>Clustal file with aligned sequences: <font color = \"#2222FF\">%s</font><br>\n",info.alignment);
	    in1 = fopen(optarg,"r");
	  } else {
	    printf("\n\n\t\t---->>>> You forgot to enter filename for the alignment! <<<<-----\n\t\t\tQUITTING!\n");
            fprintf(rep,"<p>You forgot to enter filename for the alignment</p>\n");
            fprintf(rep,"\n</body>\n</html>\n");
            fclose(rep);
	    return -1;
	  }
	  break;
      }
      case 's' : {
	  if(optarg != NULL) {
	    strcpy(info.structure,optarg);
	    in = fopen(optarg,"r");
	  } else {
	    printf("\n\n\t\t---->>>> You forgot to enter filename for the alignment! <<<<-----\n\t\t\tQUITTING!\n");
            fprintf(rep,"<p>You forgot to enter filename </p>\n");
            fprintf(rep,"\n</body>\n</html>\n");
            fclose(rep);
	    return -1;
	  }
	  break;
      }
      
      case 'o' : {    /*херня какая-то!*/
	  if(optarg != NULL) {
	    strcpy(info.outfile,optarg);
	  }else {
	    strcpy(info.outfile,optarg);
	    //out = fopen("stem.nex","w");
	  }
	  break;
      }
      case 'i' : {
      if(optarg!=NULL){
        nref=atoi(optarg);
        iflag=1;
      } else {
        printf("You should specify the number of the reference OTU for which the secondary structure is obtained!\nQUITTING!\n");
        fprintf(rep,"<p>You should specify the number of the reference OTU for which the secondary structure is obtained!</p>\n");
        fprintf(rep,"\n</body>\n</html>\n");
        fclose(rep);
        return(-1);
      }
      break;
      }
      case 'l' : {
      if(optarg!=NULL){
        minloop=atoi(optarg);
        lflag=1;
      } else {
        printf("You should specify the length of the shortest loop you would like to examine\nQUITTING!\n");
        fprintf(rep,"<p>You should specify the length of the shortest loop you would like to examine</p>\n");
        fprintf(rep,"\n</body>\n</html>\n");
        fclose(rep);
        return(-1);
      }
      break;
      }
      case 'u' : {
        uflag=1;
        break;
      }
      case 'e' : {
	eflag=atoi(optarg);
	if(eflag==0){
	  printf("eflag=0, don\'t know what to extract,\nQUITTING\n");
          fclose(rep);
	  return -1;
	}
      }
      case 't' : {
	sflag=1;
      }
      default: break;
    }
  }
  if(in==NULL && !uflag){
    printf("\n\ncould not open file %s,quitting\n",info.alignment);
    fprintf(rep,"<p>could not open file file with alignment: %s,quitting</p>",info.alignment);
    fclose(rep);
    return -1;
  }
  if(in1==NULL){
    printf("\n\ncould not open file containing secondary structure %s, quitting\n",info.structure);
    return -1;
  }
  /*это если надо только список названий с номерами, а потом - заглохнуть*/
  if(uflag==1){
    count=0;
    fgets(line,1024,in1);
    while(!feof(in1)){
      fgets(line,1024,in1);printf("\tstem positions:\t%i\n\tloop positions:\t%i\n\ttotal length\t%i\n",bo+bc,nd,i);
 
      if(strlen(line)<5 || isspace(line[0])) {
	if(count>0) break;
	count=0;
      } else {
	count++;
	sscanf(line,"%s %s",name,seq);
	printf("\t%i\t%s\n",count,name);
	if(count>nspec) nspec=count;
      }
    }
    rewind(in1);
    return(0);  
  }
  if(eflag>0){
    count=0;
    strcpy(t,"");
    fgets(line,1024,in1);
    while(!feof(in1)){
      fgets(line,1024,in1);
      if(strlen(line)<5 || isspace(line[0])) {
	//if(count>0) break;
	count=0;
      } else {
	count++;
	sscanf(line,"%s %s",name,sq);
	//printf("\t%i\t%s\n",count,name);
	if(count==eflag){
	  if(strlen(t)<1){
	    sprintf(t,"%s.fas",name);
	    tmp=fopen(t,"w");
	    fprintf(tmp,">%s\n",name);
	    fprintf(tmp,"%s",sq);
	  } else {
	    fprintf(tmp,"%s",sq);
	  }
	}
	if(count>nspec) nspec=count;
      }
    }
    fprintf(tmp,"\n");
    fclose(tmp);
    rewind(in1);
    return(0);  
  }
  printf("\tnow reading sequences\n");
  dta=malloc(sizeof(dta));
  fstat(fileno(in),&fbuf);
  ofsize=fbuf.st_size;
  ifile=malloc(ofsize*sizeof(char));
  sfile=malloc(ofsize*sizeof(char));
  pairs=malloc(ofsize*sizeof(fpoint));
  fgets(sfile,ofsize,in);
  fgets(ifile,ofsize,in);
  
  sscanf(ifile,"%s",ifile);
  dta->slen=strlen(ifile);
  dta->loop=malloc(dta->slen*sizeof(int));
  dta->brackets=malloc(dta->slen*sizeof(int));
  dta->ranks=malloc(dta->slen*sizeof(int));
  j=0;
  bo=0;
  bc=0;
  nd=0;
  r=0;
 
  for(i=0; i<dta->slen; i++){
    c=ifile[i];
    switch(c){
      case '.':
	nd++;
	dta->loop[i]=0;
	dta->ranks[i]=0;
	break;
      case ('('): 
	bo+=1;
        r+=1;
	dta->loop[i]=bo;
//	r+=1;
	dta->ranks[i]=r;
	break;
    
      case ')':
	bc+=1;
	dta->loop[i]=-bc;     /*ZIC!!!!!!*/
	dta->ranks[i]=r;
      //  bc+=1;
      r--;

	break;
    }
    dta->brackets[i]=r;
  }
  fprintf(rep, "<p><table border=\"0\" align=\"left\" cellspacing=\"1\" cellpadding=\"1\" width=\"100\%\">\n");
  fprintf(rep, "<tr><td bgcolor=\"#cccccc\"> positions in stems:</td><td align=\"right\" bgcolor=\"#008\"><font color=\"lime\">%i</font></td></tr>\n",bo+bc);
  fprintf(rep, "<tr><td bgcolor=\"#cccccc\"> positions in loops (unpaired):</td><td align=\"right\" bgcolor=\"#008\"><font color=\"lime\">%i</font></td></tr>\n",nd);
  fprintf(rep, "<tr><td bgcolor=\"#cccccc\"> total length:</td><td align=\"right\" bgcolor=\"#008\"><font color=\"lime\">%i</font></td></tr>\n",i);
  fprintf(rep,"</table></p>\n");
  seq=NULL;
  printf("\tcontrol values are \t%i\t%i\n",bo-bc,r);
  
  count=0;
  while(!feof(in1)){
    fgets(line,1024,in1);
    if(strlen(line)<5 || isspace(line[0])) {
      count=0;
    } else {
      count++;
      if(count>nspec) nspec=count;
    }
  }
  fprintf(rep,"<p> </p><p><b>Data file contains <font color = \"#2222FF\">%i</font> seqiences</b></p>\n",nspec);

  dta->snum=nspec;
  dta->ali=malloc((dta->snum+1)*sizeof(otu));
  
  for(i=0;i<dta->snum;i++){
    dta->ali[i].seq=malloc(dta->slen*sizeof(int));
    strcpy(dta->ali[i].seq,"");
    strcpy(dta->ali[i].nme,"");
  }
  rewind(in1);
  fgets(line,1024,in1);
  count=0;
  seq=malloc(dta->slen*sizeof(int));
  while(!feof(in1)){
    fgets(line,1024,in1);
    if(strlen(line)<4 || isspace(line[0])) {
      count=0;
    } else {
      sscanf(line,"%s %s",name,seq);
      if(strlen(dta->ali[count].nme)<1) strcpy(dta->ali[count].nme,name);
      strcat(dta->ali[count].seq,seq);
      count++;
    }
  }  

  fprintf(rep,"<p>STARTING MAIN NEXUS FILE <b>%s</b> \n",info.outfile);
  out = fopen(info.outfile,"w");
  dta->slen=strlen(dta->ali[0].seq);
  /*
  fprintf(out,"\nifile");
  for(i=0;i<dta->slen;i++) fprintf(out,"\t%c",ifile[i]);
  fprintf(out,"\nloop");
  for(i=0;i<dta->slen;i++) fprintf(out,"\t%i",dta->loop[i]);
  fprintf(out,"\nranks");
  for(i=0;i<dta->slen;i++) fprintf(out,"\t%i",dta->ranks[i]);
  fprintf(out,"brack\n");
  for(i=0;i<dta->slen;i++) fprintf(out,"\t%i",dta->brackets[i]);
  fprintf(out,"\n   i:");
  for(i=0;i<dta->slen;i++) fprintf(out,"\t%i",i);
  */
  fprintf(out,"#NEXUS\n[created with the loop/stem utility]\n\nbegin data;\n");
  fprintf(out,"\tdimensions ntax=%i nchar=%i;\n",dta->snum,dta->slen);
  if(rflag) fprintf(out,"\tformat datatype=rna gap=- missing=?;\n");
  else fprintf(out,"\tformat datatype=dna gap=- missing=?;\n");
  fprintf(out,"\tmatrix\n");
  for(j=0;j<dta->snum;j++){
    fprintf(out,"\t%s  \t%s\n",dta->ali[j].nme,dta->ali[j].seq);
  }
  fprintf(out,"\t;\nend;\n\n");

  count=0;
  fprintf(out,"begin mrbayes;\n");
  fprintf(out,"\t[Define pairs for the doublet model]\n");
  fprintf(out,"\tPairs   ");

  for(i=0;i<dta->slen-1;i++){
  
    if(dta->loop[i]>0){
      fprintf(out,"%i:",i+1);
      j=i+1;
      flag=1;
      count++;
      while(flag && j<dta->slen){
	
	flag=dta->brackets[j]-dta->brackets[i]+1;
        j++;
      }
      fprintf(out,"%i",j);
//      loo[j]=0;
      if(count<bo) fprintf(out,", "); else {
	fprintf(out,"; \n\n");
	break;
      }
      //ranks[j]=0;
      if (!(count % 8)) fprintf(out,"\n\t        ");
    }
  }
  
  fprintf(out,"\t[Define character sets]\n");
  fprintf(out,"\tcharset loops         = ");
  count=0;
  for(i=0;i<dta->slen;i++){
    if(dta->loop[i]==0){
      fprintf(out,"%i ",i+1);
      count++;
      if (!(count % 20)) fprintf(out,"\n\t          ");
    }
  }

  fprintf(out,";\n");
  fprintf(out,"\tcharset stems                     = ");
  count=0;
  for(i=0;i<dta->slen;i++){
    if(dta->loop[i]!=0){
      fprintf(out,"%i ",i+1); 
      count++;
      if (!(count % 20)) fprintf(out,"\n\t          ");
    }
  }
  fprintf(out,";\n\n");
  fprintf(out,"\tset nowarn=yes;\n\n");
  fprintf(out,"\t[Define partitions]\n");
  fprintf(out,"\tpartition smart = 2:stems,loops;\n");
  fprintf(out,"\tset partition = smart;\n");
  fprintf(out,"[!   Model taking into account possible co-evilution in stems --------]\n");
  fprintf(out,"\tlset applyto=(1) nucmodel=doublet;\n");
  fprintf(out,"\tlset applyto=(2) nucmodel=4by4;\n");
  fprintf(out,"\tlset nst=6;\n\n");
  fprintf(out,"\tprset ratepr=variable;\n\n");
  fprintf(out,"\tmcmcp nruns=2 nchains=4 ngen=500000 file=loops.mb;\n\tsumt burnin=2000;\n\tsump burnin=2000;\nss\n");
  fprintf(out,"[!   SIMPLE Model  ---------------------------------------------------]\n");
  fprintf(out,"\tlset applyto=(1) nucmodel=4by4;\n");
  fprintf(out,"\tlset nst=6;\n\n");
  fprintf(out,"\tprset ratepr=variable;\n\n");
  fprintf(out,"\tmcmcp nruns=4 nchains=4 ngen=500000 file=noloops.mb;\n\tsumt burnin=2000;\n\tsump burnin=2000;\nss\n");
  fprintf(out,"end;\n\n");
 
  fclose(in1);
  printf("\n\tMain nexus file DONE\n");
  fprintf(rep,"<p>Nexus file with the 2 hypotheses: <b>%s </b> done",info.outfile);
  
  if(sflag) statistics(dta->ali,0,dta->slen,dta->snum);
  /////////////////////////////////
  if(lflag){
    printf("\n\t\tNOW PREPARING INDIVIDUAL LOOPS\n");
    for(j=1;j<dta->slen-1;j++){
      if((dta->brackets[j]<dta->brackets[j-1]) &&(dta->brackets[j]<=dta->brackets[j+1])){
	i=j-1;
	loopflag=0;
	lloop=0;
	//printf("Local minimum at %i\t with ranks = %i\n",j,dta->brackets[j]);
	while(dta->brackets[i]>dta->brackets[j] && i>=0){ /*ZIC*/
	  if(dta->brackets[i]>dta->brackets[i-1]) loopflag=1;
	  if(dta->brackets[i]<dta->brackets[i-1]) lloop=1;
	  i-=1;
	}
	i +=2;
	//printf("i: %i\tj: %i\tloopflag: %i\tlloop: %i\tranks[i]: %i\tranks[j]: %i\n",i,j,loopflag,lloop,dta->ranks[i],dta->ranks[j]);
	if(loopflag && lloop && j-i>minloop){
	  esta.len=(float)(j-i);
          esta.start=j;
          esta.end=i;
          esta.ts[0]=0;
          esta.ts[1]=0;
          esta.tv[0]=0;
          esta.tv[1]=0;
	  sprintf(tfname,"results/loop%i.struct",nloops+1);
	  tmp=fopen(tfname,"w");
	  for(k=i;k<j;k++) fprintf(tmp,"%c",dta->ali[1].seq[k]);
	  fprintf(tmp,"\n");
	  for(k=i;k<j;k++) fprintf(tmp,"%c",ifile[k]);
	  fprintf(tmp,"\n");
	  fclose(tmp);
	  sprintf(tfname,"./results/stem%i.nex",nloops);
	  tmp=fopen(tfname,"w");
	  fprintf(tmp,"#NEXUS\n[created with the loop/stem utility]\n[hairpin #%i at %i -> %i]\nbegin data;\n",nloops+1,i,j);
	  fprintf(tmp,"\tdimensions ntax=%i nchar=%i;\n",dta->snum,j-i);
	  fprintf(tmp,"\tformat datatype=dna gap=- missing=?;\n");
	  fprintf(tmp,"\tmatrix\n");
	  for(k=0;k<nspec;k++){
	    fprintf(tmp,"\t%s   \t",dta->ali[k].nme);
	    for(l=i;l<j;l++){
	      fprintf(tmp,"%c",dta->ali[k].seq[l]);
	    }
	    fprintf(tmp, "\n");
	  }
	  fprintf(tmp,"\t;\nend;\n\n");
	  count=0;
	  fprintf(tmp,"begin mrbayes;\n");
	  fprintf(tmp,"\t[Define pairs for the doublet model]\n");
	  fprintf(tmp,"\tPairs   ");
	  npairs=0;
	  for(k=i;k<j;k++){
	    
	    if(dta->loop[k]>0){
	      if(npairs>0) fprintf(tmp,",");
	      fprintf(tmp,"%i:",k-i+1);
	      l=k+1;
	      flag=1;
	      count++;
	      while(flag && l<j){
		flag=dta->brackets[l]-dta->brackets[k]+1;
		l++;
	      }
	      npairs++;
	      fprintf(tmp,"%i",l-i);
	      if(count<bo-1) fprintf(tmp," "); else {
		fprintf(tmp,"; \n\n");
		break;
	      }
	      if (!(count % 8)) fprintf(tmp,"\n\t        ");
	    }
	  }
	  fprintf(tmp,";\n\n");
	  fprintf(tmp,"\t[Define character sets]\n");
	  fprintf(tmp,"\tcharset loops         = ");
	  count=0;
	  for(k=i;k<j;k++){
	    if(dta->loop[k]==0){
	      fprintf(tmp,"%i ",k-i+1);
	      count++;
	      if (!(count % 20)) fprintf(tmp,"\n\t          ");
	    }
	  }
	  fprintf(tmp,";\n");
	  fprintf(tmp,"\tcharset stems                     = ");
	  count=0;
	  for(k=i;k<j;k++){
	    if(dta->loop[k]!=0){
	      fprintf(tmp,"%i ",k-i+1); 
	      count++;
	      if (!(count % 20)) fprintf(tmp,"\n\t          ");
	    }
	  }
	  fprintf(tmp,";\n\n");
	  fprintf(tmp,"\tset nowarn=yes;\n\n");
	  fprintf(tmp,"\t[Define partitions]\n");
	  fprintf(tmp,"\tpartition smart = 2:stems,loops;\n");
	  fprintf(tmp,"\tset partition = smart;\n");
	  fprintf(tmp,"[!   Model taking into account possible co-evilution in stems --------]\n");
	  fprintf(tmp,"\tlset applyto=(1) nucmodel=doublet;\n");
	  fprintf(tmp,"\tlset applyto=(2) nucmodel=4by4;\n");
	  fprintf(tmp,"\tlset nst=6;\n\n");
	  fprintf(tmp,"\tprset ratepr=variable;\n\n");
	  fprintf(tmp,"\tmcmc nruns=2 nchains=4 ngen=500000 file=loops%i.mb;\n\tsumt burnin=2000;\n\tsump burnin=2000;\n",nloops);
	  fprintf(tmp,"[!   SIMPLE Model  ---------------------------------------------------]\n");
	  fprintf(tmp,"\tlset applyto=(1) nucmodel=4by4;\n");
	  fprintf(tmp,"\tlset nst=6;\n\n");
	  fprintf(tmp,"\tprset ratepr=variable;\n\n");
	  fprintf(tmp,"\tmcmc nruns=4 nchains=4 ngen=500000 file=noloops%i.mb;\n\tsumt burnin=2000;\n\tsump burnin=2000;\n",nloops);
	  fprintf(tmp,"end;\n\n");

	  fclose(tmp);
	  //if(sflag) statistics(dta->ali,i,j,dta->snum);
	  nloops++;
	
	}
      }
    }
  }
  fprintf(rep,"</td></tr></table>\n");
  fprintf(rep,"\n</body>\n</html>\n");
  printf("\nfinished successfully\n");

  fclose(out);
  fclose(rep);
  free(pairs);
  //if(seq!=NULL) free(seq);

  return(0);
}
