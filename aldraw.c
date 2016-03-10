/*
	Это - тоже самое, что и alidraw, но с использованием gd
	Поэтому компилировать надо скриптом cmpl
	
	*/
	
#include <gd.h>
#include <gdfonts.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int snum,slen,nnlen;
char **nmes, **chars, *indels;

int chkformat(FILE *f){
	int res,i,j,k,n,isHeader;
	char ch;
	
	rewind(f);
	i=0;j=0;isHeader=0;k=0;res=0;n=0;
	nnlen=0;
	while(!feof(f)){
		ch = fgetc(f);
		if(isHeader==1 &&  ch=='\n') {
			isHeader = 0;
			j=0;
			if(n>nnlen) nnlen=n;
		}
		if(isHeader==1 && ch!='\n') n++;
		if(ch=='>'){
			isHeader=1;
			i++;
			if(k==0) k=j;
			if(j!=k) res=1;
			n=0;
		}
		if (!isHeader && (ch=='A' || ch=='C' || ch=='G' || ch=='T' || ch=='a' || ch=='c' || ch=='g' || ch=='t' || ch=='-')) j++;
	}
	if(res==1){
		printf("\n\tSequences are not aligned!\n");
		return(1);
	}
	printf("\nThere are %i sequences %i bp long in the file,\nthe longest name is %i characters\n",i,k,nnlen);
	snum=i;
	slen=k;
	return(0);
}

int getali(FILE *f){
	int i,j,k,isH;
	char c;
	i=-1;j=0;k=0;isH=0;
	rewind(f);
	while(!feof(f)){
		c=fgetc(f);
		if(isH==1 &&  c=='\n') {
			isH = 0;
			j=0;
		}
		if(isH==1 && c!='\n' && c!='>') {
			nmes[i][k]=c;
			k++;
		}
		if (!isH && (c=='A' || c=='C' || c=='G' || c=='T' || c=='a' || c=='c' || c=='g' || c=='t' || c=='-')){
			chars[i][j]=c;
			j++;
		}
		if(c=='>'){
			isH=1;
			i++;
			k=0;
		}
	}
	
	return(0);
}

int main(int argc, char *argv[]){
    gdImagePtr image;
	int i,j,rl;
	char s[255];
    FILE *grout,*ffile;
	
    int white,black,gray1,gray2,gray3,cl;
    
    ffile=fopen(argv[1],"r");
    if(ffile==NULL){
	printf("\nFile %s does not exist in current directory\nQuitting\n",argv[1]);
	return(1);
    }
	chkformat(ffile);
	printf("\n============================================\n");
	nmes = malloc(snum*sizeof(char *));
	for(i=0;i<snum;i++) nmes[i]=malloc(nnlen*sizeof(char));
	chars = malloc(snum*sizeof(char *));
	for(i=0;i<snum;i++) chars[i]=malloc(slen*sizeof(char));
	
	getali(ffile);
	printf("\nParsing fasta file is finished\n");
    image = gdImageCreate(800,800);
    white = gdImageColorAllocate(image, 255,255,255);
    black = gdImageColorAllocate(image,0,0,0);
	gray1 = gdImageColorAllocate(image,64,64,64);
	gray2 = gdImageColorAllocate(image,192,192,192);
	gray3 = gdImageColorAllocate(image,240,240,240);
	for(i=0;i<snum;i++){
		gdImageString(image, gdFontGetSmall(),20,20+i*20,(char *) nmes[i],black);
		gdImageLine(image,30,21+i*20,160,21+i*20,gray2);
		if(slen>300) gdImageLine(image,30,41+(i+snum)*20,160,41+(i+snum)*20,gray2);
		if(slen>650) gdImageLine(image,30,61+(i+snum*2)*20,160,61+(i+snum*2)*20,gray2);
	}
	//gdImageFilledRectangle(image,150,20,750,20+snum*20,gray3);
	for(i=0;i<snum;i++){
		for(j=0;j<300;j++){
			if(j==slen) break;
			if(i==0){
				if(chars[0][j]!='-') gdImageFilledRectangle(image,150+j*2,20+20*i,151+j*2,20+20*i+19,gray1);
				else gdImageFilledRectangle(image,150+j*2,20*i+20,151+j*2,20+20*i+19,white);
			} else {
				if(chars[i][j]!=chars[0][j]) cl=gray2;
				if(chars[i][j]==chars[0][j]) cl=gray1;
				if(chars[i][j]=='-') cl=white;
				gdImageFilledRectangle(image,150+j*2,20+20*i,151+j*2,20+20*i+19,cl);
			
			}
			if(!(j % 25)) {
				gdImageLine(image,150+j*2,20*(snum+1)+2,150+j*2,20*(snum+1),black);
				sprintf(s,"%i",j);
				gdImageString(image, gdFontGetSmall(),150+j*2,20*(snum+1)+4,s,black);
			}
		}
		gdImageLine(image,150,20*(snum+1)+2,750,20*(snum+1)+2,black);
		gdImageLine(image,150,20*(snum+1)+3,750,20*(snum+1)+3,black);
		if(slen>300){
			for(j=300;j<650;j++){
				if(j==slen) break;
				if(i==0){
					if(chars[0][j]!='-') gdImageFilledRectangle(image,50+(j-300)*2,40+20*(i+snum),51+(j-300)*2,40+20*(i+snum)+19,gray1);
					else gdImageFilledRectangle(image,50+(j-300)*2,40+20*(i+snum),51+(j-300)*2,40+20*(i+snum)+19,white);
				} else {
					if(chars[i][j]!=chars[0][j]) cl=gray2;
					if(chars[i][j]==chars[0][j]) cl=gray1;
					if(chars[i][j]=='-') cl=white;
					gdImageFilledRectangle(image,50+(j-300)*2,40+20*(i+snum),51+(j-300)*2,40+20*(i+snum)+19,cl);
					
				}
				if(!(j % 25)) {
					gdImageLine(image,50+(j-300)*2,40*(snum+1)+2,50+(j-300)*2,40*(snum+1),black);
					sprintf(s,"%i",j);
					gdImageString(image, gdFontGetSmall(),50+(j-300)*2,40*(snum+1)+4,s,black);
				}
			}
			gdImageLine(image,50,40*(snum+1)+2,750,40*(snum+1)+2,black);
			gdImageLine(image,50,40*(snum+1)+3,750,40*(snum+1)+3,black);
		}
		if(slen>650){
			for(j=650;j<1000;j++){
				if(j==slen) break;
				if(i==0){
					if(chars[0][j]!='-') gdImageFilledRectangle(image,50+(j-650)*2,60+20*(i+2*snum),51+(j-650)*2,60+20*(i+2*snum)+19,gray1);
					else gdImageFilledRectangle(image,50+(j-650)*2,60+20*(i+2*snum),51+(j-650)*2,60+20*(i+2*snum)+19,white);
				} else {
					if(chars[i][j]!=chars[0][j]) cl=gray2;
					if(chars[i][j]==chars[0][j]) cl=gray1;
					if(chars[i][j]=='-') cl=white;
					gdImageFilledRectangle(image,50+(j-650)*2,60+20*(i+2*snum),51+(j-650)*2,60+20*(i+2*snum)+19,cl);
					
				}
				if(!(j % 25)) {
					gdImageLine(image,50+(j-650)*2,60*(snum+1)+2,50+(j-650)*2,60*(snum+1),black);
					sprintf(s,"%i",j);
					gdImageString(image, gdFontGetSmall(),50+(j-650)*2,60*(snum+1)+4,s,black);
				}
			}
			gdImageLine(image,50,60*(snum+1)+2,750,60*(snum+1)+2,black);
			gdImageLine(image,50,60*(snum+1)+3,750,60*(snum+1)+3,black);
		}
	}
 

    grout = fopen("alignment.png","wb");
    gdImagePng(image,grout);
    fclose(grout);
    
    gdImageDestroy(image);

    return(0);
}
