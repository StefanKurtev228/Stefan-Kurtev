#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

FILE *fp;
FILE *fp1;
int a,b,i,j,l,m,n,s,f;

int affyN,GSMn,GN,GA[50000],GNR[50000];
char G[50000][100],affys[50000][100];
int lineno;
int k;
int Gnum;
char Gne[1000][1000];
char array[1000];

char line[2000001];

char strs[3000][1000];
char filei[100];
char filen[111];
char str[111];
char files[100][100],GPL[100][100];
int fileN,fN;
int N,M;
int Gcmp[1000];
char affyS[1000],Gnm[1000];
char Gin[1000];
char pn[1000];
int E[10000][1000];
int pnum;

int PWn[2000],PW[1000][2000];

int main()
{


GN=0;
fp=fopen("AtoG/AtoG-GPL1322.txt","r");
while(fgets(line,100,fp)!=NULL){
if(strlen(line)>0){
sscanf(line,"%s %d %s",&affyS,&k,&Gnm);
if(k>GN)GN=k;
sprintf(G[k],"%s",Gnm);
}}
fclose(fp);




fp1=fopen("PATHWAYSFLY.txt","w");

pnum=0;
fp=fopen("Pathway/c2.cp.v7.0.symbols.gmt","r");
while(fgets(line,2000,fp)!=NULL){
	



s=1;j=-1;for(i=0;i<=strlen(line);i++){	
if(j<500)j=j+1;
if((j>=0)&&(line[i]=='\t')){s=s+1;strs[s-1][j]='\0';j=-1;}
if(j>=0)strs[s][j]=line[i];
if(line[i]=='\n')strs[s][j]='\0';
}


for(i=1;i<=s;i++){if(strlen(strs[i])==0){s=i;goto SSSS;}}SSSS:;


if((strs[1][0]!='K')&&(strs[1][0]!='R'))goto SKIP;

pnum=pnum+1;


PWn[pnum]=0;
for(i=3;i<=s;i++){
for(k=1;k<=GN;k++){
if(strcmp(strs[i],G[k])==0){PWn[pnum]=PWn[pnum]+1; PW[PWn[pnum]][pnum]=k; goto END;}}END:;
}
printf("%d %d\n",s-2,PWn[pnum]);
fprintf(fp1,"%s\t%d\t",strs[1],PWn[pnum]);
for(k=1;k<=PWn[pnum];k++)fprintf(fp1,"%d\t",PW[k][pnum]);
fprintf(fp1,"\n");


SKIP:;

}

fclose(fp);

fclose(fp1);

return 0;
}
