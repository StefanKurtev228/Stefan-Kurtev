#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

FILE *fp;

int N,i,j,t,k,n,m,GN,GNF,s,lineno;

char line[100];
char strs[3000][1000];



char DRsym[30000][1000],Hgsym[30000][1000];
char G[30000][1000],GF[30000][1000];
char c;
int AA[1000];

int HOMhfN[25000],HOMhf[100][25000];
int HOMfhN[25000],HOMfh[100][25000];

int HOMhfsc[100][25000];
int HOMfhsc[100][25000];

int DIOPT;

int main()
{

for(i=0;i<=255;i++){
AA[i]=i;}

for(i=97;i<=122;i++){
AA[i]=i-32;}

for(i=0;i<=255;i++){c=(char)(i);
if((int)c>=0)k=AA[(int)c];
}

fp=fopen("AtoG/AtoG-GPL570.txt","r");
GN=0;
while(fgets(line,100,fp)!=NULL){
if (strlen(line)>0){
sscanf(line, "%s %d %s", &strs[1],&k,&strs[2]); 
sprintf(G[k],"%s",strs[2]); if (k>GN) GN=k;
}
}
fclose(fp);

for(i=1;i<=GN;i++)HOMhfN[i]=0;

fp=fopen("AtoG/AtoG-GPL1322.txt","r");
GNF=0;
while(fgets(line,100,fp)!=NULL){
if (strlen(line)>0){
sscanf(line, "%s %d %s", &strs[1],&t,&strs[2]); 
sprintf(GF[t],"%s",strs[2]); if (t>GNF) GNF=t;
}
}
fclose(fp);

for(i=1;i<=GNF;i++)HOMfhN[i]=0;

N=0;
lineno=-10000;
fp=fopen("AtoG/dmel_human_orthologs_disease_fb_2019_04.tsv", "r");
while(fgets(line,10000,fp)!=NULL){
lineno=lineno+1;

s=1;j=-1;for(i=0;i<=strlen(line);i++){	
j=j+1;
if((j>=0)&&(line[i]=='\t')){s=s+1;strs[s-1][j]='\0';j=-1;}
if(j>=0)strs[s][j]=line[i];
if(line[i]=='\n')strs[s][j]='\0';
}



j=-1;for(i=0;i<=strlen(strs[1]);i++){if(strs[1][i]!='\"'){j=j+1;strs[1][j]=strs[1][i];}}

if(strcmp(strs[1],"Finished report_human_orthologs: Thu Aug  8 20:41:39 2019" )==0) {lineno=-1000000;}

if(strcmp(strs[1],"##Dmel_gene_ID")==0){lineno=0;}

if(lineno>0){
	
for(m=0;m<=strlen(strs[2]);m++){k=(int)strs[2][m];if(k>=0){k=AA[k];strs[2][m]=(char)k;}}
for(m=0;m<=strlen(strs[5]);m++){k=(int)strs[5][m];if(k>=0){k=AA[k];strs[5][m]=(char)k;}}


n=0;for(i=1;i<=GN;i++){if(strcmp(G[i],strs[5])==0){n=i;goto QQQQ;}}QQQQ:;

m=0;for(i=1;i<=GNF;i++){if(strcmp(GF[i],strs[2])==0){m=i;goto QQQQ1;}}QQQQ1:;


if(n*m!=0){
HOMhfN[n]=HOMhfN[n]+1;HOMhf[HOMhfN[n]][n]=m;
HOMfhN[m]=HOMfhN[m]+1;HOMfh[HOMfhN[m]][m]=n;
}



}}
fclose(fp);

fp=fopen("GSE-I/OTHOLOG_HtoFLY.txt","w");
for(i=1;i<=GN;i++){if(HOMhfN[i]>0){
fprintf(fp,"%d\t%d\t",i,HOMhfN[i]);
for(j=1;j<=HOMhfN[i];j++)fprintf(fp,"%d\t",HOMhf[j][i]);fprintf(fp,"\n");
}}
printf("%d\n",GNF);
fp=fopen("GSE-I/OTHOLOG_FLYtoH.txt","w");
for(i=1;i<=GNF;i++){if(HOMfhN[i]>0){
fprintf(fp,"%d\t%d\t",i,HOMfhN[i]);
for(j=1;j<=HOMfhN[i];j++)fprintf(fp,"%d\t",HOMfh[j][i]);fprintf(fp,"\n");
}}

return 0;

}


















