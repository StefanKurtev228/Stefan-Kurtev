#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

FILE *fp;

int a,b,i,j,k,l,m,n,k,s,f;

char array[1000];

char line[2000001];

char strs[3000][1000];

double E[3000][42000];

int EI[3000][30000];

int affyN,GSMn,GN,GA[50000],GNR[50000];
char G[50000][100],affys[50000][100];
int lineno;
double Nexp;

double ave,num,sd;

double zcut=2.0;

char files[100][100],GPL[100][100];
int fileN,fN;
double phi,chisq,R;

char filei[100];
char filen[111];
char str[111];

double GSE_I();

double ANALYSE_GSE_I();
double CTidot[4],CTdoti[4],CT[4][4];
int N,M;




int main()
{
	
	
fileN=0;
fp=fopen("GSE/LIST.txt","r");
while(fgets(line,100,fp)!=NULL){
if(strlen(line)>0){
fileN=fileN+1;
sscanf(line,"%s %s",&files[fileN],&GPL[fileN]);	
}}
fclose(fp);

for(fN=1;fN<=fileN;fN++)printf("%d %s %s\n",fN,files[fN],GPL[fN]);

fp=fopen("GSE-I/LIST-TOP-SLICE.txt","w");
fclose(fp);

for(fN=1;fN<=fileN;fN++)GSE_I();


/*
ANALYSE_GSE_I();
*****/

return 0;
}

double GSE_I()
{
	
GSMn=0;

sprintf(filei,"AtoG/AtoG-%s.txt",GPL[fN]); printf("%s\n",filei);
fp=fopen(filei,"r");
affyN=0;
while(fgets(line,100,fp)!=NULL){
if(strlen(line)>0){
affyN=affyN+1; sscanf(line,"%s %d %s",&affys[affyN],&GA[affyN],&G[affyN]);GNR[GA[affyN]]=affyN;if(GA[affyN]>GN)GN=GA[affyN];
}}
fclose(fp);

printf("affyN %d\n",affyN);





for(i=1;i<=affyN;i++){
for(j=1;j<=2999;j++){
E[j][i]=0;
}}



for(i=1;i<=GN;i++){
for(j=1;j<=2999;j++){
EI[j][i]=0;
}}



lineno=-10000000;
sprintf(filei,"GSE/%s.txt",files[fN]);printf("%s\n",filei);
fp=fopen(filei,"r");
while(fgets(line,2000000,fp)!=NULL){


lineno=lineno+1;

if(strlen(line)>100000)goto DROPTHIS;

s=1;j=-1;for(i=0;i<=strlen(line);i++){	
j=j+1;
if((j>=0)&&(line[i]=='\t')){s=s+1;strs[s-1][j]='\0';j=-1;}
if((j<500)&&(j>=0))strs[s][j]=line[i];
if(line[i]=='\n')strs[s][j]='\0';
}


printf("%s %d\n",strs[1],s);


j=-1;for(i=0;i<=strlen(strs[1]);i++){if(strs[1][i]!='\"'){j=j+1;strs[1][j]=strs[1][i];}}
if(strcmp(strs[1],"!series_matrix_table_end")==0){lineno=-1000000;}



if(strcmp(strs[1],"ID_REF")==0){lineno=0;}

if(lineno>0){

m=0;
for(i=1;i<=affyN;i++){if(strcmp(affys[i],strs[1])==0){m=i;goto POP;}}POP:;

if(m!=0){
GSMn=0;for(i=2;i<=s;i++){GSMn=GSMn+1;sscanf(strs[i],"%lf",&E[GSMn][m]);}
}

}

DROPTHIS:;
}

fclose(fp);

printf("GSMn = %d\n",GSMn);


for(k=1;k<=affyN;k++){

num=0;ave=0;sd=0;

for(i=1;i<=GSMn;i++){
num=num+1;
ave=ave+E[i][k];
sd=sd+E[i][k]*E[i][k];
}


ave=ave/num;
if(sd/num-ave*ave>0)sd=sqrt(sd/num-ave*ave);
else sd=0;

if(sd>0){for(i=1;i<=GSMn;i++)E[i][k]=(E[i][k]-ave)/sd;}
else {for(i=1;i<=GSMn;i++)E[i][k]=0;}

for(i=1;i<=GSMn;i++){	
if(E[i][k]>zcut)EI[i][GA[k]]=EI[i][GA[k]]+1;
if(E[i][k]<-zcut)EI[i][GA[k]]=EI[i][GA[k]]-1;
}


}

for(i=1;i<=GN;i++){
for(j=1;j<=GSMn;j++){
k=0;
if(EI[j][i]>0)k=1;
if(EI[j][i]<0)k=2;
EI[j][i]=k;
}}






fp=fopen("GSE-I/LIST-TOP-SLICE.txt","a");
for(j=1;j<=GSMn;j++){
	
k=0;for(i=1;i<=GN;i++){if(EI[j][i]!=0)k=k+1;}

if(k>0){
fprintf(fp,"%d\t",k);
for(i=1;i<=GN;i++){
if(EI[j][i]==1)fprintf(fp,"%d\t",i);
if(EI[j][i]==2)fprintf(fp,"%d\t",-i);}
fprintf(fp,"\n");}
}
fclose(fp);


return 0;
}

