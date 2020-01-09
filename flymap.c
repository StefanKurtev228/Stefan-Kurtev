#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

FILE *fp;
FILE *fp1;
int i,j;
int Fg,k,n,t,m;
char Pn[1000];
char line[100];

int FH[30000];

int PP[10000];


int main()
{

fp=fopen("GSE-I/OTHOLOG_FLYtoH.txt","r");
for(i=1;i<=10000000;i++){
fscanf(fp,"%d",&Fg);if(feof(fp))goto UUUU;
fscanf(fp,"%d",&k);if(feof(fp))goto UUUU;
for(j=1;j<=k;j++){fscanf(fp,"%d",&n);FH[n]=Fg;}
printf("%d\n",Fg);
}
UUUU:;
fclose(fp);

fp1=fopen("FLYPATHS.txt","w");
fp=fopen("PATHWAYS.txt","r");
for(i=1;i<=100000;i++){
fscanf(fp,"%s",&Pn);if(feof(fp))goto END;
fscanf(fp,"%d",&t);if(feof(fp))goto END;

k=0;
for(j=1;j<=t;j++){fscanf(fp,"%d",&n);if(FH[n]!=0){
for(m=1;m<=k;m++){if(FH[n]==PP[m])goto DROPit;}
k=k+1;PP[k]=FH[n];
DROPit:;
}}

fprintf(fp1,"%s\t%d\t",Pn,k);

for(j=1;j<=k;j++)fprintf(fp1,"%d\t",PP[j]);

fprintf(fp1,"\n");
}
END:;
fclose(fp);







return 0;
}
