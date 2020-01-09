#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

FILE *fp;

char filei[1000];

char str[1000];

int a,b,i,j,k,l,m,n,s;

int GN;
int T[30001];
int Tp[30001],Tm[30001];
int Tpp[2001][30001],Tpm[2001][30001],Tmp[2001][30001],Tmm[2001][30001];
int BASKET=2000;
int Bn;
int maxG,minG;
double phi;
double CTdotp,CTdotm,CTpdot,CTmdot;
int max1,max2;

int n1,n2,n3,n4;

int main()
{








for(j=1;j<=30000;j++){
Tp[j]=0;Tm[j]=0;
}

GN=0;
fp=fopen("GSE-I/LIST-TOP-SLICEFLY.txt","r");
for(i=1;i<=1000000;i++){
fscanf(fp,"%d",&k);if(feof(fp))goto BBBB;	
for(j=1;j<=k;j++){fscanf(fp,"%d",&n);m=abs(n);if(m>GN)GN=m;
if(n>0)Tp[m]=Tp[m]+1;
if(n<0)Tm[m]=Tm[m]+1;
}
}
BBBB:;
fclose(fp);

fp=fopen("GSE-I/LIST-TpTmFLY.txt","w");
for(i=1;i<=GN;i++){
fprintf(fp,"%d\t%d\t%d\n",i,Tp[i],Tm[i]);
printf("%d %d %d\n",i,Tp[i],Tm[i]);
}
fclose(fp);

for(Bn=1;Bn<=1+GN/BASKET;Bn++){

minG=(Bn-1)*BASKET;
maxG=Bn*BASKET;

for(i=1;i<=BASKET;i++){
for(j=1;j<=GN;j++){
Tpp[i][j]=0;
Tpm[i][j]=0;
Tmp[i][j]=0;
Tmm[i][j]=0;
}
}

max1=0;max2=0;

printf("Genes from %d to %d\n",minG+1,maxG);

fp=fopen("GSE-I/LIST-TOP-SLICEFLY.txt","r");
for(i=1;i<=1000000;i++){
fscanf(fp,"%d",&k);if(feof(fp))goto DDDD;	
for(j=1;j<=k;j++){fscanf(fp,"%d",&T[j]);}
for(n=1;n<=k;n++){
for(m=1;m<=k;m++){
if((abs(T[n])>minG)&&(abs(T[n])<=maxG)){
a=abs(T[n])-minG;
b=abs(T[m]);

if((T[n]>0)&&(T[m]>0))Tpp[a][b]=Tpp[a][b]+1;
if((T[n]>0)&&(T[m]<0))Tpm[a][b]=Tpm[a][b]+1;
if((T[n]<0)&&(T[m]>0))Tmp[a][b]=Tmp[a][b]+1;
if((T[n]<0)&&(T[m]<0))Tmm[a][b]=Tmm[a][b]+1;
}}}

}
DDDD:;
fclose(fp);

sprintf(filei,"GSE-I/LIST-TppTpmTmpTmmFLY-%d.txt",Bn);
fp=fopen(filei,"w");
for(i=1;i<=BASKET;i++){
k=0;
for(j=1;j<=GN;j++){if(Tpp[i][j]+Tpm[i][j]+Tmp[i][j]+Tmm[i][j]>0)k=k+1;}
if(k>0){
fprintf(fp,"%d\t%d\t",minG+i,k);
for(j=1;j<=GN;j++){
if(Tpp[i][j]+Tpm[i][j]+Tmp[i][j]+Tmm[i][j]>0)fprintf(fp,"%d\t%d\t%d\t%d\t%d\t",j,Tpp[i][j],Tpm[i][j],Tmp[i][j],Tmm[i][j]);}
fprintf(fp,"\n");
}}
fclose(fp);

}

















thend:;
return 0;
}
