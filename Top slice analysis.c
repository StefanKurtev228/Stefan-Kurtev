
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

FILE *fp;
FILE *fp1;

char filei[1000];
char fileo[1000];

char str[1000];

int a,b,i,j,k,l,m,n,s;

int GN;
int Tpp,Tpm,Tmp,Tmm;
int Bn;
int maxG,minG;
double phi;
double CTdotp,CTdotm,CTpdot,CTmdot;
int max1,max2;

int n1,n2,n3,n4;

int N;


int Zi[30000];
double ZZ[30000];

double AVE,NUM;



int main()
{



fp1=fopen("GSE-I/ZSCORE.txt","w");

N=0;
for(Bn=1;Bn<=12;Bn++){
sprintf(filei,"GSE-I/LIST-TppTpmTmpTmm-%d.txt",Bn);
fp=fopen(filei,"r");
for(i=1;i<=1000000;i++){
fscanf(fp,"%d %d",&n,&k);if(feof(fp))goto WWWW;



for(j=1;j<=k;j++){
fscanf(fp,"%d %d %d %d %d",&m,&Tpp,&Tpm,&Tmp,&Tmm);

CTdotp=Tpp+Tmp;
CTdotm=Tpm+Tmm;
CTpdot=Tpp+Tpm;
CTmdot=Tmp+Tmm;

phi=0;
if(CTdotp*CTdotm*CTpdot*CTmdot>0){
phi=(Tpp*Tmm-Tpm*Tmp)/sqrt(CTdotp*CTdotm*CTpdot*CTmdot);N=N+1;
}

printf("%d %d %f %f\n",n,m,phi*sqrt(Tpp+Tpm+Tmp+Tmm),phi);
}
}




}
WWWW:;
fclose(fp);

printf("%d\n",Bn);

}
fprintf(fp1,"%d\n",N);
fclose(fp1);




thend:;
return 0;
}
