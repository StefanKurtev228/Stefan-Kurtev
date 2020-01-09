#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

FILE *fp;
FILE *fp1;
int R,f,h,k,n,m,i,j;
int Hpn,Qn,Fpn;
char line[100];
int Hn[30000],Fn[100][30000];
int x;
double Z,Zp;

int PAIRn;
int PAIRi[1000000],PAIRj[1000000];
double PAIRsc[1000000];
double Cpn,Cnp,Cnn,Cpdot,Cdotp,Cndot,Cdotn;
int Pa,Pb;
int N;
double Nf=13796;
double Neff;
double Zscore,phi;

double prob;


int main()
{

fp=fopen("GSE-I/OTHOLOG_HtoFLY.txt","r");
for(n=1;n<=1000000;n++){
fscanf(fp,"%d",&n);if(feof(fp))goto XXXX1;
fscanf(fp,"%d",&m);if(feof(fp))goto XXXX1;
Hn[n]=m;
for(i=1;i<=Hn[n];i++){
fscanf(fp,"%d",&Fn[i][n]);
}
}
XXXX1:;
fclose(fp);




PAIRn=0;
fp=fopen("GSE-I/ZSCORE.txt","r");
while(fgets(line,100,fp)!=NULL){
if(strlen(line)>0){
sscanf(line,"%d %d %lf",&n,&m,&Z);if(m>n)printf("mistake\n");
for(i=1;i<=Hn[n];i++){
for(j=1;j<=Hn[m];j++){
if(Fn[i][n]>Fn[j][m]){Pb=Fn[j][m];Pa=Fn[i][n];}
else {Pa=Fn[j][m];Pb=Fn[i][n];}

for(k=1;k<=PAIRn;k++){
if((PAIRi[k]==Pa)&&(PAIRj[k]==Pb)&&(PAIRsc[k]>Z)){PAIRsc[k]=Z;goto SKIP;}
}
PAIRn=PAIRn+1;PAIRi[PAIRn]=Pa;PAIRj[PAIRn]=Pb;PAIRsc[PAIRn]=Z;
SKIP:;
}
}
}
}
fclose(fp);

printf("number of pairs %d\n",PAIRn);



fp=fopen("FlYP.txt","w");
for(i=1;i<=PAIRn;i++){
fprintf(fp,"%d\t%d\t%f\n",PAIRi[i],PAIRj[i],PAIRsc[i]);
}
fclose(fp);

Qn=0;
Fpn=0;
fp1=fopen("MATCHEDPAIRS.txt","w");
fp=fopen("GSE-I/ZSCOREFLY.txt","r");
while(fgets(line,100,fp)!=NULL){
if(strlen(line)>0){
sscanf(line,"%d %d %lf",&i,&j,&Zp);if(j>i)printf("mistake\n");
Fpn=Fpn+1;
for(x=1;x<=PAIRn;x++){
if((PAIRi[x]==i)&&(PAIRj[x]==j)){Qn=Qn+1;
fprintf(fp1,"%d\t%d\t%lf\n",i,j,Zp);
printf("%d\t%d\t%lf\n",i,j,Zp);
}

}
}
}
fclose(fp);
fclose(fp1);

printf("pairs matched %d %d %d\n",Qn,PAIRn,Fpn);

Neff=Nf*Nf/2;printf("%f\n",Neff);

prob=(double)Fpn/Neff;


printf("probability of getting a matched pair is %f expected number is %f\n",prob,(double)PAIRn*prob);

printf("standard deviation %f\n",sqrt(Neff*prob*(1-prob)));

printf("n-sigma = %f\n", (double)(Qn-(double)PAIRn*prob)/sqrt((double)PAIRn*prob*(1-prob)));



thend:;

return 0;

}








