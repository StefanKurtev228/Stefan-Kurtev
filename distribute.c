#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
char line[1000];
FILE *fp;
FILE *fp1;
int k,i,j,n,m,x,r,t,q,Tpp,Tpm,Tmp,Tmm,Pn;
double phi,phiran;
int GN;
int CTdotp,CTdotm,CTpdot,CTmdot;
int Zn[10000],phin[10000];
char str[1000][10000];
double phisum[10000],Zsum[10000],phiavg[10000],Zavg[10000];
int PWn,PWN[3000],PW[3000][25000];
char filei[1000];
char filen[1000];
int Bn;
int B;
int BASKET=10;
double Zs,Zsran;
int D[1000],DK,DKran,Dran[1000];

int PWrand[3000][25000];
int CTdotpran,CTdotmran,CTpdotran,CTmdotran;

double PLOT();
double X[1000];

int main()

{

PLOT();goto thend;	
	
PWn=0;GN=0;
fp=fopen("PATHWAYS.txt","r");
for(i=1;i<=1000;i++){
fscanf(fp,"%s %d",&str[i],&k);if(feof(fp))goto SKIP;
PWn=PWn+1;
PWN[PWn]=k;
for(j=1;j<=k;j++){fscanf(fp,"%d",&n);PW[PWn][n]=1;if(n>GN)GN=n;}
}
SKIP:;
fclose(fp);

for(q=1;q<=PWn;q++){
Zsum[q]=0;
phisum[q]=0;
phin[q]=0;
Zn[q]=0;
}

srand((unsigned int)time(NULL));

for(i=1;i<=PWn;i++){
for(j=1;j<=PWN[i];j++){
k=1+rand()%GN;
PWrand[i][k]=1;
}}



for(Bn=1;Bn<=12;Bn++){

printf("%d\n",Bn);
sprintf(filei,"GSE-I/LIST-TppTpmTmpTmm-%d.txt",Bn);
fp=fopen(filei,"r");
for(i=1;i<=100000;i++){
fscanf(fp,"%d %d",&n,&k);if(feof(fp))goto WWWW;
for(j=1;j<=k;j++){
fscanf(fp,"%d %d %d %d %d",&m,&Tpp,&Tpm,&Tmp,&Tmm);



for(q=1;q<=PWn;q++){

if((PWrand[q][n]*PWrand[q][m]==1)&&(n>m)){
CTdotpran=Tpp+Tmp;
CTdotmran=Tpm+Tmm;
CTpdotran=Tpp+Tpm;
CTmdotran=Tmp+Tmm;

phiran=0;if(CTdotpran*CTdotmran*CTpdotran*CTmdotran>0)phiran=(Tpp*Tmm-Tpm*Tmp)/sqrt(CTdotpran*CTdotmran*CTpdotran*CTmdotran);
Zsran=phiran*sqrt(Tpp+Tpm+Tmp+Tmm);
DKran=(int)(fabs(Zsran)*10);
Dran[DKran]=Dran[DKran]+1;
}
if((PW[q][n]*PW[q][m]==1)&&(n>m)){


CTdotp=Tpp+Tmp;
CTdotm=Tpm+Tmm;
CTpdot=Tpp+Tpm;
CTmdot=Tmp+Tmm;

phi=0;if(CTdotp*CTdotm*CTpdot*CTmdot>0)phi=(Tpp*Tmm-Tpm*Tmp)/sqrt(CTdotp*CTdotm*CTpdot*CTmdot);
phisum[q]=phisum[q]+phi;
phin[q]=phin[q]+1;
Zs=phi*sqrt(Tpp+Tpm+Tmp+Tmm);
Zsum[q]=Zsum[q]+Zs;
Zn[q]=Zn[q]+1;
DK=(int)(fabs(Zs)*10);
D[DK]=D[DK]+1;
}
}

}}

WWWW:;
fclose(fp);
}

fp=fopen("DIS.txt","w");
for(i=0;i<=110;i++){
fprintf(fp,"%f\t%d\n",0.1*i,D[i]);
printf("%f\t%d\n",0.1*i,D[i]);
}
fclose(fp);


fp=fopen("DISran.txt","w");
for(i=0;i<=110;i++){
fprintf(fp,"%f\t%d\n",0.1*i,Dran[i]);
printf("%f\t%d\n",0.1*i,Dran[i]);
}
fclose(fp);






thend:;


return 0;
}



double PLOT()
{
int N;
N=0;
fp=fopen("FLYDIS.txt","r");
while(fgets(line,100,fp)!=NULL){
N=N+1;
sscanf(line,"%lf %d",&X[N],&D[N]);
}
fclose(fp);

N=0;
fp=fopen("FLYDISran.txt","r");
while(fgets(line,100,fp)!=NULL){
N=N+1;
sscanf(line,"%lf %d",&X[N],&Dran[N]);
}
fclose(fp);

fp=fopen("plot.txt","w");
for(i=3;i<=N-1;i++)fprintf(fp,"%f\t%d\t%d\n",X[i],(D[i-1]+D[i]+D[i+1])/3,(Dran[i-1]+Dran[i]+Dran[i+1])/3);
fclose(fp);

fp=fopen("plot.plt","w");
fprintf(fp,"plot \"plot.txt\" using 1:2 with lines, \"plot.txt\" using 1:3 with lines\n");
fclose(fp);

system("\"\"C:/Program Files/gnuplot/bin/wgnuplot\" -persist \"plot.plt\"\"");

return 0;

}




