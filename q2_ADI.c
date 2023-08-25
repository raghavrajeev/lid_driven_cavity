#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main(){
int m=31,n=31;//Grid
double xl=1.0,yl=1.0,dx,dy,beta,beta_sq;
double phi1[n][m],phi2[n][m];
dx=xl/(m-1);//step size
dy=yl/(n-1);
beta=dx/dy;
beta_sq=beta*beta;

//Boundary conditions
for (int i=0;i<m;i++){
    for(int j=0;j<n;j++){
        if(j==n-1){//top boundary
            phi2[j][i]=0.0;
        }
        else if(j==0){//bottom boundary
            phi2[j][i]=1.0;
        }
        else if(i==0){//left boundary
            phi2[j][i]=1.0;
        }
        else if(i==m-1){//right boundary
            phi2[j][i]=1.0;
        }
        else{
            phi2[j][i]=0.5;
        }
    }
}

FILE *fp1,*fp2;
fp1=fopen("q2_ADI.dat","w");//For storing solution of phi
fp2=fopen("q2_ADI_error.dat","w");//For storing error at each iteration
fprintf(fp1,"x \t y \t phi \n");
fprintf(fp2,"Iteration \t Error \n");

//Iterative solving
double A,B[m],C[m],D[m],P[m],Q[m],error=1.0,iter=0,A_2,B_2[n],C_2[n],D_2[n],P_2[n],Q_2[n];

while(error>1e-6){
    for (int i=0;i<m;i++){
        for (int j=0;j<n;j++){
            phi1[j][i]=phi2[j][i];//Copying of matrix from previous iteration
        }
    }
    for(int i=1;i<m-1;i++){

        A_2=-2*(1+beta_sq);
        B_2[0]=1;
        C_2[0]=1;
        D_2[0]=-beta_sq*(phi2[0][i-1]+phi1[0][i+1]);
        P_2[0]=-B_2[0]/A_2;//Assigning P and Q at bottom most point
        Q_2[0]=(D_2[0]-C_2[0])/A_2;

        for(int j=1;j<n-1;j++){//Assigning P and Q at each point in line
            D_2[j]=-beta_sq*(phi2[j][i-1]+phi1[j][i+1]);
            B_2[j]=1;
            C_2[j]=1;
            P_2[j]=-B_2[j]/(A_2+C_2[j]*P_2[j-1]);
            Q_2[j]=(D_2[j]-C_2[j]*Q_2[j-1])/(A_2+C_2[j]*P_2[j-1]);
        }

        for(int j=n-2;j>=1;j--){
            phi2[j][i]=P_2[j]*phi2[j+1][i]+Q_2[j];//Back solving phi
        }

    }
    for(int j=1;j<n-1;j++){
        A=-2*(1+beta_sq);
        B[0]=1;
        C[0]=1;
        D[0]=-beta_sq*(phi2[j-1][0]+phi2[j+1][0]);
        P[0]=-B[0]/(A);//Assigning P and Q at left most point
        Q[0]=(D[0]-C[0])/(A);
        for(int i=1;i<m-1;i++){//Assigning P and Q at each point in line
            D[i]=-beta_sq*(phi2[j-1][i]+phi2[j+1][i]);
            B[i]=1;
            C[i]=1;
            P[i]=-B[i]/(A+C[i]*P[i-1]);
            Q[i]=(D[i]-C[i]*Q[i-1])/(A+C[i]*P[i-1]);
        }

        for(int i=m-2;i>=1;i--){
            phi2[j][i]=P[i]*phi2[j][i+1]+Q[i];//Back solving phi
        }
    }

    error=0;

    for (int i=0;i<m;i++){
        for (int j=0;j<n;j++){
            error+=pow((phi2[j][i]-phi1[j][i]),2);//Error calculation
        }
    }
    iter=iter+1;
    error=sqrt(error/((m-2)*(n-2)));//Calculation of error
    printf("Iteration %.1lf  ",iter);
    printf("Error %.8lf\n",error);
    fprintf(fp2,"%.1lf\t%.8lf\n",iter,error);//Storing error at each iteration

}

for (int i=0;i<m;i++){
    for (int j=0;j<n;j++){
        fprintf(fp1,"%.2f \t %.2f \t %lf \n",i*dx,j*dy,phi2[j][i]);//Storing solution at each grid point
    }
}
fclose(fp1);
fclose(fp2);
return 0;
}
