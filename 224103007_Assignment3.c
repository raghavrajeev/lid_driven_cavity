#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main(){

//Dimensions of cavity
double H=1.0,L=1.0,Re;
int m,n;

//Input for Re, m and n
printf("Enter Re: ");
scanf("%lf",&Re);
printf("\n");

printf("Enter M (horizontal grid size): ");
scanf("%d",&m);
printf("\n");

printf("Enter M (vertical grid size): ");
scanf("%d",&n);
printf("\n");

double psi1[n][m],w1[n][m],psi2[n][m],w2[n][m],v[n][m],u[n][m],beta,beta_sq,dx=H/(m-1.0),dy=H/(n-1.0);
beta=dx/dy;
beta_sq=beta*beta;

//Boundary and initial conditions for velocity and streamfunction
for(int i=0;i<m;i++){
    for(int j=0;j<n;j++){
        if(i==0){
            u[j][i]=0.0;
            v[j][i]=0.0;
            psi2[j][i]=0.0;
        }
        else if(i==m-1){
            u[j][i]=0.0;
            v[j][i]=0.0;
            psi2[j][i]=0.0;
        }
        else if(j==0){
            u[j][i]=0.0;
            v[j][i]=0.0;
            psi2[j][i]=0.0;
        }
        else if(j==n-1){
            u[j][i]=1.0;
            v[j][i]=0.0;
            psi2[j][i]=0.0;
        }
        else{
            u[j][i]=0.0;
            v[j][i]=0.0;
            psi2[j][i]=0.0;
        }
    }
}

//Initial value for vorticity at boundary
for(int i=0;i<m;i++){
    for(int j=0;j<n;j++){
        if(j==0){
            w2[j][i]=(2.0/(dy*dy))*(psi2[j][i]-psi2[j+1][i]);}
        else if(j==n-1){
            w2[j][i]=(2.0/(dy*dy)*(psi2[j][i]-psi2[j-1][i]))-(2.0*u[j][i]/dy);}
        else if(i==0){
            w2[j][i]=(2.0/(dx*dx))*(psi2[j][i]-psi2[j][i+1]);}
        else if(i==m-1){
            w2[j][i]=(2.0/(dx*dx))*(psi2[j][i]-psi2[j][i-1]);}
        else{
            w2[j][i]=0.0;}
    }
}


double error_psi=1.0,error_w=1.0,iter=0;


while(error_psi>1e-6||error_w>1e-6){
    //Copying of matrix from previous iteration
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            w1[j][i]=w2[j][i];
            psi1[j][i]=psi2[j][i];
        }
    }

    //Solving for streamfunction
    for(int j=1;j<m-1;j++){
        for(int i=1;i<n-1;i++){
            psi2[j][i]=(1.0/(2.0*(1+beta_sq)))*(psi2[j][i+1]+psi2[j][i-1]+beta_sq*(psi2[j+1][i]+psi2[j-1][i])+dx*dx*w2[j][i]);
        }
    }

    //Solving for omega
    for(int j=1;j<m-1;j++){
        for(int i=1;i<n-1;i++){
            w2[j][i]=(1.0/(2.0*(1+beta_sq)))*(     (1.0-(psi2[j+1][i]-psi2[j-1][i])*beta*Re/4)*w2[j][i+1]  +  (1.0+(psi2[j+1][i]-psi2[j-1][i])*beta*Re/4)*w2[j][i-1]
            +    (1.0+(psi2[j][i+1]-psi2[j][i-1])*Re/(4*beta))*beta_sq*w2[j+1][i]   +    (1.0-(psi2[j][i+1]-psi2[j][i-1])*Re/(4*beta))*beta_sq*w2[j-1][i]    );
        }
    }

    //Omega at boundary
    for(int j=0;j<m;j++){
        for(int i=0;i<n;i++){
            if(j==0){
                w2[j][i]=(2.0/(dy*dy))*(psi2[j][i]-psi2[j+1][i]);}
            else if(j==n-1){
                w2[j][i]=((2.0/(dy*dy))*(psi2[j][i]-psi2[j-1][i]))-(2.0*u[j][i]/dy);}
            else if(i==0){
                w2[j][i]=(2.0/(dx*dx))*(psi2[j][i]-psi2[j][i+1]);}
            else if(i==m-1){
                w2[j][i]=(2.0/(dx*dx))*(psi2[j][i]-psi2[j][i-1]);}
        }
    }

    error_psi=0;
    error_w=0;

    //Error calculation
    for(int i=1;i<m-1;i++){
        for(int j=1;j<n-1;j++){
            error_psi=error_psi+pow(psi2[j][i]-psi1[j][i],2);
            error_w=error_w+pow(w2[j][i]-w1[j][i],2);
        }
    }

    error_psi=sqrt(error_psi/((m-2)*(n-2)));
    error_w=sqrt(error_w/((m-2)*(n-2)));

    iter+=1;

    printf("Iterations=%.lf\tpsi_error=%.9lf\tomega_error=%.9lf\n",iter,error_psi,error_w);

}

//Calculating velocity inside the cavity
for(int i=1;i<m-1;i++){
    for(int j=1;j<n-1;j++){
        u[j][i]=(1/(2*dy))*(psi2[j+1][i]-psi2[j-1][i]);
        v[j][i]=(-1/(2*dx))*(psi2[j][i+1]-psi2[j][i-1]);
    }
}


//Printing output to files
FILE *fp;

fp=fopen("streamfunction.dat","w");
for(int i=0;i<m;i++){
    for(int j=0;j<n;j++){
        fprintf(fp,"%lf\t%lf\t%lf\n",i*dx,j*dy,psi2[j][i]);
    }
}
fclose(fp);


fp=fopen("velocity.dat","w");
for(int i=0;i<m;i++){
    for(int j=0;j<n;j++){
        fprintf(fp,"%lf\t%lf\t%lf\t%lf\n",i*dx,j*dy,u[j][i],v[j][i]);
    }
}
fclose(fp);

fp=fopen("u_cent.dat","w");
int mid_x=m/2;
for(int j=0;j<n;j++){
    fprintf(fp,"%lf\t%lf\n",u[j][mid_x],j*dy);
}
fclose(fp);

int mid_y=n/2;
fp=fopen("v_cent.dat","w");
for(int i=0;i<m;i++){
    fprintf(fp,"%lf\t%lf\n",v[mid_y][i],i*dx);
}
fclose(fp);

return 0;
}
