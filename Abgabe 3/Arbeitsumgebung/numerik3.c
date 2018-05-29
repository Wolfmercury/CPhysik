#include <stdio.h>
#include <stdlib.h>
#include <math.h>
double planck = 6.626070040*pow(10,-34);
double c =  299792458.0;
double k = 1.38064852*pow(10,-23);
double nu_2 = 2*pow(10,13);
double nu_3 = 7.04*pow(10,13);
double Gamma = 3*pow(10,10);
double S_2 = 2.45*pow(10,-11);
double S_3 = 2.74*pow(10,-10);
double epsilon_S = 5.5*pow(10,-6);
double T_S = 5750.0;
double N_0 = 8.1*pow(10,25);

double rho(double ny, double T){
    return 8*M_PI*planck*pow(ny,3.0)/(pow(c,3.0)*(exp(ny/T)-1));
}

double sigma(double ny){
	return (S_2*Gamma/(M_PI*pow((ny-nu_2),2)+pow(Gamma,2)))+(S_3*Gamma/(M_PI*pow((ny-nu_3),2)+pow(Gamma,2)));
}
double f(double ny,double eta){
	return exp(-N_0*eta*sigma(ny));
}
double Z_I(double T,double Genauigkeit){
	double Z = 1.0;
    double Wert = 0.0;
    double m = 1;
    while ((Z-Wert)/Z > Genauigkeit){
		if (m == 1){
			Z = 0;
		}
        Wert = Z;
        Z += (48*M_PI*pow(k,4)*pow(T,4))/(pow(c,3)*pow(planck,3)*pow(m,4));
        m += 1;
	}
	return 1/Z;
}
double Simpson(nu,eta,T){
    return rho(nu,T)*(1-f(eta,nu));
}
double Integral(double l, double n,double T, double eta){
    double h=l/n;
    int i=1;
    double H;
    double sum=0.0;
    while (i<(n/h)){
            H = i*h;
            if(i%2==0){
                sum=sum+4*Simpson(H,eta,T);
            }else{
                sum=sum+2*Simpson(H,eta,T);
            }
            i=i+1;
    }
    return sum*(h/3);
}

int main(){
    return 0;
}
