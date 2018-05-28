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

double x_i(double k, double i){
    return k*i;
}
double rho(double ny, double T){
    return 8*M_PI*planck*pow(ny,3.0)/(pow(c,3.0)*(exp(ny/T)-1));
}
double rho_Abl(double ny, double T){
    return 8*M_PI*planck*pow(ny,3.0)*planck*ny*exp(planck*ny/(k*T))/(pow(c,3)*k*T*T*pow((exp(planck*ny/(k*T))-1),2));
}
double sigma(double ny){
	return (S_2*Gamma/(M_PI*pow((ny-nu_2),2)+pow(Gamma,2)))+(S_3*Gamma/(M_PI*pow((ny-nu_3),2)+pow(Gamma,2)));
}
double f(double ny){
	return exp(-N_0*sigma(ny));
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
double Z_I_Abl(double T,double Genauigkeit){
	double sum = 1.0;
    double Wert = 0.0;
    double m = 1;
    while ((sum-Wert)/sum > Genauigkeit){
		if (m == 1){
			sum = 0;
		}
        Wert = sum;
        sum += (48*M_PI*pow(k,4))/(pow(c,3)*pow(planck,3)*pow(m,4));
        m += 1;
	}
	return -1/(pow(T,5)*sum);
}
double Integral(double l, double n,double T){
    double h=l/n;
    double i=0.0;
    double sum=0.0;
    while (x_i(h, i)<n){
        sum=sum+h*rho(x_i(h,i),T)*(1-f(x_i(h,i)));
        i=i+1.0;
    }
    return sum;
}
double Integral_Abl(double l, double n,double T){
    double h=l/n;
    double i=0.0;
    double sum=0.0;
    while (x_i(h, i)<n){
        sum=sum+h*rho_Abl(x_i(h,i),T)*(1-f(x_i(h,i)));
        i=i+1.0;
    }
    return sum;
}
double Newton (double eps, double l, double n, double Tstart,double Genauigkeit){
    double dT;
    double T=Tstart;
    while (dT<eps){
        dT=Integral(l,n,T)*Z_I(T,Genauigkeit)/((Integral_Abl(l,n,T)*Z_I(T,Genauigkeit))+Z_I_Abl(T,Genauigkeit)*Integral(l,n,T));
        T=T-dT;
    }
    return T;
}
int main(){
    return 0;
}
