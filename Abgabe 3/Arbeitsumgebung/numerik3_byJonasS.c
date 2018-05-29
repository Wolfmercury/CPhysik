// gcc -std=c99 -Wpedantic -Wall numerik3.c -o numerik3.exe -lm
// ./numerik3.exe
// Jonas Schneegans, Jan Szyszko



#include <stdio.h>
#include <math.h>


// Z: analytische loesung des integrals ueber rho als funktion der temperatur
// input: t (temperatur)
double Z( double t) {
  // konstanten
  const double h = 6.62606979e-34;        // plancksches wirkungsquantum in J*s
  const double c = 299792458;             // lichtgeschwindigkeit in m/s
  const double pi = 3.14159265359;        // kreiszahl
  const double k = 1.38064852e-23;        // boltzmann konstante in J/K
  // ergebnis
  double erg = (8*pow(k,4)*pow(pi,5)*pow(t,4))/(15*pow(h,3)*pow(c,3));
  return erg;
}

// rho: spektrale energiedichte als funktion der frequenz und
//	der temperatur
// input: v (freqeunz)
//	  t (temperatur)
double rho ( double v, double t) {
  // konstanten
  const double h = 6.62606979e-34;        // plancksches wirkungsquantum in J*s
  const double c = 299792458;             // lichtgeschwindigkeit in m/s
  const double pi = 3.14159265359;        // kreiszahl
  const double k = 1.38064852e-23;        // boltzmann konstante in J/K
  // ergebnis
  double erg = ((8*pi*h)/pow(c,3))*pow(v,3)/(exp((h*v)/(k*t))-1);
  return erg;
}

// sigma: funktion zur bestimmung des wirkungsquerschnittes sigma als
//	  funktion der frequenz
// input: x (frequenz)
double sigma ( double x) {
  // konstanten
  const double pi = 3.14159265359;        // kreiszahl
  const double gamma = 3.0e10;            // konstante aus sigma in 1/s
  const double v2 = 2.0e13;               // ny1, konstante aus sigma in 1/s
  const double v3 = 7.04e13;              // ny2, konstante aus sigma in 1/s
  const double s2 = 2.45e-11;             // konstante aus sigma in m^2/s
  const double s3 = 2.74e-10;             // konstante aus sigma in m^2/s
  // ergebnis
  double erg = (s2*gamma)/(pi*pow((x-v2),2)+pow(gamma,2)) + (s3*gamma)/(pi*pow((x-v3),2)+pow(gamma,2));
  return erg;
}



// integrand: funktion, welche zur bestimmung von epsilon, numerisch
//	      integriert werden soll ( rho*(1-exp(-No*n*sigma)) )
// input: v (frequenz)
//	  t (temperatur)
//	  n (CO2 seaulendichte als anteil von No)
double integrand ( double v, double t, double n) {
  // konstanten
  const double No = 8.1e25;               // aktuelle CO2 sauelendichte in 1/m^2
  // ergebnis
  double erg = rho(v,t)*(1-exp(-No*n*sigma(v)));
  return (1/Z(t))*erg;
}



// simpson: funktion zur numerischen bestimmung des integrals
//          ueber f*rho mit den simpsonverfahren
// input: m (ungerade zahl der zur approximation verwendeten intervalle)
//        a (untere grenze des integrals)
//        b (obere grenze des integrals)
//	  t (temperatur)
//	  n (CO2 seaulendichte als anteil von No)
// output: ergebnis (approximation des integrals)
double simpson ( int m, double a, double b, double t, double n) {
  double sum = 0;               // speichert die zwischenergebnisse
  double hb = fabs(b-a)/m;      // breite eines intervalls
  double xi1 = a;               // ungerade stuetzstellen
  double xi2 = a+hb;            // gerade stuetzstellen
  int i;                        // laufindex
  // summiert die beitraege, die mit 2 multipliziert auf
  for ( i=1; i<m-2; i+=2) {
    xi1 = xi1 + 2*hb;
    sum = sum + 2*integrand(xi1,t,n);
  }
  // summiert die beitraege, die mit 4 multipliziert auf
  for ( i=2; i<m; i+=2) {
    sum = sum + 4*integrand(xi2,t,n);
    xi2 = xi2 + 2*hb;
  }
  sum = sum + integrand(a,t,n) + integrand(b-hb,t,n);
  double ergebnis = (hb/3)*sum;
  return ergebnis;
}

// F: funktion die mit dem fixpunktverfahren geloest werden soll
// input: t (temperatur)
//        n (CO2 seaulendichte als anteil von No)
double F (double t, double n) {
  int m = 100001;
  // konstanten
  const double Ts = 5750;                 // Temperatur der sonnenoberflaeche in K
  const double epsilons = 5.5e-6;         // emmisivitaet der sonne (in diesen model)
  // intervall wurde so geaehlt, dass der integrand außerhalb
  // naeherrungsweise verschwindet
  double a = 1;
  double b = 1e16;
  double erg = Ts*sqrt(sqrt((epsilons*(2-simpson(m,a,b,Ts,n)))/(2-simpson(m,a,b,t,n))));
  return erg;
}

// fixpunkt: funktion zum lösen der betrachten gleichung
//           mithilfe des fixpunktverfahrens
// input: start (startpunkt)
//        epsilon (genauigkeit)
//	  n (CO2 saeulendichte als anteil von No)
// output: ergebnis (mittlere erdtemperatur als funktion von n
//	   und in kelvin)
double fixpunkt ( double n, double start, double epsilon) {
  double sum1 = start;		// speichert das ergebnis von i
  double sum2 = start;          // speichert das ergebnis von i-1
  int i;
  double ergebnis;
  for( i=0; ; i++) {
    sum1 = F(sum1,n);
    // wenn die genauigkeit erreicht ist wird sum1 = ergebnis gesetzt
    if ( fabs(sum1-sum2) < epsilon) {
      ergebnis = sum1;
      break;
    }
    // wenn die genauigkeit nicht erreicht ist wird der wert in sum2 gespeichert
    sum2 = sum1;
  }
  return ergebnis;
}


// dreipunkt: funktion zur berechnung der ableitung mit der 3-punktformel
// input: x (stelle, an der die ableitung berechnet wird)
//        h (schrittbreite)
// output: y (wert der ableitung an der stelle x)
double dreipunkt(double x, double h) {
  double n = 1;			// CO2 saeulendichte als anteil von No
  double y = (F(x+h,n) - F(x-h,n))/(2*h);
  return y;
}



// main-funkton ruft nur noch fixpunkt auf und gibt den wert
// fuer n=1 aus (in grad celsius) und prueft ob die ableitung
// bei TE konvergieren kann (abgeleitet wird nach der Dreipunkformel)
int main(void) {
  double ergebnis = fixpunkt(1,300,1e-6);
  double abl = dreipunkt( ergebnis, 1e-3);	// ableitung von bei ergebnis
  printf("Die mittlere Erdtemperatur in grad Celsius und mit n=1\nlautet nach diesem Modell:\n%f\n", ergebnis-273.15);
  printf("Die Ableitung von F bei dieser Temperatur betraegt ca:\n%e\n", abl);
  return 0;
}


