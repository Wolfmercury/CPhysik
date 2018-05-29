# coding: utf8
'''
// python3 numerik3.py
// Kilian Kranz, Achim Bruchertseifer
'''
import numpy as np

#Konstanten:
planck = 6.626070040e-34    #planckkonstante
c =  299792458              #lichtgeschwindigkeit [m/s]
k = 1.38064852e-23          #Boltzmann const [J K^-1]
nu_2 = 2e13                 #von sigma  [1/s]
nu_3 = 7.04e13              #von sigma  [1/s]
Gamma = 3e10                #Gamma von sigma [1/s]
S_2 = 2.45e-11              #austrahlung dings [m^2/s]
S_3 = 2.74e-10              #austrahlung dingsn [m^2/s]
epsilon_S = 5.5e-6          #emessivitaet
T_S = 5750                  #Temperatur der Sonne [K]
N_0 = 8.1e25                #Anzahl an momentanen CO2 molekuelen in luftsaule mit 1 m^2 querschnitt in [m^-2]


class treibhaus:
    #input: Temperatur T, Wellenlaenge nu
    #output: rho
    #berechnung der spektralen Energiedichte fuer unpolarisierte Strahlung
    def rho(nu,T):
        Rechnung = ((8*np.pi*planck)/(c**3))*nu**3/(np.exp((planck*nu)/(k*T))-1)
        #print(nu)
        return Rechnung   
    
    #input: Frequent nu
    #output: Wirkungsquerschnitt sigma
    #Formel für Wirkungsquerschnitt
    def sigma(nu):
        Wert = (S_2/np.pi)*Gamma/((nu-nu_2)**2+Gamma**2) + (S_3/np.pi)*Gamma/((nu-nu_3)**2+Gamma**2)
       # print(nu,Wert)
        return Wert
    
    #input: n_schlange (hier eta) und nu
    #output: Wahscheinlichkeit fuer das Entweichen eines Photons der Frequenz nu (f)
    #funktion im Integral    
    def f(eta,nu):
        Wert = np.exp(-(eta*N_0*treibhaus.sigma(nu)))
        return Wert
        
    #Für das Integral nach Simpsons benötigte hilfsfunktion
    def Simpson_Fkt(nu,eta,T):
        Wert = treibhaus.rho(nu,T)*(1-treibhaus.f(eta,nu))
        #print(treibhaus.rho(nu,T))
        return Wert

    #input; l: Breite des Intervalls, n: Anzahl Stützstellen (Gerade), T: Temperatur, eta: n schlange
    #output: Integral
    #numerische Berechnung des Integrals nach Simpson
    def Integral_S(l,n,T,eta):
        Summe = 0   #Gesammte Summe
        h = abs(l-1)/n  #breite eines intervalls
        xi1 = 1         #für die ungeraden Intervallschritte
        xi2 = 1+h       #für die geraden Intervallschritte
        u = 0           #laufindex
        while u<n:
            if u % 2 == 0:  #gerade komponenten
                xi1 += 2*h
                Summe += 2*treibhaus.Simpson_Fkt(xi1,eta,T)
            else:           #ungerade
                Summe += 4*treibhaus.Simpson_Fkt(xi2,eta,T)
                xi2 += 2*h
            u += 1
        Summe += treibhaus.Simpson_Fkt(1,eta,T) +treibhaus.Simpson_Fkt(l-h,eta,T)   #Randterme
        return (h/3)*Summe
         
    #input: Temperatur T und gewünschte Genauigkeit
    #Output: 1/analytische Lösung des Integrals über die spektrale Energiedichte rho
    #analytisch geslöstes Integral wird bestimmt mit einer Genauigkeit, die gegeben wird(aber 1/Z)
    def Z_Integral(Genauigkeit,T):
        Z = 1       #division by zero umgehen in while Schleife Kriterium  
        Wert = 0    #wird aktuallisiert für die Durchläufe
        m = 1       #Startwert der Näherung
        while (Z-Wert)/Z > Genauigkeit:
            if m == 1:   #devision by zero umgehung
                Z = 0
            Wert = Z
            Z += (48*np.pi*k**4*T**4)/(c**3*planck**3*m**4) #Summand der analytischen Lösung
            m += 1
        return 1/Z
        
    
    #input: l: Breite des Intervalls, n: Anzahl Stützstellen (Gerade), T: Temperatur,Genauigkeit, eta: n_schlange 
    #output: epsilon von gewünschter Temperatur und n_schlange #
    #Funktion Epsilon, die das produktder analytischen Lösung von 1/Z und die numerischen nach Simpson berechnet    
    def epsilon(l,n,T,eta,Genauigkeit):
        Wert = treibhaus.Z_Integral(Genauigkeit,T)*treibhaus.Integral_S(l,n,T,eta)
        return Wert
    
    #input: l: Breite des Intervalls, n: Anzahl Stützstellen (Gerade), Tstart: Starttemperatur,Genauigkeit, eta: n_schlange 
    #output: mittlere Erdtemperatur als Funktion von n schlange (tilde) in [K]
    #Fixpunktgleichung zum Lösen der Gleichung mit einen Unbekannten    
    def Fixpunkt(l,n,Tstart,Genauigkeit,eta):
        Summe1 = Tstart         #Starttemperatur
        Summe2 = Tstart
        while True:
            Summe1 = T_S*(epsilon_S*(2-treibhaus.epsilon(l,n,T_S,eta,Genauigkeit))/(2-treibhaus.epsilon(l,n,Summe1,eta,Genauigkeit)))**(1/4)
            if abs(Summe1-Summe2) < Genauigkeit:    #genauigkeit Erreicht
                break
            Summe2 = Summe1
        return Summe1
        

#nicht schön aber erfüllt hier den Zweck:
print("\n\nDie mittlere Erdtemperatur in °C mit n=1 ist somit:",treibhaus.Fixpunkt(1e16,100001,300,1e-9,1)-273.25)


