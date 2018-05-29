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
        Rechnung = ((8*np.pi*planck*nu**3)/(c**3)) * 1/((np.exp((planck*nu)/(k*T)))-1)
        return Rechnung   
    
    #input: Frequent nu
    #output: Wirkungsquerschnitt sigma
    #Formel für Wirkungsquerschnitt
    def sigma(nu):
        Wert = (S_2/np.pi)*Gamma/((nu-nu_2)**2+Gamma**2) + (S_3/np.pi)*Gamma/((nu-nu_3)**2+Gamma**2)
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
        return Wert

    #input; l: Breite des Intervalls, n: Anzahl Stützstellen (Gerade), T: Temperatur, eta: n schlange
    #output: Integral
    #numerische Berechnung des Integrals nach Simpson
    def Integral_S(l,n,T,eta):
        h = l / n           #Schrittweite
        Summe = 0   
        i = h               #laufindex
        m = 0
        while m < n:    #bis zur vorletzten Schrittweite
            if m % 2 == 0:  #Gerader Koeffizient    (für Gewichtung)
                Summe += 4*treibhaus.Simpson_Fkt(i,eta,T)
            else:           #ungerade Koeffizient
                Summe += 2*treibhaus.Simpson_Fkt(i,eta,T)
            i += h          #um Schrittweise höhere Index
            m += 1
        Summe += treibhaus.Simpson_Fkt(l,eta,T) #aufsummieren der Stützwerte am Beginn und am Ende der Fkt. (Beginn hier 0)
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
    #output: epsilon von gewünschter Temperatur und n_schlange 
    #Funktion Epsilon, die das produktder analytischen Lösung von 1/Z und die numerischen nach Simpson berechnet    
    def epsilon(l,n,T,eta,Genauigkeit):
        Wert = treibhaus.Z_Integral(Genauigkeit,T)*treibhaus.Integral_S(l,n,T,eta)
        return Wert
    
    #input: l: Breite des Intervalls, n: Anzahl Stützstellen (Gerade), Tstart: Starttemperatur,Genauigkeit, eta: n_schlange 
    #output: mittlere Erdtemperatur als Funktion von n schlange (tilde) in [K]
    #Fixpunktgleichung zum Lösen der Gleichung mit einen Unbekannten    
    def Fixpunkt(l,n,Tstart,Genauigkeit,eta):
        TE = Tstart #damit while schleife kein division by zero error
        Wert = 0
        T = Tstart  #Starttemperatur
        while (TE-Wert)/TE > Genauigkeit:
            Wert = TE
            TE = T_S*(epsilon_S*(2-treibhaus.epsilon(l,n,T_S,eta,Genauigkeit))/(2-treibhaus.epsilon(l,n,T,eta,Genauigkeit)))**(1/4)
            T = TE
        return TE
        
        

print(treibhaus.Fixpunkt(100,500,300,1e-9,1))

#eta:1: Zimmertemperatur
#eta:2: groeßenordnung e8 ca

