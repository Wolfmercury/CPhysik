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
    #input: Zahlenwert x, Zahlenwert y
    #output: Multiplikation davon
    #rechnet triviale rechnung
    def x_i(x,y):
        return x*y
    
    #input: Temperatur T, Wellenlaenge nu
    #output: rho
    #berechnung der spektralen Energiedichte fuer unpolarisierte Strahlung
    def rho(nu,T):
        Rechnung = ((8*np.pi*planck*nu**3)/c**3) * 1/(np.exp((planck*nu)/(k*T))-1)
        return Rechnung   
    
    #input: Temperatur T, Wellenlaenge nu
    #output: rho
    #berechnung der Ableitung der spektralen Energiedichte fuer unpolarisierte Strahlung
    def rho_diff(nu,T):
        Rechnung = (8*np.pi*planck*nu**3)/c**3 * planck*nu*np.exp(planck*nu/(k*T))/(k*T**2*(np.exp(planck*nu/(k*T))-1)**2)
        return Rechnung
    
    #input: Frequent nu
    #output: Wirkungsquerschnitt sigma
    #Formel für Wirkungsquerschnitt
    def sigma(nu):
        Wert = (S_2/np.pi*(Gamma/(((nu-nu_2)**2)+Gamma**2))) + (S_3/np.pi*(Gamma/(((nu-nu_3)**2)+Gamma**2)))
        return Wert
    
    #input: n_schlange (hier eta) und nu
    #output: Wahscheinlichkeit fuer das Entweichen eines Photons der Frequenz nu (f)
    #funktion im Integral    
    def f(eta,nu):
        return np.exp(-(eta*N_0*treibhaus.sigma(nu)))
        
    #input:Genauigkeit n, Temperaut T, länge l des Intervalls
    #output: Wert für bestimmte Genauigkeit
    #analyitisches Integral (epsilon(T))
    def Integral(l,n,T,eta):
        h = l/n     #stützstellen breite
        i = 1       #laufindex
        summe = 0   #gewünschte numerische Summe des Integrals
        while treibhaus.x_i(h,i) < n:
            summe += h*treibhaus.rho(treibhaus.x_i(h,i),T)*(1-treibhaus.f(eta,treibhaus.x_i(h,i)))
            i += 1
            #print(summe)
        return summe

    #Für das Integral nach Simpsons benötigte hilfsfunktion
    def Simpson_Fkt(nu,eta,T):
        Wert = treibhaus.rho(nu,T)*(1-treibhaus.f(eta,nu))
        return Wert

    #input; l: Breite des Intervalls, n: Anzahl Stützstellen (Gerade), T: Temperatur, eta: n schlange
    #output: Integral
    #numerische Berechnung des Integrals nach Simpson
    def Integral_S(l,n,T,eta):
        h = l / n       #Schrittweite
        #Summe = treibhaus.Simpson_Fkt(0,eta,T)
        Summe = 0
        i = h       #laufindex
        while i < n:
            if i % 2 == 0:  #Gerader Koeffizient    (für Gewichtung)
                Summe += 4*treibhaus.Simpson_Fkt(i,eta,T)
            else:           #ungerade Koeffizient
                Summe += 2*treibhaus.Simpson_Fkt(i,eta,T)
            i += h  #um Schrittweise höhere Index
        Summe += treibhaus.Simpson_Fkt(n,eta,T)
        return Summe
        


    #input:Genauigkeit n, Temperaut T, länge l des Intervalls
    #output: Wert für bestimmte Genauigkeit
    #Ableitung des analyitisches Integral, für das Newton Kriterium (epsilon(T)')
    def Integral_diff(l,n,T,eta):
        h = l/n     #stützstellen breite
        i = 1       #laufindex
        summe = 0   #gewünschte numerische Summe des Integrals
        while treibhaus.x_i(h,i) < n:
            summe += h*treibhaus.rho_diff(treibhaus.x_i(h,i),T)*(1-treibhaus.f(eta,treibhaus.x_i(h,i)))
            i += 1
        return summe
        
    #input: Temperatur T und gewünschte Genauigkeit
    #Output: 1/analytische Lösung des Integrals über die spektrale Energiedichte rho
    #analytisch geslöstes Integral wird bestimmt mit einer Genauigkeit, die gegeben wird(aber 1/Z)
    def Z_Integral(Genauigkeit,T):
        Z = 1       #division by zero umgehen in while Schleife Kriterium  
        Wert = 0    #wird aktuallisiert für die Durchläufe
        m = 1       #Startwert der Näherung
        while (Z-Wert)/Z > Genauigkeit:
            if m == 1:
                Z = 0
            Wert = Z
            Z += (48*np.pi*k**4*T**4)/(c**3*planck**3*m**4)
            m += 1
        return 1/Z
        
    #Analztisches Integral von Z ableiten
    def Z_Integral_diff(Genauigkeit,T):
        Z = 1       #division by zero umgehen in while Schleife Kriterium  
        Wert = 0    #wird aktuallisiert für die Durchläufe
        m = 1       #Startwert der Näherung
        while (Z-Wert)/Z > Genauigkeit:
            if m == 1:
                Z = 0
            Wert = Z
            Z += (192*np.pi*k**4*T**3)/(c**3*planck**3*m**4)
            m += 1
        return 1/Z
        
    #input: Startwerte: Genauigkeit i.A., Startschrittweite l, n ka, eta
    #output:
    #Newton Kriterium zum lösen der Gleichung mit einer Unbekannten
    def Newton(l,n,Tstart,Genauigkeit,eta):
        dT = 1      #Newtonkriterium
        T = Tstart
        while dT > Genauigkeit:
            print(treibhaus.Integral_S(l,n,T,eta)*treibhaus.Z_Integral(Genauigkeit,T))
            dT = treibhaus.Integral_S(l,n,T,eta)*treibhaus.Z_Integral(Genauigkeit,T) / (treibhaus.Integral_diff(l,n,T,eta)*treibhaus.Z_Integral(Genauigkeit,T)+treibhaus.Z_Integral_diff(Genauigkeit,T)*treibhaus.Integral_S(l,n,T,eta))
            T -= dT
            #print(T)
        return T

    
    def epsilon(l,n,T,eta,Genauigkeit):
        Wert = treibhaus.Z_Integral(Genauigkeit,T)*treibhaus.Integral_S(l,n,T,eta)
        return Wert
    
    
    
    def Fixpunktformel(l,n,Tstart,Genauigkeit,eta):
        TE = 1      #=1 damit while schleife kein division by zero error
        Wert = 0
        T = Tstart
        
        while (TE-Wert)/TE > Genauigkeit:
            if T == Tstart: #devision by zero umgehen
                TE = 0
            Wert = TE
            TE = T_S*((epsilon_S*(2-treibhaus.epsilon(l,n,T_S,eta,Genauigkeit)))/(2-treibhaus.epsilon(l,n,T,eta,Genauigkeit)))**(1/4)
            print(TE-Wert)
            T = TE
            print(T)
        return TE
        


#print(treibhaus.Integral_S(20,100,1,1))
#print(treibhaus.Z_Integral(1e-8,1))
#print(treibhaus.Integral_S(20,100,1,1)*treibhaus.Z_Integral(1e-8,1))

print(treibhaus.Fixpunktformel(1,10,1,1e-6,1))

#print(treibhaus.Newton(1,10,1,1e-6,1))

#BEI SUMME von INTEGRAL NUMERISCH DIE WHILE SCHLEIFE SO äNDERN WIE SONST AUCH
#eta:1: Zimmertemperatur
#eta:2: groeßenordnung e8 ca

