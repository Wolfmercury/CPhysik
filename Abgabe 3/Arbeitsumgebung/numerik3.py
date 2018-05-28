# python (ggf. "python3")
# coding: utf8
import numpy as np
'''
// python3 Numerik1.py
// Kilian Kranz, Achim Bruchertseifer


'''
#Konstanten:
planck = 6.626070040e-34    #planckkonstante
c =  299792458              #lichtgeschwindigkeit m/s
k = 1.38064852e-23          #Boltzmann const J K^-1
nu_2 = 2e13                 #von sigma das 1/s
nu_3 = 7.04e13              #von sigma das 1/s
Gamma = 3e10                #GAmme von sigma 1/s
S_2 = 2.45e-11              #austrahlung dings m`2/s
S_3 = 2.74e-10              #austrahlung dings m`2/s
epsilon_S = 5.5e-6          #emessivitaet
T_S = 5750                  #Temperatur der Sonne in K
N_0 = 8.1e25                #Anzahl an momentanen CO2 molekuelen in luftsaule mit 1 m^2 querschnitt in m^-2


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
        Rechnung = (8*np.pi*planck*nu**3)/c**3 * 1/(np.exp((planck*nu)/(k*T))-1)
        return Rechnung   
    
    #input: Temperatur T, Wellenlaenge nu
    #output: rho
    #berechnung der Ableitung der spektralen Energiedichte fuer unpolarisierte Strahlung
    def rho_diff(nu,T):
        Rechnung = (8*np.pi*planck*nu**3)/c**3 * planck*nu*np.exp(planck*nu/(k*T))/(k*T**2*(np.exp(planck*nu/(k*T))-1)**2)
        return Rechnung
    
    #input:nu
    #output:sigma
    #Formel für Wirkungsquerschnitt
    def sigma(nu):
        Wert = S_2/np.pi*(Gamma/(((nu-nu_2)**2)+Gamma**2))+S_3/np.pi*(Gamma/(((nu-nu_3)**2)+Gamma**2))
        return Wert
    

    #input: n_schlange (hier eta) und nu
    #output: f
    #funktion im Integral    
    def f(eta,nu):
        return np.exp(-eta*N_0*treibhaus.sigma(nu))


    #input:Genauigkeit n, Temperaut T, länge l des Intervalls
    #output: Wert für bestimmte Genauigkeit
    #analyitisches Integral (epsilon(T))
    def Integral(l,n,T,eta):
        h = l/n     #stützstellen breite
        i = 1       #laufindex
        summe = 0   #gewünschte analytische Summe des Integrals
        while treibhaus.x_i(h,i) < n:
            summe += h*treibhaus.rho(treibhaus.x_i(h,i),T)*(1-treibhaus.f(eta,treibhaus.x_i(h,i)))
            i += 1
        return summe

    #input:Genauigkeit n, Temperaut T, länge l des Intervalls
    #output: Wert für bestimmte Genauigkeit
    #Ableitung des analyitisches Integral, für das Newton Kriterium (epsilon(T)')
    def Integral_diff(l,n,T,eta):
        h = l/n     #stützstellen breite
        i = 1       #laufindex
        summe = 0   #gewünschte analytische Summe des Integrals
        while treibhaus.x_i(h,i) < n:
            summe += h*treibhaus.rho_diff(treibhaus.x_i(h,i),T)*(1-treibhaus.f(eta,treibhaus.x_i(h,i)))
            i += 1
        return summe
        
    #input: Temperatur T und gewünschte Genauigkeit
    #Output: 1/analytische Lösung des Integrals über die spektrale Energiedichte rho
    #analytisch geslöstes Integral wird bestimmt mit einer Genauigkeit, die gegeben wird(aber 1/Z)
    def Z_Integral(Genauigkeit,T):
        Z = 1   #division by zero umgehen in while Schleife Kriterium  
        Wert = 0  #wird aktuallisiert für die Durchläufe
        m = 1   #Startwert der Näherung
        while (Z-Wert)/Z > Genauigkeit:
            if m == 1:
                Z = 0
            Wert = Z
            Z += (48*np.pi*k**4*T**4)/(c**3*planck**3*m**4)
            m += 1
        return 1/Z
        
    #1/Z abgeleitet
    def Z_Integral_diff(Genauigkeit,T):
        Z = 1   #division by zero umgehen in while Schleife Kriterium  
        Wert = 0  #wird aktuallisiert für die Durchläufe
        m = 1   #Startwert der Näherung
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
        while dT>Genauigkeit:
            dT = treibhaus.Integral(l,n,T,eta)*treibhaus.Z_Integral(Genauigkeit,T) / (treibhaus.Integral_diff(l,n,T,eta)*treibhaus.Z_Integral(Genauigkeit,T)+treibhaus.Z_Integral_diff(Genauigkeit,T)*treibhaus.Integral(l,n,T,eta))
            T -= dT
        return T



