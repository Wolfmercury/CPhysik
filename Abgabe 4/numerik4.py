# coding: utf8
import numpy as np
'''
// python3 numerik4.py
// Kilian Kranz, Achim Bruchertseifer
'''

'''
TODO:
    (o) lösen der Phi Gleichung
    (o) integrieren der Phi gleichung (falls ein eigener Schritt)
    (o) damit r(phi(t_i)) kriegen
    (o) r dimensionslos machen
    (o) T bestimmen
    (o) letzte Theta gleichung lösen
    (o) Diagramm plotten (besondere)
'''

epsilon = 0.1

def f(y0,y1):
    Rechnung = -2*(epsilon*np.sin(y0))/(epsilon*np.cos(y0)+1)*y1**2
    return Rechnung

#input:
#output:
#Runge Kutter Verfahren 4           (vektorweise mit y_0 oben, y_1 unten)
def RK4(h):
    Wert_y0 = 0
    Wert_y1 = 0
    #Startwerte/Anfagsbed.
    y0 = 2
    y1 = 5
    while (y0-Wert_y0)/y0 > 1e-6 or (y1-Wert_y1)/y1 > 1e-6:
        Wert_y0 = y0
        Wert_y1 = y1
        y0_neu = y0 + RK4_y0(h,y0,y1) #neu, da y1 alten y0 Wert 'braucht'
        y1 += RK4_y1(h,y0,y1)         #RK Schritt fuer y1
        y0 = y0_neu                   #nun darf es gleich werden
        #print(y0)
    return (y0,y1)
        
def RK4_y0(h,y0,y1):
    k1 = h*f(y0,y1)
    k2 = h*(y1+k1/2)    #f0 = y1
    k3 = h*(y1+k2/2)
    k4 = h*(y1+k3)
    
    return (1/6*(k1+2*k2+2*k3+k4))

def RK4_y1(h,y0,y1):
    k1 = h*y1
    k2 = h*f(y0+k1/2,y1)
    k3 = h*f(y0+k2/2,y1)
    k4 = h*f(y0+k3,y1)
    
    return (1/6*(k1+2*k2+2*k3+k4))

print(RK4(0.0025))
