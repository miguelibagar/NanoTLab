import numpy as np
import matplotlib.pyplot as plt
import statistics as st

##################################
##################################

pi = 3.14159265359
kB = 1.380649e-23

def importa(filename,fs,ncycles,num_files,A,fosc,kappa):
    dttemp = 1/fosc    
    tiempo_total = dttemp*ncycles
            
    npoints_per_file = int(fs*tiempo_total) # Numero total de puntos por archivo
    npoints_per_cycle = int(fs*dttemp) # Numero de puntos por ciclo
    
    ########
    
    #### Vector de puntos a medir
    pts_med = np.zeros(ncycles+1)
    pts_med[0] = 0 
    
    medias = np.zeros(npoints_per_cycle)
    momento2 = np.zeros(npoints_per_cycle)
    
    # Calculadora de puntos a medir
    m = 1
    while m < ncycles+1:
        pts_med[m] = pts_med[m-1] + npoints_per_cycle
        m += 1
    
    cFile = 1
    while cFile <= num_files:
        ##### Lectura del fichero numero cFile
        fid = open(filename + str(cFile) + '.txt', 'r')
        
        m = 0
        x = np.zeros(npoints_per_file)
        for i in range(101):
            fid, next(fid)
        while m < npoints_per_file:
            line = fid.readline()
            line_split = line.split()
            ls1 = line_split[0]
            x[m] = float(ls1)
            
            m = m + 1
        
        #####
        
        x = x/kappa
        
        ### Serie de equilibrio
        m = 0
        while m < int(ncycles):
            tin = int(pts_med[m])
            tfi = int(pts_med[m+1])
            
            xsec = x[tin:tfi]
            
            medias = medias + xsec
            momento2 = momento2 + np.power(xsec,2)
            
            m += 1
        
        print('Archivo ' + str(cFile))
        cFile += 1
    
    ts = np.arange(npoints_per_cycle)/fs
    x0 = A*np.sin(2*pi*fosc*ts)
    dtx0 = 2*pi*fosc*A*np.cos(2*pi*fosc*ts)
    
    medias = medias/(num_files*ncycles) 
    momento2 = momento2/(num_files*ncycles)
    
    return medias,momento2,x0,dtx0,ts

##############################
##############################

dir_name = 'C:/Users/Miguel Ibáñez/Desktop/ruben/'

#### Datos generales
fs = 50000  # Hz                # Frecuencia de sampleo
ncycles = 200
num_files = 5

fosc = 10 # Hz
A = 4e-7
kappa = 1.075e-5
T = 299.15

rp = 0.5e-6
eta = 0.0008701
#eta = 1e-3
gamma = 6*pi*eta*rp

medias,momento2,x0,dtx0,ts = importa(dir_name,fs,ncycles,num_files,A,fosc,kappa)



medias = medias - x0
varianzas = momento2 + np.power(x0,2) - 2*x0*medias



dtw = -kappa*dtx0*medias - kappa*x0*dtx0
dtq = kappa*(kB*T-kappa*varianzas)/gamma

pr_dtw = np.mean(dtw)
pr_dtq = np.mean(dtq)

print('Trabajo = ' + str(pr_dtw))
print('Calor = ' + str(pr_dtq))

w2 = np.power(2*pi*fosc,2)
A2 = np.power(A,2)
gamma2 = np.power(gamma,2)
kappa2 = np.power(kappa,2)

teorico = -0.5*((gamma*w2*A2) / (1 + w2*gamma2/kappa2))
print('Teorico = ' + str(teorico))



