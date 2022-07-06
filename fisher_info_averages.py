import numpy as np
import matplotlib.pyplot as plt
import statistics as st

pi=3.14159265359 

##################################
##################################

def trayectoria(filename,fs,ncycles,num_files,dttemp,ntrans,margen):
    #tiempo_total = dttemp*(2*ncycles)
    tiempo_totalR = dttemp*(2*ncycles+1)
            
    #npoints_per_file = int(fs*tiempo_total) # Numero total de puntos por archivo
    npoints_per_fileR = int(fs*tiempo_totalR)
    npoints_per_semicycle = int(fs*dttemp) # Numero de puntos por semiciclo
    npoints_per_cycle = 2*npoints_per_semicycle # Numero de puntos por ciclo
    npoints_margin = int(fs*margen)
       
    ntrans = int(ntrans)
    
    ########
        
    pts_med_eq = np.zeros(2*ncycles)
    pts_med_eq[0] = npoints_margin + npoints_per_semicycle #Eq
    pts_med_eq[1] = npoints_per_cycle - 21 #Eq
    
    pts_med_transO = np.zeros(2*ncycles)
    pts_med_transO[0] = int(dttemp*fs+3)
    pts_med_transO[1] = pts_med_transO[0] + ntrans
    
    pts_med_transR = np.zeros(2*ncycles)
    pts_med_transR[0] = int(2*dttemp*fs+5)
    pts_med_transR[1] = pts_med_transR[0] + ntrans
    
    npts_eq = int(npoints_per_semicycle-npoints_margin-21)
    eq_series = np.zeros(ncycles*num_files*npts_eq)
    
    trans_array = np.zeros((ncycles*num_files,ntrans))
    trans_arrayR = np.zeros((ncycles*num_files,ntrans))
    
    # Calculadora de puntos a medir
    m = 2
    while m < 2*ncycles:
        pts_med_eq[m] = pts_med_eq[m-2] + npoints_per_cycle
        pts_med_transO[m] = pts_med_transO[m-2] + npoints_per_cycle
        pts_med_transR[m] = pts_med_transR[m-2] + npoints_per_cycle
        m += 1   
    
    mm = 0
    cFile = 1
    while cFile <= num_files:
        ##### Lectura del fichero numero cFile
        fid = open(filename + str(cFile) + '.out', 'r') 
        
        m = 0
        x = np.zeros(npoints_per_fileR)
        for i in range(8):
            fid, next(fid)
        while m <= npoints_per_fileR-1:
            line = fid.readline()
            line_split = line.split()
            ls1 = line_split[1]
            x[m] = float(ls1)
            
            m = m + 1
        
        #####
        
        # Serie de equilibrio
        m = 0
        while m < int(2*ncycles):
            tin = int(pts_med_eq[m])
            tin1 = int(pts_med_transO[m])
            tin3 = int(pts_med_transR[m])
            tfi = int(pts_med_eq[m+1])
            tfi1 = int(pts_med_transO[m+1])
            tfi3 = int(pts_med_transR[m+1])
       
            xsec_eq = x[tin:tfi]
            xsec_tr = x[tin1:tfi1]
            xsec_trR = x[tin3:tfi3]
            
            i0=int(mm*npts_eq)
            i1=int((mm+1)*npts_eq)
            eq_series[i0:i1] = xsec_eq
            trans_array[mm,:] = xsec_tr
            trans_arrayR[mm,:] = xsec_trR
            
            m=m+2
            mm=mm+1
        
        
        cFile += 1
        
    return trans_array,trans_arrayR,eq_series
    
def histogramas(bins,serie):
    hist,bins_hist = np.histogram(serie,bins=bins,density=True)
    
    return hist,bins_hist

def integrales(mDistr,nbins,ntrans,Dx,fs):
    info = np.zeros(ntrans-2)
    info2 = np.zeros(ntrans-2)

    n=1
    while n<ntrans-1:
        sum1 = 0
        sum2= 0
        
        xtv = mDistr[n,:]
        xt0v = mDistr[n-1,:]
        xt1v = mDistr[n+1,:]
        
        m=0
        while m<nbins-1:
            xt = xtv[m]
            xt0 = xt0v[m]
            xt1 = xt1v[m]
            
            if xt*xt0*xt1 == 0:
                m += 1
                continue
            
            logxt0 = np.log(xt0)
            logxt1 = np.log(xt1)
            
            sum1 = sum1 + (0.25*Dx*fs*fs)*xt*np.power(logxt1-logxt0,2)
            sum2 = sum2 + (0.5*Dx*fs)*xt*(logxt1-logxt0)
            m += 1
        
        info[n-1] = sum1 - sum2*sum2
        n += 1
        
    n=1
    while n<ntrans-1:
        integr2 = 0
        
        xtv = mDistr[n,:]
        xt0v = mDistr[n-1,:]
        xt1v = mDistr[n+1,:]
        
        m=0
        while m<nbins-1:
            xt = xtv[m]
            xt0 = xt0v[m]
            xt1 = xt1v[m]
            
            if xt*xt0*xt1==0:
                m += 1
                continue
            
            logxt0 = np.log(xt0)
            logxt1 = np.log(xt1)
            
            integr2 = integr2 + Dx*np.power(0.5*fs*(logxt1-logxt0),2)*xt
            
            m += 1
        
        info2[n-1] = integr2
        
        n += 1    
        
    return info,info2



def protocolo(dir_name,opt,fs,ncycles,num_files,dttemp,ntrans,nbins,kappa,Nc,Neq,margen):
    dir_name = dir_name + opt + '/'
    
    trans_array,trans_arrayR,eq_series = \
        trayectoria(dir_name,fs,ncycles,num_files,dttemp,ntrans,margen)
        
    bins= np.linspace(eq_series.min()-0.5e-8,eq_series.max()+0.5e-8,num=nbins)
    Dx = bins[1]-bins[0]

    # Histogramas transitorios
    i=0

    
    mDistr = np.zeros((ntrans,nbins-1))
    #eMatriz = np.zeros((ntrans,nbins))
    
    mDistrR = np.zeros((ntrans,nbins-1))
    #eMatrizR = np.zeros((ntrans,nbins))
    
    while i < ntrans:
        trans_series = trans_array[:,i]
        hist_trans,bins_trans = histogramas(bins,trans_series)
        #eHist_trans = np.sqrt(hist_trans / (Nc*Dx))
        
        mDistr[i,:] = hist_trans
        #eMatriz[i,:] = eHist_trans
        
        trans_series = trans_arrayR[:,i]
        hist_trans,bins_trans = histogramas(bins,trans_series)
        #eHist_trans = np.sqrt(hist_trans / (Nc*Dx))
        
        mDistrR[i,:] = hist_trans
        #eMatrizR[i,:] = eHist_trans
        
        i += 1

    infoF,infoF2 = integrales(mDistr,nbins,ntrans,Dx,fs)
    infoR,infoR2 = integrales(mDistrR,nbins,ntrans,Dx,fs)
            
    return infoF,infoF2,infoR,infoR2

def dibujos(serieHeat,serieCool,eHeat,eCool,ts,magnitude,opt,fname,N):
    plt.figure()
    plt.semilogx(ts[0:N],serieHeat[0:N],'.r')
    plt.semilogx(ts[0:N],serieCool[0:N],'.b')
    plt.legend(['Cold-warm','Hot-warm'])
    plt.xlabel("Time (s)")
    plt.ylabel(magnitude + " (a.u.)")
    plt.title("Evolution of " + magnitude + " (" + opt + " process)")
    plt.savefig(dir_name + fname + '-' + opt + '.png')

    plt.figure()
    ax2 = plt.axes()
    ax2.set_xscale("log")
    ax2.errorbar(ts[0:N],serieHeat[0:N],yerr=eHeat[0:N], fmt = '.r')
    ax2.errorbar(ts[0:N],serieCool[0:N],yerr=eCool[0:N], fmt = '.b')
    plt.legend(['Cold-warm','Hot-warm'])
    ax2.set_xlabel("Time (s)")
    plt.ylabel(magnitude + " (a.u.)")
    plt.title("Evolution of " + magnitude + " (" + opt + " process)")
    plt.savefig(dir_name + fname + '-' + opt + 'error.png')
    

##############################
##############################

dir_name = 'C:/lab/t4/'

#### Datos generales
fs = 50000  # Hz                # Frecuencia de sampleo
dttemp = 0.005 # s
ncycles = 2399
num_files = 10
margen = 0.003
ntrans = 100
#kappa = 41.64e-6 #tanda7
#kappa = 14.04e-6
kappa = 79.40e-6


Nc = ncycles*num_files
Neq = int(Nc*fs*(dttemp-margen)-21)
nbins = 60


##############################
##############################
        
#### CASO DIRECTO
print('DIRECTO')

ts = np.arange(ntrans-2)/fs

infoF = np.zeros(ntrans-2)
infoF2 = np.zeros(ntrans-2)
infoR = np.zeros(ntrans-2)
infoR2 = np.zeros(ntrans-2)

infoFI = np.zeros(ntrans-2)
infoFI2 = np.zeros(ntrans-2)
infoRI = np.zeros(ntrans-2)
infoRI2 = np.zeros(ntrans-2)

contador = 0
while contador<num_files:
    infoFv,infoF2v,infoRv,infoR2v = \
        protocolo(dir_name,'directo',fs,ncycles,1,dttemp,ntrans,nbins,kappa,Nc,Neq,margen)
        
    infoFIv,infoFI2v,infoRIv,infoRI2v = \
        protocolo(dir_name,'inverso',fs,ncycles,1,dttemp,ntrans,nbins,kappa,Nc,Neq,margen)  
    
    infoF = infoF + infoFv
    infoF2 = infoF2 + infoF2v
    infoR = infoR + infoRv
    infoR2 = infoR2 + infoR2v

    infoFI = infoFI + infoFIv
    infoFI2 = infoFI2 + infoFI2v
    infoRI = infoRI + infoRIv
    infoRI2 = infoRI2 + infoRI2v
    
    contador += 1
    print('Archivo ' + str(contador))
    
 

infoF = infoF/num_files
infoF2 = infoF2/num_files
infoR = infoR/num_files
infoR2 = infoR2/num_files

infoFI = infoFI/num_files
infoFI2 = infoFI2/num_files
infoRI = infoRI/num_files
infoRI2 = infoRI2/num_files
 
fich = open('infoAv.out','w')

for i in range(ntrans-2):
    fich.write(str(ts[i]) + ' ' + str(infoF[i]) + ' ' + str(infoF2[i]) + '\n')
    
fich.close()   
    
# LF = np.zeros(ntrans-2)
# LR = np.zeros(ntrans-2)
# LFI = np.zeros(ntrans-2)
# LRI = np.zeros(ntrans-2)

# n=0
# while n<ntrans-2:
#     LF[n] = np.sum(np.sqrt(infoF[0:n])/fs)
#     LR[n] = np.sum(np.sqrt(infoR[0:n])/fs)
#     LFI[n] = np.sum(np.sqrt(infoFI[0:n])/fs)
#     LRI[n] = np.sum(np.sqrt(infoRI[0:n])/fs)
    
#     n +=1

# plt.figure()
# plt.plot(ts,LF,'.',ts,LR,'.')

# plt.figure()
# plt.plot(ts,LFI,'.',ts,LRI,'.')

# plt.figure()
# plt.plot(ts,infoF2,'.',ts,infoR2,'.')
# a = 0.182-1
# w = 2*pi*863.3
# ysF = 2*w*w*a*a*np.exp(-4*w*ts) / (1+a*np.exp(-2*w*ts))
# ysR = 2*w*w*a*a*np.exp(-4*w*ts) / (1+a-a*np.exp(-2*w*ts))
# plt.plot(ts,ysF,'--')
# plt.plot(ts,ysR,'--')

# plt.figure()
# plt.plot(ts,infoFI2,'.',ts,infoRI2,'.')