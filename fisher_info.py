import numpy as np
import matplotlib.pyplot as plt
import statistics as st

##################################
##################################

def trayectoria(filename,fs,ncycles,num_files,dttemp,ntrans,margen):
    tiempo_total = dttemp*(2*ncycles)
    tiempo_totalR = dttemp*(2*ncycles+1)
            
    npoints_per_file = int(fs*tiempo_total) # Numero total de puntos por archivo
    npoints_per_fileR = int(fs*tiempo_totalR)
    npoints_per_semicycle = int(fs*dttemp) # Numero de puntos por semiciclo
    npoints_per_cycle = 2*npoints_per_semicycle # Numero de puntos por ciclo
    npoints_margin = int(fs*margen)
       
    ntrans = int(ntrans)
    
    ########
    
    #### Vector de puntos a medir
    pts_med_low = np.zeros(2*ncycles)
    pts_med_low[0] = npoints_margin 
    pts_med_low[1] = npoints_per_semicycle - 21 #Eq
    
    pts_med_eq = np.zeros(2*ncycles)
    pts_med_eq[0] = npoints_margin + npoints_per_semicycle #Eq
    pts_med_eq[1] = npoints_per_cycle - 21 #Eq
    
    pts_med_transO = np.zeros(2*ncycles)
    pts_med_transO[0] = int(dttemp*fs-5)
    pts_med_transO[1] = pts_med_transO[0] + ntrans
    
    pts_med_transR = np.zeros(2*ncycles)
    pts_med_transR[0] = int(2*dttemp*fs-5)
    pts_med_transR[1] = pts_med_transR[0] + ntrans
    
    npts_eq = int(npoints_per_semicycle-npoints_margin-21)
    eq_series = np.zeros(ncycles*num_files*npts_eq)
    low_series = np.zeros(ncycles*num_files*npts_eq)
    
    trans_array = np.zeros((ncycles*num_files,ntrans))
    trans_arrayR = np.zeros((ncycles*num_files,ntrans))
    
    #### Vector de temperaturas
    temp1 = np.zeros(ncycles*num_files)
    temp2 = np.zeros(ncycles*num_files)
    
    # Calculadora de puntos a medir
    m = 2
    while m < 2*ncycles:
        pts_med_eq[m] = pts_med_eq[m-2] + npoints_per_cycle
        pts_med_transO[m] = pts_med_transO[m-2] + npoints_per_cycle
        pts_med_transR[m] = pts_med_transR[m-2] + npoints_per_cycle
        pts_med_low[m] = pts_med_low[m-2] + npoints_per_cycle
        m += 1
    
    m = 2
    while m < 2*ncycles-1:
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
            tin2 = int(pts_med_low[m])
            tin3 = int(pts_med_transR[m])
            tfi = int(pts_med_eq[m+1])
            tfi1 = int(pts_med_transO[m+1])
            tfi2 = int(pts_med_low[m+1])
            tfi3 = int(pts_med_transR[m+1])
       
            xsec_eq = x[tin:tfi]
            xsec_tr = x[tin1:tfi1]
            xsec_trR = x[tin3:tfi3]
            xsec_low = x[tin2:tfi2]
            
            temp1[mm] = st.variance(xsec_low)
            temp2[mm] = st.variance(xsec_eq)
            
            i0=int(mm*npts_eq)
            i1=int((mm+1)*npts_eq)
            eq_series[i0:i1] = xsec_eq
            trans_array[mm,:] = xsec_tr
            trans_arrayR[mm,:] = xsec_trR
            low_series[i0:i1] = xsec_low
            
            m=m+2
            mm=mm+1
        
        print('Archivo ' + str(cFile))
        cFile += 1
        
    ## Borrar ahora
    ts = np.arange(ncycles*num_files*npts_eq)/50000
        
    fich = open(dir_name + 'trace.txt','w')
            
    for i in range(ncycles*num_files*npts_eq):
        fich.write(str(ts[i]) + '     \t\t' + str(eq_series[i]) + \
                       '     \t\t' + str(low_series[i]))
        fich.write('\n')

    fich.close()
        
    return eq_series,trans_array,trans_arrayR,low_series,temp1,temp2
    
def histogramas(bins,serie):
    hist,bins_hist = np.histogram(serie,bins=bins,density=True)
    
    return hist,bins_hist

def integrales(mDistr,nbins,ntrans,Dx,fs):
    info = np.zeros(ntrans-2)
    
    
    
    m=0
    n=1
    while n<ntrans-2:
        integr = 0
        
        xtv = mDistr[n,:]
        xt0v = mDistr[n-1,:]
        xt1v = mDistr[n+1,:]
        
        m=0
        while m<nbins-1:
            xt = xtv[m]
            xt0 = xt0v[m]
            xt1 = xt1v[m]
            
            if xt == 0:
                m += 1
                continue
            
            integr = integr + Dx*np.power(0.5*fs*(xt1-xt0),2)/xt
            
            m += 1
        
        info[n] = integr
        
        n += 1
        
    return info



def protocolo(dir_name,opt,fs,ncycles,num_files,dttemp,ntrans,margen,nbins,kappa,Nc,Neq):
    dir_name = dir_name + opt + '/'
    
    eq_series,trans_array,trans_arrayR,low_series,temp1D,temp2D = \
        trayectoria(dir_name,fs,ncycles,num_files,dttemp,ntrans,margen)
        
    bins= np.linspace(eq_series.min()-0.5e-8,eq_series.max()+0.5e-8,num=nbins)
    Dx = bins[1]-bins[0]

    # Histogramas transitorios
    i=0
    divDV = np.zeros(ntrans)
    eDivDV = np.zeros(ntrans)
    divDVr = np.zeros(ntrans)
    eDivDVr = np.zeros(ntrans)

    energDV = np.zeros(ntrans) # NEW!
    eEnergDV = np.zeros(ntrans)
    energDVr = np.zeros(ntrans)
    eEnergDVr = np.zeros(ntrans)
    
    SDV = np.zeros(ntrans)
    eSDV = np.zeros(ntrans)
    SDVr = np.zeros(ntrans)
    eSDVr = np.zeros(ntrans)
    
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

    infoF = integrales(mDistr,nbins,ntrans,Dx,fs)
    infoR = integrales(mDistrR,nbins,ntrans,Dx,fs)
            
    return infoF,infoR

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

dir_name = 'C:/Users/DELL/Desktop/LAB/23062022/set1/t1/'

#### Datos generales
fs = 50000  # Hz                # Frecuencia de sampleo
dttemp = 0.005 # s
ncycles = 2399
num_files = 10
margen = 0.003
ntrans = 55
#kappa = 41.64e-6 #tanda7
#kappa = 14.04e-6
kappa = 79.40e-6


Nc = ncycles*num_files
Neq = int(Nc*fs*(dttemp-margen)-21)
nbins = 75


##############################
##############################
        
#### CASO DIRECTO
print('DIRECTO')

ts = np.arange(ntrans-2)/fs

infoF,infoR = \
    protocolo(dir_name,'directo',fs,ncycles,num_files,dttemp,ntrans,margen,nbins,kappa,Nc,Neq) 

