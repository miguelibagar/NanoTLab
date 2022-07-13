import numpy as np
import matplotlib.pyplot as plt

#############################################################################
#############################################################################
#############################################################################
#################### FUNCION GENERAL DE IMPORTACION #########################
#############################################################################
#############################################################################
#############################################################################

def importa(filename,fs,ncycles,num_files,dttemp,ntrans,margen,contador):
    #tiempo_total = dttemp*(2*ncycles)
    tiempo_totalR = dttemp*(2*ncycles+1)
            
    #npoints_per_file = int(fs*tiempo_total) # Numero total de puntos por archivo
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
    
    noise = np.zeros((ncycles*num_files,ntrans))
    noiseR = np.zeros((ncycles*num_files,ntrans))
    
    # Calculadora de puntos a medir
    m = 2
    while m < 2*ncycles:
        pts_med_eq[m] = pts_med_eq[m-2] + npoints_per_cycle
        pts_med_transO[m] = pts_med_transO[m-2] + npoints_per_cycle
        pts_med_transR[m] = pts_med_transR[m-2] + npoints_per_cycle
        pts_med_low[m] = pts_med_low[m-2] + npoints_per_cycle
        m += 1
    
    mm = 0
    cFile = 1
    while cFile <= num_files:
        ##### Lectura del fichero numero cFile
        fid = open(filename + str(cFile+contador) + '.out', 'r') 
        
        m = 0
        x = np.zeros(npoints_per_fileR)
        y = np.zeros(npoints_per_fileR) # Noisy signal
        
        for i in range(8):
            fid, next(fid)
        while m <= npoints_per_fileR-1:
            line = fid.readline()
            line_split = line.split()
            ls1 = line_split[1]
            ls13 = line_split[12]
            x[m] = float(ls1)
            y[m] = float(ls13)
            
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
            
            ysec_tr = y[tin1:tfi1]
            ysec_trR = y[tin3:tfi3]
            
            i0=int(mm*npts_eq)
            i1=int((mm+1)*npts_eq)
            eq_series[i0:i1] = xsec_eq
            trans_array[mm,:] = xsec_tr
            trans_arrayR[mm,:] = xsec_trR
            low_series[i0:i1] = xsec_low
            
            noise[mm,:] = ysec_tr
            noiseR[mm,:] = ysec_trR
            
            m=m+2
            mm=mm+1
        
        print('Archivo ' + str(contador+1))
        cFile += 1
        
    return eq_series,trans_array,trans_arrayR,low_series,noise,noiseR

#############################################################################
#############################################################################
############################ HISTOGRAMAS ####################################
#############################################################################
#############################################################################
    
def histogramas(bins,serie):
    hist,bins_hist = np.histogram(serie,bins=bins,density=True)
    
    return hist,bins_hist

#############################################################################
#############################################################################
############################ KULLBACK-LEIBLER ###############################
#############################################################################
#############################################################################

def kullback_leibler(distr1,eDistr1,distrRef,eDistrRef,nbins,Dx):
    kld = 0
    eKld = 0
    energ = 0
    eEnerg = 0
    
    m=0
    while m<nbins-1:
        xt = distr1[m]
        ext = eDistr1[m]
        xeq = distrRef[m]
        exeq = eDistrRef[m]
        
        if xt == 0 or xeq == 0:
            m += 1
                        
            continue
        
        kld = kld + Dx*xt*np.log(xt/xeq)
        eKld = eKld + np.power(ext*(1+np.log(xt/xeq)),2) + \
                 np.power(exeq*xt/xeq,2)
        energ = energ - Dx*xt*np.log(xeq) # NEW!
        eEnerg = eEnerg + np.power(ext*(np.log(xeq)),2) + \
                   np.power(exeq*xt/xeq,2)
        m += 1
    
    S = energ - kld
    
    eKld = np.sqrt(eKld)*Dx
    eEnerg = np.sqrt(eEnerg)*Dx
    eS = eEnerg + eKld
    
    return kld,eKld,energ,eEnerg,S,eS

#############################################################################
#############################################################################
#############################################################################
#############################################################################
############################ PROTOCOLO GODEC ################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################

def protocolo_godec(dir_name,opt,fs,ncycles,num_files,dttemp,ntrans,margen,nbins,kappa,Nc,Neq,contador):
    dir_name = dir_name + opt + '/'
    
    eq_series,trans_array,trans_arrayR,low_series,_,_ = \
        importa(dir_name,fs,ncycles,num_files,dttemp,ntrans,margen,contador)
    
    ## Creacion del vector de bins
    max1 = trans_array.max()
    max2 = trans_arrayR.max()
    max3 = eq_series.max()
    max4 = low_series.max()
    
    min1 = trans_array.min()
    min2 = trans_arrayR.min()
    min3 = eq_series.min()
    min4 = low_series.min()
    
    xmax = np.max(np.array([max1,max2,max3,max4]))
    xmin = np.min(np.array([min1,min2,min3,min4]))
    
    bins = np.linspace(xmin-0.05e-8,xmax+0.05e-8,num=nbins)
    Dx = bins[1]-bins[0]    
    
    # bins= np.linspace(eq_series.min()-0.5e-8,eq_series.max()+0.5e-8,num=nbins)
    # Dx = bins[1]-bins[0]

    hist_eq,bins_eq = histogramas(bins,eq_series)
    eHist_eq = np.sqrt(hist_eq / (Neq*Dx))

    hist_low,bins_low = histogramas(bins,low_series)
    eHist_low = np.sqrt(hist_low / (Neq*Dx))

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

    while i < ntrans:
        trans_series = trans_array[:,i]
        hist_trans,bins_trans = histogramas(bins,trans_series)
        eHist_trans = np.sqrt(hist_trans / (Nc*Dx))
        
        kld,eKld,energ,eEnerg,S,eS = \
            kullback_leibler(hist_trans,eHist_trans,hist_eq,eHist_eq,nbins,Dx)
        
        divDV[i] = kld
        eDivDV[i] = eKld
        energDV[i] = energ
        eEnergDV[i] = eEnerg
        SDV[i] = S
        eSDV[i] = eS
        
        trans_series = trans_arrayR[:,i]
        hist_trans,bins_trans = histogramas(bins,trans_series)
        eHist_trans = np.sqrt(hist_trans / (Nc*Dx))
        
        kld,eKld,energ,eEnerg,S,eS = \
            kullback_leibler(hist_trans,eHist_trans,hist_eq,eHist_eq,nbins,Dx)
        
        divDVr[i] = kld
        eDivDVr[i] = eKld
        energDVr[i] = energ
        eEnergDVr[i] = eEnerg
        SDVr[i] = S
        eSDVr[i] = eS
        
        i += 1

    return divDV,eDivDV,divDVr,eDivDVr,energDV,eEnergDV,energDVr,\
        eEnergDVr,SDV,eSDV,SDVr,eSDVr

#############################################################################
#############################################################################
############################ MEDIAS GODEC ###################################
#############################################################################
#############################################################################

def medias_godec(dir_name,opt,fs,ncycles,num_files,dttemp,ntrans,margen,nbins,kappa,Nc,Neq):
    divDV = np.zeros(ntrans)
    eDivDV = np.zeros(ntrans)
    divDVr = np.zeros(ntrans)
    eDivDVr = np.zeros(ntrans)
    
    energDV = np.zeros(ntrans)
    eEnergDV = np.zeros(ntrans)
    energDVr = np.zeros(ntrans)
    eEnergDVr = np.zeros(ntrans)
    
    SDV = np.zeros(ntrans)
    eSDV = np.zeros(ntrans)
    SDVr = np.zeros(ntrans)
    eSDVr = np.zeros(ntrans)
    
    contador = 0
    while contador < num_files:
        divDVV,eDivDVV,divDVVr,eDivDVVr,energDVV,eEnergDVV,energDVVr,eEnergDVVr,SDVV,\
            eSDVV,SDVVr,eSDVVr =\
            protocolo_godec(dir_name,opt,fs,ncycles,1,dttemp,ntrans,margen,nbins,kappa,Nc,Neq,contador) 
        
        divDV = divDV + divDVV
        eDivDV = eDivDV + eDivDVV
        divDVr = divDVr + divDVVr
        eDivDVr = eDivDVr + eDivDVVr
        
        energDV = energDV + energDVV
        eEnergDV = eEnergDV + eEnergDVV
        energDVr = energDVr + energDVVr
        eEnergDVr = eEnergDVr + eEnergDVVr
        
        SDV = SDV + SDVV
        eSDV = eSDV + eSDVV
        SDVr = SDVr + SDVVr
        eSDVr = eSDVr + eSDVVr
        
        contador += 1
        
    divDV = divDV/num_files
    eDivDV = eDivDV/num_files
    divDVr = divDVr/num_files
    eDivDVr = eDivDVr/num_files
    
    energDV = energDV/num_files
    eEnergDV = eEnergDV/num_files
    energDVr = energDVr/num_files
    eEnergDVr = eEnergDVr/num_files
    
    SDV = SDV/num_files
    eSDV = eSDV/num_files
    SDVr = SDVr/num_files
    eSDVr = eSDVr/num_files

    return divDV,eDivDV,divDVr,eDivDVr,energDV,eEnergDV,energDVr,eEnergDVr,SDV,eSDV,SDVr,eSDVr

#############################################################################
#############################################################################
#############################################################################
#############################################################################
######################## ANALISIS DE CUMULANTES #############################
#############################################################################
#############################################################################
#############################################################################
#############################################################################

def cumulantes(trans_array,ntrans):
    media = np.zeros(ntrans)
    varianza = np.zeros(ntrans)
    cum3 = np.zeros(ntrans)
    cum4 = np.zeros(ntrans)
    
    i=0
    while i < ntrans:
        trans_series = trans_array[:,i]
        
        media[i] = np.mean(trans_series)
        varianza[i] = np.var(trans_series)
        cum3[i] = np.mean(np.power(trans_series-media[i],3))
        cum4[i] = np.mean(np.power(trans_series-media[i],4)) - 3*varianza[i]*varianza[i]  
        
        i+=1
    
    return media,varianza,cum3,cum4

#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################

def cum_posicion(dir_name,opt,fs,ncycles,num_files,dttemp,ntrans,margen):
    dir_name = dir_name + opt + '/'
    
    _,trans_array,trans_arrayR,_,noise,noiseR = \
        importa(dir_name,fs,ncycles,num_files,dttemp,ntrans,margen,0)
    
    media,varianza,cum3,cum4 = cumulantes(trans_array,ntrans)
    mediaR,varianzaR,cum3R,cum4R = cumulantes(trans_arrayR,ntrans)
    
    mediaN,varianzaN,_,_ = cumulantes(noise,ntrans)
    mediaNR,varianzaNR,_,_ = cumulantes(noiseR,ntrans)
    
    return media,mediaR,varianza,varianzaR,cum3,cum3R,cum4,cum4R,mediaN,mediaNR,varianzaN,varianzaNR

#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################

def kld_cota_relax(dir_name,opt,fs,ncycles,num_files,dttemp,ntrans,margen,nbins,kappa,Nc,Neq):
    dir_name = dir_name + opt + '/'
    
    eq_series,trans_array,trans_arrayR,low_series,_,_ = \
        importa(dir_name,fs,ncycles,num_files,dttemp,ntrans,margen,0)
    
    ## Creacion del vector de bins
    max1 = trans_array.max()
    max2 = trans_arrayR.max()
    max3 = eq_series.max()
    max4 = low_series.max()
    
    min1 = trans_array.min()
    min2 = trans_arrayR.min()
    min3 = eq_series.min()
    min4 = low_series.min()
    
    xmax = np.max(np.array([max1,max2,max3,max4]))
    xmin = np.min(np.array([min1,min2,min3,min4]))
    
    bins = np.linspace(xmin-0.05e-8,xmax+0.05e-8,num=nbins)
    Dx = bins[1]-bins[0]    
    
    # bins= np.linspace(eq_series.min()-0.5e-8,eq_series.max()+0.5e-8,num=nbins)
    # Dx = bins[1]-bins[0]

    hist_eq,bins_eq = histogramas(bins,eq_series)
    eHist_eq = np.sqrt(hist_eq / (Neq*Dx))

    hist_low,bins_low = histogramas(bins,low_series)
    eHist_low = np.sqrt(hist_low / (Neq*Dx))

    # Histogramas transitorios
    i=0
    divDV = np.zeros(ntrans)
    eDivDV = np.zeros(ntrans)
    divDVr = np.zeros(ntrans)
    eDivDVr = np.zeros(ntrans)

    while i < ntrans:
        trans_series = trans_array[:,i]
        hist_trans,bins_trans = histogramas(bins,trans_series)
        eHist_trans = np.sqrt(hist_trans / (Nc*Dx))
        
        kld,eKld,_,_,_,_ = \
            kullback_leibler(hist_low,eHist_low,hist_trans,eHist_trans,nbins,Dx)
        
        divDV[i] = kld
        eDivDV[i] = eKld
        
        trans_series = trans_arrayR[:,i]
        hist_trans,bins_trans = histogramas(bins,trans_series)
        eHist_trans = np.sqrt(hist_trans / (Nc*Dx))
        
        kld,eKld,_,_,_,_ = \
            kullback_leibler(hist_eq,eHist_eq,hist_trans,eHist_trans,nbins,Dx)
        
        divDVr[i] = kld
        eDivDVr[i] = eKld
        
        i += 1

    return divDV,eDivDV,divDVr,eDivDVr

#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################

def dibujos_godec(serieHeat,serieCool,eHeat,eCool,ts,magnitude,opt,fname,N,ymax):
    plt.figure()
    plt.semilogx(ts[0:N],serieHeat[0:N],'.r')
    plt.semilogx(ts[0:N],serieCool[0:N],'.b')
    plt.legend(['Cold-warm','Hot-warm'])
    plt.xlabel("Time (s)")
    plt.ylabel(magnitude + " (a.u.)")
    plt.ylim(0,ymax)
    plt.title("Evolution of " + magnitude + " (" + opt + " process)")
    plt.savefig(dir_name + fname + '-' + opt + '-means.png')

    plt.figure()
    ax2 = plt.axes()
    ax2.set_xscale("log")
    ax2.errorbar(ts[0:N],serieHeat[0:N],yerr=eHeat[0:N], fmt = '.r')
    ax2.errorbar(ts[0:N],serieCool[0:N],yerr=eCool[0:N], fmt = '.b')
    plt.legend(['Cold-warm','Hot-warm'])
    ax2.set_xlabel("Time (s)")
    plt.ylabel(magnitude + " (a.u.)")
    plt.ylim(0,ymax)
    plt.title("Evolution of " + magnitude + " (" + opt + " process)")
    plt.savefig(dir_name + fname + '-' + opt + 'error-means.png')
    

#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################

dir_name = 'C:/lab/t4/'

#### Datos generales
fs = 50000  # Hz                # Frecuencia de sampleo
dttemp = 0.005 # s
ncycles = 2399
num_files = 10
margen = 0.003
ntrans = 65
kappa = 41.64e-6 #tanda7
#kappa = 14.04e-6
#kappa = 79.40e-6
ymax = 0.45


Nc = ncycles*num_files
Neq = int(Nc*fs*(dttemp-margen)-21)
nbins = 75


#############################################################################
#############################################################################
#############################################################################
#############################################################################
        
#### HEATING - GODEC
print('DIRECTO')

divDV,eDivDV,divDVr,eDivDVr,energDV,eEnergDV,energDVr,eEnergDVr,SDV,eSDV,SDVr,eSDVr =\
    medias_godec(dir_name,'directo',fs,ncycles,num_files,dttemp,ntrans,margen,nbins,kappa,Nc,Neq)

#############################################################################
#############################################################################

#### COOLING - GODEC
print('\n **********\n ********** \n')
print('INVERSO')

divIV,eDivIV,divIVr,eDivIVr,energIV,eEnergIV,energIVr,eEnergIVr,SIV,eSIV,SIVr,eSIVr =\
    medias_godec(dir_name,'inverso',fs,ncycles,num_files,dttemp,ntrans,margen,nbins,kappa,Nc,Neq)

#############################################################################
#############################################################################
#############################################################################
#############################################################################

### DIBUJOS GODEC
# Directo
ts = np.arange(ntrans)/fs
n = 60

# KLD
dibujos_godec(divDV,divIV,eDivDV,eDivIV,ts,"Kullback-Leibler divergence","ordinary","KLD",ntrans,ymax)
dibujos_godec(divDVr,divIVr,eDivDVr,eDivIVr,ts,"Kullback-Leibler divergence","reciprocal","KLD",ntrans,ymax)

# Energia

#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################

ntrans_cum = 200
ts = np.arange(ntrans)/fs

media,mediaR,varianza,varianzaR,cum3,cum3R,cum4,cum4R,mediaN,mediaNR,varianzaN,varianzaNR =\
    cum_posicion(dir_name,'directo',fs,ncycles,num_files,dttemp,ntrans_cum,margen)

plt.figure()
plt.plot(ts,varianza)

#############################################################################
#############################################################################
############################## COTA RELAX ###################################
#############################################################################
#############################################################################

## Bath
sigmaf = np.mean(varianza[80:])

dtSbath = -(0.25*fs/sigmaf)*(varianza[2:]-varianza[:-2])
ts = (1 + np.arange(len(dtSbath)))/fs

Sbath = np.zeros(len(dtSbath))
for i in range(len(dtSbath)):
    Sbath[i] = np.sum(dtSbath[0:i])/fs
    
plt.plot(ts,Sbath)

## System
Ssys = SDV - SDV[0]
Ssys = Ssys[1:-1]

## Total
entropy = Ssys + Sbath[0:ntrans-2]
plt.figure()
plt.plot(ts[0:ntrans-2],entropy)

#### Divergencia inversa
divD_inv,eDivD_inv,divDr_inv,eDivDr_inv =\
kld_cota_relax(dir_name,'directo',fs,ncycles,num_files,dttemp,ntrans,margen,nbins,kappa,Nc,Neq)

plt.figure()
plt.plot(ts[0:ntrans-2],entropy,'.')
plt.plot(ts[0:ntrans-2],divD_inv[0:ntrans-2],'.')
plt.legend(['Entropia','KLD inversa'])

## Max vel relax
entropy_mean = entropy/ts[0:ntrans-2]
plt.figure()
plt.plot(ts[0:ntrans-2],ts[0:ntrans-2],'--')
plt.plot(ts[0:ntrans-2],divD_inv[0:ntrans-2]/entropy_mean,'.')


# #### PRINTEAR SERIES DE DIV
# fich = open(dir_name + 'kld_evol_ordinario-means.txt','w')
# fich.write('Time (s)\t\tDiv_heat\t\tDiv_cool \n')
    
# for i in range(len(divDV)):
#     fich.write(str(ts[i]) + '     \t\t' + str(divDV[i]) + '\t' + \
#                 str(eDivDV[i]) + '\t' +  str(divIV[i]) + '\t' + str(eDivIV[i]))
#     fich.write('\n')

# fich.close()

# ####

# fich = open(dir_name + 'kld_evol_reciproco-means.txt','w')
# fich.write('Time (s)\t\tDiv_heat\t\tDiv_cool \n')
    
# for i in range(len(divDVr)):
#     fich.write(str(ts[i]) + '     \t\t' + str(divDVr[i]) + '\t' + \
#                 str(eDivDVr[i]) + '\t' +  str(divIVr[i]) + '\t' + str(eDivIVr[i]))
#     fich.write('\n')

# fich.close()


# ####



# fich = open(dir_name + 'technical_data-means.txt','w')

# fich.write('Stiffness (S.I.) = ' + str(kappa) + '\n'  + '\n')

# fich.write('COLD-WARM PROCESS' + '\n')
# fich.write('Relative temp. (K) = ' + str(tempRelD) + ' +- ' + str(eTempRelD) + '\n')
# fich.write('Experimental KLd (a.u.) = '+str(divD)  + '+-' + str(eDivD) + '\n')
# fich.write('KLd from relative temp. (a.u.) = ' + str(KLD_dir)  + '+-' + str(eKLD_dir) + '\n')

# fich.write('\n')

# fich.write('HOT-WARM PROCESS' + '\n')
# fich.write('Relative temp. (K) = ' + str(tempRelI) + ' +- ' + str(eTempRelI) + '\n')
# fich.write('Experimental KLd (a.u.) = '+str(divI)  + '+-' + str(eDivI) + '\n')
# fich.write('KLd from relative temp. (a.u.) = ' + str(KLD_inv) + '+-' + str(eKLD_inv) + '\n')

# fich.close()