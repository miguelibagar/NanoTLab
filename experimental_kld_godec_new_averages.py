import numpy as np
import matplotlib.pyplot as plt
import statistics as st

##################################
##################################

def trayectoria(filename,fs,ncycles,num_files,dttemp,ntrans,margen,contador):
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
        fid = open(filename + str(cFile+contador) + '.out', 'r') 
        
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
        
        print('Archivo ' + str(contador+1))
        cFile += 1
        
    return eq_series,trans_array,trans_arrayR,low_series,temp1,temp2
    
def histogramas(bins,serie):
    hist,bins_hist = np.histogram(serie,bins=bins,density=True)
    
    return hist,bins_hist

def integrales(distr1,eDistr1,distrRef,eDistrRef,nbins,Dx):
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

def temps(sigma,kappa,nbinsT,title,option):
    kB = 1.38064852E-23   # J/K
    
    temp1D = kappa*sigma/kB

    tempI = np.mean(temp1D)
    eTempI = np.sqrt(st.variance(temp1D))    
    
    # bins1 = np.linspace(temp1D.min()-10,temp1D.max()+10,num=nbinsT)
    
    
    # hist1,bins1 = np.histogram(temp1D,bins=bins1)
    # plt.figure()
    # plt.hist(bins1[:-1], bins1, weights=hist1,density=True)

    # xs = np.linspace(temp1D.min(),temp1D.max(),num=200)
    # plt.plot(xs,np.exp(-np.power(xs-tempI,2)/(2*eTempI**2)) / (np.sqrt(2*3.14159)*eTempI))
    # plt.xlabel("Effective temperature (K)")
    # plt.title("Histogram of " + title + " temperatures (" + option + ")")
    # plt.savefig(dir_name + 'hist_' + title + '_' + option + '-means.png')
    
    return tempI,eTempI

def protocolo(dir_name,opt,fs,ncycles,num_files,dttemp,ntrans,margen,nbins,kappa,Nc,Neq,contador):
    dir_name = dir_name + opt + '/'
    
    eq_series,trans_array,trans_arrayR,low_series,temp1D,temp2D = \
        trayectoria(dir_name,fs,ncycles,num_files,dttemp,ntrans,margen,contador)
        
    bins= np.linspace(eq_series.min()-0.5e-8,eq_series.max()+0.5e-8,num=nbins)
    Dx = bins[1]-bins[0]

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
            integrales(hist_trans,eHist_trans,hist_eq,eHist_eq,nbins,Dx)
        
        divDV[i] = kld
        eDivDV[i] = eKld
        energDV[i] = energ
        eEnergDV[i] = eEnerg
        SDV[i] = S
        eSDV[i] = eS
        
        trans_series = trans_arrayR[:,i]
        hist_trans,bins_trans = histogramas(bins,trans_series)
        
        kld,eKld,energ,eEnerg,S,eS = \
            integrales(hist_trans,eHist_trans,hist_eq,eHist_eq,nbins,Dx)
        
        divDVr[i] = kld
        eDivDVr[i] = eKld
        energDVr[i] = energ
        eEnergDVr[i] = eEnerg
        SDVr[i] = S
        eSDVr[i] = eS
        
        i += 1

    # Histogramas globales
    divD,eDivD,_,_,_,_ = integrales(hist_low,eHist_low,hist_eq,eHist_eq,nbins,Dx)

    #### Temperaturas
    if opt=='directo':
        tempCD,eTempCD = temps(temp1D,kappa,30,'cold','direct')
        tempWD,eTempWD = temps(temp2D,kappa,30,'warm','direct')
    
    if opt=='inverso':
        tempCD,eTempCD = temps(temp1D,kappa,30,'hot','inverse')
        tempWD,eTempWD = temps(temp2D,kappa,30,'warm','inverse')

    tempRelD = tempCD/tempWD
    eTempRelD = np.sqrt(np.power(tempCD/tempWD,2) + \
                        np.power(tempCD*eTempWD/(tempWD**2),2))

    KLD_dir = 0.5*(tempRelD - 1 - np.log(tempRelD))
    eKLD_dir = 0.5*np.abs(1-1/tempRelD)*eTempRelD
            
    return divDV,eDivDV,divDVr,eDivDVr,energDV,eEnergDV,energDVr,\
        eEnergDVr,SDV,eSDV,SDVr,eSDVr,divD,eDivD,tempRelD,eTempRelD,\
            KLD_dir,eKLD_dir,hist_trans,bins_trans

def dibujos(serieHeat,serieCool,eHeat,eCool,ts,magnitude,opt,fname,N):
    plt.figure()
    plt.semilogx(ts[0:N],serieHeat[0:N],'.r')
    plt.semilogx(ts[0:N],serieCool[0:N],'.b')
    plt.legend(['Cold-warm','Hot-warm'])
    plt.xlabel("Time (s)")
    plt.ylabel(magnitude + " (a.u.)")
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
    plt.title("Evolution of " + magnitude + " (" + opt + " process)")
    plt.savefig(dir_name + fname + '-' + opt + 'error-means.png')
    

##############################
##############################

dir_name = 'C:/Users/DELL/Desktop/LAB/23062022/set2/t6 (incompleta)/'

#### Datos generales
fs = 50000  # Hz                # Frecuencia de sampleo
dttemp = 0.005 # s
ncycles = 2399
num_files = 4
margen = 0.003
ntrans = 65
kappa = 41.64e-6 #tanda7
#kappa = 14.04e-6
#kappa = 79.40e-6


Nc = ncycles*num_files
Neq = int(Nc*fs*(dttemp-margen)-21)
nbins = 75


##############################
##############################
        
#### CASO DIRECTO
print('DIRECTO')

divDV = np.zeros(ntrans)
eDivDV = np.zeros(ntrans)
divDVr = np.zeros(ntrans)
eDivDVr = np.zeros(ntrans)

divD = 0
eDivD = 0

contador = 0
while contador < num_files:
    divDVV,eDivDVV,divDVVr,eDivDVVr,energDVV,eEnergDVV,energDVVr,eEnergDVVr,SDVV,\
        eSDVV,SDVVr,eSDVVr,\
        divDD,eDivDD,tempRelD,eTempRelD,KLD_dir,eKLD_dir,hist_trans,bins_trans = \
        protocolo(dir_name,'directo',fs,ncycles,1,dttemp,ntrans,margen,nbins,kappa,Nc,Neq,contador) 
    
    divDV = divDV + divDVV
    eDivDV = eDivDV + eDivDVV
    divDVr = divDVr + divDVVr
    eDivDVr = eDivDVr + eDivDVVr
    
    divD = divD + divDD
    eDivD = eDivD + eDivDD
    
    # plt.figure()
    # plt.hist(bins_trans[:-1], bins_trans, weights=hist_trans,density=True)
    
    contador += 1

divDV = divDV/num_files
eDivDV = eDivDV/num_files
divDVr = divDVr/num_files
eDivDVr = eDivDVr/num_files

divD = divD/num_files
eDivD = eDivD/num_files

print('\n ********** \n')
print('Distancia (directo) = '+str(divD) + '+-' + str(eDivD))

print('Temp. relativa = ' + str(tempRelD) + ' +- ' + str(eTempRelD))

KLD_dir = 0.5*(tempRelD - 1 - np.log(tempRelD))
eKLD_dir = 0.5*np.abs(1-1/tempRelD)*eTempRelD
print('KLD a esas temperaturas (directo) = ' + str(KLD_dir) + '+-' + str(eKLD_dir))


##############################################################
##############################################################
##############################################################

#### CASO INVERSO
print('\n **********\n ********** \n')
print('INVERSO')

divIV = np.zeros(ntrans)
eDivIV = np.zeros(ntrans)
divIVr = np.zeros(ntrans)
eDivIVr = np.zeros(ntrans)

divI = 0
eDivI = 0

contador = 0
while contador < num_files:
    divIVV,eDivIVV,divIVVr,eDivIVVr,energIVV,eEnergIVV,energIVVr,eEnergIVVr,SIVV,\
        eSIVV,SIVVr,eSIVVr,\
        divII,eDivII,tempRelI,eTempRelI,KLD_inv,eKLD_inv,hist_trans,bins_trans = \
        protocolo(dir_name,'inverso',fs,ncycles,1,dttemp,ntrans,margen,nbins,kappa,Nc,Neq,contador) 
    
    divIV = divIV + divIVV
    eDivIV = eDivIV + eDivIVV
    divIVr = divIVr + divIVVr
    eDivIVr = eDivIVr + eDivIVVr
    
    divI = divI + divII
    eDivI = eDivI + eDivII
    
    # plt.figure()
    # plt.hist(bins_trans[:-1], bins_trans, weights=hist_trans,density=True)
    
    contador += 1

divIV = divIV/num_files
eDivIV = eDivIV/num_files
divIVr = divIVr/num_files
eDivIVr = eDivIVr/num_files

divI = divI/num_files
eDivI = np.sqrt(eDivI)/num_files

print('\n ********** \n')
print('Distancia (inverso) = '+str(divI) + '+-' + str(eDivI))

print('Temp. relativa = ' + str(tempRelI) + ' +- ' + str(eTempRelI))

KLD_inv = 0.5*(tempRelI - 1 - np.log(tempRelI))
eKLD_inv = 0.5*np.abs(1-1/tempRelI)*eTempRelI
print('KLD a esas temperaturas (inverso) = ' + str(KLD_inv) + '+-' + str(eKLD_inv))

##############################################################
##############################################################
##############################################################

### DIBUJOS
# Directo
ts = np.arange(ntrans)/fs
n = 60

# KLD
dibujos(divDV,divIV,eDivDV,eDivIV,ts,"Kullback-Leibler divergence","ordinary","KLD",100)
dibujos(divDVr,divIVr,eDivDVr,eDivIVr,ts,"Kullback-Leibler divergence","reciprocal","KLD",100)

# Energia

##############################################################
##############################################################
##############################################################

# #### PRINTEAR SERIES DE DIV
fich = open(dir_name + 'kld_evol_ordinario-means.txt','w')
fich.write('Time (s)\t\tDiv_heat\t\tDiv_cool \n')
    
for i in range(len(divDV)):
    fich.write(str(ts[i]) + '     \t\t' + str(divDV[i]) + '\t' + \
                str(eDivDV[i]) + '\t' +  str(divIV[i]) + '\t' + str(eDivIV[i]))
    fich.write('\n')

fich.close()

####

fich = open(dir_name + 'kld_evol_reciproco-means.txt','w')
fich.write('Time (s)\t\tDiv_heat\t\tDiv_cool \n')
    
for i in range(len(divDVr)):
    fich.write(str(ts[i]) + '     \t\t' + str(divDVr[i]) + '\t' + \
                str(eDivDVr[i]) + '\t' +  str(divIVr[i]) + '\t' + str(eDivIVr[i]))
    fich.write('\n')

fich.close()


####



fich = open(dir_name + 'technical_data-means.txt','w')

fich.write('Stiffness (S.I.) = ' + str(kappa) + '\n'  + '\n')

fich.write('COLD-WARM PROCESS' + '\n')
fich.write('Relative temp. (K) = ' + str(tempRelD) + ' +- ' + str(eTempRelD) + '\n')
fich.write('Experimental KLd (a.u.) = '+str(divD)  + '+-' + str(eDivD) + '\n')
fich.write('KLd from relative temp. (a.u.) = ' + str(KLD_dir)  + '+-' + str(eKLD_dir) + '\n')

fich.write('\n')

fich.write('HOT-WARM PROCESS' + '\n')
fich.write('Relative temp. (K) = ' + str(tempRelI) + ' +- ' + str(eTempRelI) + '\n')
fich.write('Experimental KLd (a.u.) = '+str(divI)  + '+-' + str(eDivI) + '\n')
fich.write('KLd from relative temp. (a.u.) = ' + str(KLD_inv) + '+-' + str(eKLD_inv) + '\n')

fich.close()