# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 13:37:41 2019

@author: Willem
"""



import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm, lognorm, exponweib
from math import pi

#plt.close('all')

num = 4





#--------------------------Beam design-------------------------------------
Fdesign =   60000.0 #71500*9.81 #
Fapp = 0.65*Fdesign # Applied force N
#h	=	0.19	#1.0	Total depth of beam	m
b	=	0.08	#0.6 	Breadth of beam	m
L = 1.0 #11.5#	1.0 #	Length of beam	m
S = 0.2 #2.0	 #0.2 #	Span of applied force	m
cover = 0.025

fyc=40.0*1e6
fys= 450.0*1e6
Es	=	2.08E+11	#	Young's modulus of steel	Pa
Ec	=	2.3e10#33.2E+9	#	Young's modulus of concrete	Pa

centre = 27.7 # [14.0342,27.7]
slope = 8.858 # [3.5,8.858]
means = []

#cert = 0.1 #choose certainty/ies of P-F
plt.figure()
for h  in [0.19,1.0]:
    heff = h- cover
    cert = 0.1 #choose certainty/ies of P-F
    M = Fdesign/2.0*(L-S)/2.0
    print('M=',M)
    K = M/(b*heff**2*fyc)
    print('K=',K)
    if K>0.156:
        print ('Not singly')
    #    exit()
    
    z = heff*(0.5+np.sqrt(0.25-K/0.9))
    print('z=',z)
    Ast_ten = M/(0.87*fys*z)
    #Ast_ten = np.pi*(12.0e-3/2)**2
    Dst_ten = 2.0*np.sqrt(Ast_ten/(num*pi))
#    Dst_ten = 0.0085
    print ('dia of rebars = ', Dst_ten)
    
    #    print 'D_{st}=', Dst_ten
    
    #-----------------------S-N curve of rebar:--------------------------------
    SD_Tilly = 0.0472 #For Sr
#    SD_Tilly = 0.01 #For Sr
    def Tilly_reverse(sig): #in MPa   
        sig = np.log10(sig)
        return centre-slope*sig
    def TillySD_reverse(sig): #in MPa
        sig = np.log10(sig)
        return centre+10**(2*SD_Tilly)-slope*sig
    
    def ppf(x,sig): #Tilly's CDF
        mean = 10**Tilly_reverse(sig)
        sd = (np.log(10**TillySD_reverse(sig))-np.log(mean))/2.0
        if sig>650.0:
            return 1.0
        else:    
            return lognorm.ppf(x,sd,scale=mean)
    
    
    #---------------------------CD calc----------------------------------------
    def cd(n):# crack depth from bottom
        #Adjust steel area
        Aten	=	n*pi*(Dst_ten/2.0)**2	#	Area of tension rebar	m^2
        A_ten	=	Es/Ec*Aten	#	Adjusted tension rebar area	m^2
        #	Find the neutral axis	
        A	=	-b/2	#		
        B	=	A_ten + h*b	#		
        C	=	-(A_ten*(h-heff)+b/2*h**2)	#		
        cd	=	(-B + np.sqrt(B**2-4*A*C))/(2*A)	#	Neutral axis	
        return cd
    #    print "cd for ",num,"rebars =" , cd(num)
    
    #---------------------------stress calc----------------------------------------
    
    def Stresscalc(n): #
        def cd(n):# crack depth from bottom
            #Adjust steel area
            Aten	=	n*pi*(Dst_ten/2.0)**2	#	Area of tension rebar	m^2
            A_ten	=	Es/Ec*Aten	#	Adjusted tension rebar area	m^2
            #	Find the neutral axis	
            A	=	-b/2	#		
            B	=	A_ten + h*b	#		
            C	=	-(A_ten*(h-heff)+b/2*h**2)	#		
            cd	=	(-B + np.sqrt(B**2-4*A*C))/(2*A)	#	Neutral axis	
            return cd
        cd = cd(n)
        #Adjust steel area
        Aten	=	n*pi*(Dst_ten/2)**2	#	Area of tension rebar	m^2
        A_ten	=	Es/Ec*Aten	#	Adjusted tension rebar area	m^2    
        #Find second moment of area		
        y1	=	cd-(h-heff)	#					
        y2	=	(h-cd)/2	#	
        A1	=	A_ten	#		
        A2	=	(h-cd)*b	#			
        I1	=	1.0/12.0*A_ten*Dst_ten**2	#	
        I2	=	1.0/12.0*b*(h-cd)**3		
        Itot	=	I1+y1**2*A1+I2+y2**2*A2     #	Second moment of area	m^4
        #Find stress in tension rebar					
        Mmax	=	Fapp/2.0*(L-S)/2.0	#		Nm
        sig_ten_	=	(Mmax*(y1)/Itot)/10**6	#		Pa
        sig_ten	=	Es/Ec * sig_ten_#		MPa 
        return sig_ten
    #    print 'Stress with num number of rebar:' , Stresscalc(num)
    #    print 'Stress with num-1 number of rebar:' , Stresscalc(num-1)
    
    print ('stress at',num, '=', Stresscalc(num))
    print ('stress at',num-1, '=', Stresscalc(num-1))
    print ('stress at',num-2, '=', Stresscalc(num-2))
    
    #for ploting (delete)
    means.append(Tilly_reverse(Stresscalc(num)))
    
    
    
    #------------------------Weibull--------------------------------------
    
    
    def weib(x,lamb,k):
        return (k / lamb) * (x / lamb)**(k-1) * np.exp(-(x/lamb)**k)
    
    
    def weibcumdist(x,lamb,k): #weibull CDF 
        return 1-np.exp(-(x/lamb)**k)
    
    def weibcumdist_rev(x,lamb,k):
        return lamb*(-np.log(1-x))**(1/k)
    
    
    #----------------------------Monte Carlo---------------------------------
    steps = 5000 
    
    ns= []
    n1s=[]
    nf=[]
    for i in range(num-1):
        ns.append([])
    
    for k in range(steps):
        #Generate num random numbers between 0  &  1
        rand = []        
        for i in range(num):
            rand.append(np.random.random())
        rand = sorted(rand)
    
        N=np.zeros((num,num))    
        for i in range(num):
            for j in range(num):
                N[i][j]=ppf(rand[j],Stresscalc(num-i))
    #                print ppf(rand[num-1],Stresscalc(num-(num-1)))
        n=np.zeros(num) #maak num length vec
        for i in range(num): # 0,1,2,3... num-1
            if i == 0:
                n[0]=N[i][i]
            else:
                som = 0.0
                for j in range(i-1):
                    som = som + (n[j+1]-n[j])/N[j][i]
                n[i]=n[i-1]+(1.0-n[0]/N[0][i]-som)*N[i][i]
        n1s.append(n[0]) 
        nf.append(n[num-1])
        for i in range(num-1):
            ns[i].append(n[i+1]-n[i]) #n_k+1 -n_k
    
    print(ppf(0.5,Stresscalc(num)))
    
    #plt.figure()
    #plt.hist(nf, bins=100)
    #plt.xscale('log')
    
    
    
    
    #----------------------------Fit Weibull over n_k+1 - n_k----------------
    
    for i in range(num-1): 
        weibull_params = exponweib.fit(ns[i], floc=0, f0=1)
        k_shape, lamb_scale = weibull_params[1], weibull_params[3]
    #    
        #--------plot Weibul fit----------------
#        plt.figure()
#        ax = plt.subplot(111)
#        plt.ylabel('Number of failures')
#        ax2 = plt.twinx()
#        ax.hist(ns[i], bins=100)
#        shape,loc,scale = lognorm.fit(ns[i])
#        xlog = np.logspace(0, 6, 500)
#        xweib = np.linspace(0,10**6,500)  
#        ax2.plot(xweib, weib(xweib, lamb_scale, k_shape),'g')
#        plt.grid()
#        plt.xlabel('N')
#        plt.ylabel('Probability desity')
#        plt.title('Weibull probability density fuction of $n_{' +str(i+2)+ '}-n_{'+str(i+1)+'}$')
#        
#        #--------plot CPD of Weibul----------------
#        plt.figure()
#        plt.plot(xweib,weibcumdist(xweib,lamb_scale,k_shape))
#        plt.grid()
#        plt.xlabel('N')
#        plt.ylabel('Probability of failure')
#        plt.title('Cumulative distribution function of $n_{' +str(i+2)+ '}-n_{'+str(i+1)+'}$')
    
    
    
    #print FailureRate(100,lamb_scale,k_shape)
    #8/1000000 kans om by 1 cycle te breek
    #8/1000000/100 kans om die eerste dag te breek (assume 100 cycles per dag)
    
    
    #--------------plot Failure curve---------------
          
    plotx=[]
    plotx.append(0)
    plotx.append(0)
    
    nplot =  20000#220000#np.average(n1s) #cycles to failure of 1st rebar
    plotx.append(nplot)
    plotx.append(nplot)
    for i in range(num-1):
        nplot = nplot + \
        weibcumdist_rev(cert,exponweib.fit(ns[i], floc=0, f0=1)[3],exponweib.fit(ns[i], floc=0, f0=1)[1])
        
        plotx.append(nplot)
        plotx.append(nplot)
    ploty=[]
    for i in range(0,num):
        i = i + 1
        ploty.append(cd(i))
        ploty.append(cd(i))
    ploty.append(0)
    ploty = ploty[::-1]
    ploty.append(h)
    plt.plot(plotx,ploty,linewidth=2)

    plt.legend()
    plt.axis([0,max(plotx)*1.1,0,h*1.05])
    plt.xlabel('Cycles')
    plt.ylabel('Crack depth (CD) [m]')
     
    #plt.savefig('C:\\Users\Willem\\Dropbox\\Meesters\\Analytical Model and Calcs\\Exercise model\\'+'Fail'+str(num)+'_'+ str(h) + '_' + str(b)+'.png')
        
    
    
#    #--------------plot P-F curve---------------
#    plt.figure()
#    
#    plotx=[]
#    plotx.append(0)
#    plotx.append(0)
#    
#    nplot =1#ppf(cert,Stresscalc(num)) #np.average(n1s) #220000#ppf(cert,Stresscalc(num)) #cycles to failure of 1st rebar
#    plotx.append(nplot)
#    plotx.append(nplot)
#    for i in range(num-1):
#        nplot = nplot + \
#        weibcumdist_rev(cert,exponweib.fit(ns[i], floc=0, f0=1)[3],exponweib.fit(ns[i], floc=0, f0=1)[1])
#        
#        plotx.append(nplot)
#        plotx.append(nplot)
#    ploty=[]
#    for i in range(0,num):
#        i = i + 1
#        ploty.append(1-cd(i)/h+cd(num)/h)
#        ploty.append(1-cd(i)/h+cd(num)/h)
#    ploty.append(1)
#    ploty = ploty[::-1]
#    ploty.append(0)
#    pflabel='Certainty of Failure ='+ str(cert)
#    plt.plot(plotx,ploty,linewidth=2,label=pflabel)
#    plt.grid()
#    
#    plt.axis([0,max(plotx)*1.1,0,1.3])
#    plt.xlabel('Cycles')
#    plt.ylabel('Condition')
#    #plt.legend()
#    title= 'P-F interval for '+str(num)+' rebars'
#    plt.title(title)
#    #    ##plt.savefig('C:\\Users\Willem\\Dropbox\\Meesters\\Analytical Model and Calcs\\Exercise model\\'+'PF'+str(num)+'_'+ str(h) + '_' + str(b)+'.png')
#    plt.legend()
#    plt.grid()   
#    plt.axis([0,max(plotx)*1.1,0,h*1.05])
#    plt.xlabel('Cycles')
#    plt.ylabel('Crack depth (CD) [m]')
#    title= ' H = '+ str(h) + 'm & B = ' + str(b)+'m \n Certainty = 90% & ' +str(num)+' rebars \n '#$\sigma_{yc}$ = 30, 40 MPa & $\sigma_{ys}$ = 450 MPa'
#    plt.title(title) 
plt.grid()  
plt.show()
##