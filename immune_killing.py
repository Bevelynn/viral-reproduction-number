
from scipy.special import binom as binom_coef
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve

plt.rcParams.update({'font.size': 18})
plt.rcParams.update({'axes.labelsize': 'large'})

'''
################
parameter values
################
'''

p = 1000

tauI = 1.25
nI = 10
deltI = nI/tauI

beta = 10**-4
T0=10**5
c = 7
theta = beta*T0/(beta*T0 + c)

start=0

'''
############################################################################
functions to calculate the burst size and reproduction number distributions
############################################################################
'''
#using the model of viral release by budding

def burst_size_prob(p, nuI, b, bursting):
    '''calculates the probability of burst size b'''
    #if using the model of viral release by budding, set bursting = False
    #if using the model of viral release by bursting, set bursting = True
    
    if bursting == False:
        #probability of burst size b
        prob = sum([binom_coef(k+b-1,b)*(p/(p+deltI+nuI))**b*deltI**(k-1)*nuI/(p+deltI+nuI)**k for k in range(1, nI)])
        prob += binom_coef(nI+b-1,b)*(p/(p+deltI+nuI))**b*deltI**(nI-1)*(deltI+nuI)/(p+deltI+nuI)**nI
    else:
        if b==0:
            prob = 1-(deltI/(deltI+nuI))**nI*(1-((deltI+nuI)/(p+deltI+nuI))**nI)
        else:
            prob = binom_coef(nI+b-1,b)*(p/(p+deltI+nuI))**b*(deltI/(p+deltI+nuI))**nI
    return prob

def pmf_B(p, nuI, max_b, bursting):
    '''gives the burst size distribution for values of b = 0, ..., max_b'''
    return np.array([burst_size_prob(p, nuI, b, bursting) for b in range(max_b+1)])



def pmf_R(nuI, max_r, bursting):
    '''gives the reproduction number distribution for values of r = 0, ..., max_r'''
    #assuming a constant number of T0 uninfected target cells
    return pmf_B(theta*p, nuI, max_r, bursting)

def mean_R(nuI, bursting):
    '''calculates the mean reproduction number'''
    if nuI == 0:
        mean = theta*p*tauI
    elif bursting == False:
        mean = (1/nuI)*theta*p*(1-(deltI/(deltI+nuI))**nI)
    else:
        mean = tauI*theta*p*(deltI/(deltI+nuI))**(nI+1)
    return mean
 
    
'''
#####################################
plotting the burst size distribution
#####################################
'''
    
#the chosen values of nu to plot the results for
nus=[0, 0.25, 0.5, 1, 1.6]
nu_labels = ['0', '0.25', '0.5', '1', '1.6']

max_b=2500
js=range(max_b+1)


#budding case
plt.figure()
for nu in nus:
    #calculate the burst size distribution
    B_dist=pmf_B(p, nu, max_b, False)
    print('Prob(B < max_b)=',sum(B_dist))
    print('mean burst size =', mean_R(nu, False)/theta)
    plt.fill_between(js, [0]*len(js), B_dist, alpha=0.4, label=r'$\nu$='+str(nu))
plt.xlabel('b (burst size)')
plt.ylabel(r'$\mathbb{P}(B = b)$')
plt.legend(fontsize = 15)
plt.savefig('B_dist_immune_killing_budding.png', bbox_inches='tight')


#bursting case
f1, ax1 = plt.subplots(1,2,figsize=(20,5))
plt.subplots_adjust(left=None, bottom=None, right=None, top=1.24, wspace=0.4, hspace=0.4)

plt.subplot(1, 2, 2)
#create a list to store the probabilities of zero burst size for different values of nuI
probs_0=[]
for nuI in nus:
    #calculate the burst size distribution
    B_dist=pmf_B(p, nuI, max_b, True)
    print('Prob(B < max_b)=',sum(B_dist))
    print('mean burst size =', mean_R(nuI, True)/theta)
    probs_0.append(B_dist[0])
    plt.fill_between(js[1:], [0]*(len(js)-1), B_dist[1:]/(1-B_dist[0]), alpha=0.4, label=r'$\nu_I$='+str(nuI))
plt.xlabel('b (burst size)')
plt.ylabel(r'$\mathbb{P}(B = b | b \neq 0)$')
plt.legend()

plt.subplot(1,2,1)
for i, nuI in enumerate(nus):
    plt.bar(i, probs_0[i], alpha=0.4)
plt.xlabel(r'$\nu_I$')
plt.ylabel('$\mathbb{P}(B=0)$')
plt.xticks(range(len(nus)), nu_labels)
plt.savefig('B_dist_immune_killing_bursting.png', bbox_inches='tight')


    
'''
#############################################
plotting the reproduction number distribution
#############################################
'''
    
max_r=int(theta*max_b)
Rjs=range(max_r+1)

#budding case
plt.figure()
for nu in nus:
    #calculate the reproduction number distribution
    R_dist=pmf_R(nu, max_r, False)
    print('Prob(R < max_r)=',sum(R_dist))
    print('mean reproduction number =', mean_R(nu, False))
    plt.fill_between(Rjs, [0]*len(Rjs), R_dist, alpha=0.4, label=r'$\nu_I$='+str(nu))
plt.xlabel('r (secondary infections)')
plt.ylabel(r'$\mathbb{P}(R = r)$')
plt.legend(fontsize = 15)
plt.savefig('R_dist_immune_killing_budding.png', bbox_inches='tight')

#bursting case
f1, ax1 = plt.subplots(1,2,figsize=(20,5))
plt.subplots_adjust(left=None, bottom=None, right=None, top=1.24, wspace=0.4, hspace=0.4)

plt.subplot(1, 2, 2)
#create a list to store the probabilities of zero secondary infections for different values of nI
probs_0R=[]
for nuI in nus:
    #calculate the reproduction number distribution
    R_dist=pmf_R(nuI, max_r, True)
    print('Prob(R < max_r)=',sum(R_dist))
    print('mean reproduction number =', mean_R(nuI, True))
    probs_0R.append(R_dist[0])
    plt.fill_between(Rjs[1:], [0]*(len(Rjs)-1), R_dist[1:]/(1-R_dist[0]), alpha=0.4, label=r'$\nu_I$='+str(nuI))
plt.xlabel('r (secondary infections)')
plt.ylabel(r'$\mathbb{P}(R = r | r \neq 0)$')
plt.legend()

plt.subplot(1,2,1)
for i, nuI in enumerate(nus):
    plt.bar(i, probs_0R[i], alpha=0.4)
plt.xlabel(r'$\nu_I$')
plt.ylabel('$\mathbb{P}(R=0)$')
plt.xticks(range(len(nus)), nu_labels)
plt.savefig('R_dist_immune_killing_bursting.png', bbox_inches='tight')


'''
###########################################################
functions to calculate the probability of viral extinction
###########################################################
'''

def pgf_nb(s, k, prob):
    #pgf of a negative binomal distrbution with parameters k and prob
    return ((1-prob)/(1-prob*s))**k

def pgf(s, nu, bursting):
    #pgf of reproduction number distribution
    #called \pi(s) in paper
    if bursting == False:
        #budding case
        #Eq. (2.14) in paper with \nu_E=0
        prob=tauI*theta*p/(tauI*theta*p + nI + tauI*nu)
        first_sum = sum([(nI/(nI+tauI*nu))**(k-1)*(tauI*nu/(nI+tauI*nu))*pgf_nb(s, k, prob) for k in range(1, nI)])
        return first_sum + (nI/(nI+tauI*nu))**(nI-1)*pgf_nb(s, nI, prob)
    else:
        #bursting case
        #Eq. (2.15) in paper with \nu_E=0
        return 1 - (nI/(nI + tauI*nu))**(nI)*(1 - ((nI + tauI*nu)/(nI + tauI*nu + (1-s)*tauI*theta*p))**nI)
 
def roots(s, nu, bursting):
    #function to obtain s-pi(s)
    #we want to find roots of this function, i,e., the fixed points of the pgf,
    #since the smallest root gives the probability that a population
    #of infected cells starting with 1 infected cell will go extinct
    return s-pgf(s, nu, bursting)

def prob_i(nu, i, bursting):
    #probability of virus extinction given you start with i infected cells
    #i.e., probability that all i will go extinct
    return (fsolve(roots, start, args=(nu, bursting)))**i

'''
###################################
plot for probability of extinction
###################################
'''
Is=np.logspace(0,2,100)

#budding case
plt.figure()
for nu in nus:
    probs=np.array([prob_i(nu, I, False) for I in Is])
    plt.plot(Is,probs,label=r'$\nu_I$='+str(nu))
plt.gca().set_xscale("log")
plt.legend(loc=0)
plt.xlabel('Initial number of infected cells')
plt.ylabel('Probability of extinction')

#bursting case
plt.figure()
for nu in nus:
    probs=np.array([prob_i(nu, I, True) for I in Is])
    print('prob extinction starting from one infected cell =', probs[0])
    plt.plot(Is,probs,label=r'$\nu_I$='+str(nu))
plt.gca().set_xscale("log")
plt.legend(loc=0)
plt.xlabel('Initial number of infected cells')
plt.ylabel('Probability of extinction')
plt.savefig('prob_extinction_immune_killing_bursting.png', bbox_inches='tight')