
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
nuI = 1.6

c = 7

beta = 10**-4
T0=10**5
c = 7
theta = beta*T0/(beta*T0 + c)

start=0

'''
###########################################################################
functions to calculate the burst size and reproduction number distributions
###########################################################################
'''

def burst_size_prob(p, nI, b, bursting):
    '''calculates the probability of burst size b'''
    #if using the model of viral release by budding, set bursting = False
    #if using the model of viral release by bursting, set bursting = True
    
    deltI = nI/tauI
    if bursting == False:
        prob = sum([binom_coef(k+b-1,b)*(p/(p+deltI+nuI))**b*deltI**(k-1)*nuI/(p+deltI+nuI)**k for k in range(1, nI)])
        prob += binom_coef(nI+b-1,b)*(p/(p+deltI+nuI))**b*deltI**(nI-1)*(deltI+nuI)/(p+deltI+nuI)**nI
    elif b==0:
        prob = 1-(deltI/(deltI+nuI))**nI*(1-((deltI+nuI)/(p+deltI+nuI))**nI)
    else:
        prob = binom_coef(nI+b-1,b)*(p/(p+deltI+nuI))**b*(deltI/(p+deltI+nuI))**nI
    return prob

def pmf_B(p, nI, max_b, bursting):
    '''gives the burst size distribution for values of b = 0, ..., max_b'''
    return np.array([burst_size_prob(p, nI, b, bursting) for b in range(max_b+1)])



def pmf_R(nI, max_r, bursting):
    '''gives the reproduction number distribution for values of r = 0, ..., max_r'''
    #assuming a constant number of T0 uninfected target cells
    return pmf_B(theta*p, nI, max_r, bursting)

def mean_R(nI, bursting):
    '''calculates the mean reproduction number'''
    deltI = nI/tauI
    if nuI == 0:
        mean = theta*p*tauI
    elif bursting == False:
        mean = (1/nuI)*theta*p*(1-(deltI/(deltI+nuI))**nI)
    else:
        mean = tauI*theta*p*(deltI/(deltI+nuI))**(nI+1)
    return mean

'''
##########################################################
plotting the burst size distribution for the budding model
##########################################################
'''
    
#the chosen values of nu to plot the results for
nIs=[1, 2, 5, 10, 20]
nI_labels = ['1', '2', '5', '10', '20']

max_b=2000
js=range(max_b+1)

plt.figure()
for nI in nIs:
    #calculate the burst size distribution
    B_dist=pmf_B(p, nI, max_b, False)
    print('Prob(B < max_b)=',sum(B_dist))
    print('mean burst size =', mean_R(nI, False)/theta)
    plt.fill_between(js, [0]*len(js), B_dist, alpha=0.4, label=r'$n_I$='+str(nI))
plt.xlabel('b (burst size)')
plt.ylabel(r'$\mathbb{P}(B = b)$')
plt.legend()
plt.savefig('B_dist_infectious_stages_budding.png', bbox_inches='tight')


    
'''
###################################################################
plotting the reproduction number distribution for the budding model
###################################################################
'''
    
max_r=int(theta*max_b)
Rjs=range(max_r+1)

plt.figure()
for nI in nIs:
    #calculate the reproduction number distribution
    R_dist=pmf_R(nI, max_r, False)
    print('Prob(R < max_r)=',sum(R_dist))
    print('mean reproduction number =', mean_R(nI, False))
    plt.fill_between(Rjs, [0]*len(Rjs), R_dist, alpha=0.4, label=r'$n_I$='+str(nI))
plt.xlabel('r (secondary infections)')
plt.ylabel(r'$\mathbb{P}(R = r)$')
plt.legend()
plt.savefig('R_dist_infectious_stages_budding.png', bbox_inches='tight')




'''
###########################################################
plotting the burst size distribution for the bursting model
###########################################################
'''

f1, ax1 = plt.subplots(1,2,figsize=(20,5))
plt.subplots_adjust(left=None, bottom=None, right=None, top=1.24, wspace=0.4, hspace=0.4)

plt.subplot(1, 2, 2)
#create a list to store the probabilities of zero burst size for different values of nI
probs_0=[]
for nI in nIs:
    #calculate the burst size distribution
    B_dist=pmf_B(p, nI, max_b, True)
    print('Prob(B < max_b)=',sum(B_dist))
    print('mean burst size =', mean_R(nI, True)/theta)
    probs_0.append(B_dist[0])
    plt.fill_between(js[1:], [0]*(len(js)-1), B_dist[1:]/(1-B_dist[0]), alpha=0.4, label=r'$n_I$='+str(nI))
plt.xlabel('b (burst size)')
plt.ylabel(r'$\mathbb{P}(B = b | b \neq 0)$')
plt.legend()

plt.subplot(1,2,1)
for i, nI in enumerate(nIs):
    plt.bar(i, probs_0[i], alpha=0.4)
plt.xlabel(r'$n_I$')
plt.ylabel('$\mathbb{P}(B=0)$')
plt.xticks(range(len(nIs)), nI_labels)
plt.savefig('B_dist_infectious_stages_bursting.png', bbox_inches='tight')


    
'''
####################################################################
plotting the reproduction number distribution for the bursting model
####################################################################
'''

f1, ax1 = plt.subplots(1,2,figsize=(20,5))
plt.subplots_adjust(left=None, bottom=None, right=None, top=1.24, wspace=0.4, hspace=0.4)

plt.subplot(1, 2, 2)
#create a list to store the probabilities of zero secondary infections for different values of nI
probs_0R=[]
for nI in nIs:
    #calculate the reproduction number distribution
    R_dist=pmf_R(nI, max_r, True)
    print('Prob(R < max_r)=',sum(R_dist))
    print('mean reproduction number =', mean_R(nI, True))
    probs_0R.append(R_dist[0])
    plt.fill_between(Rjs[1:], [0]*(len(Rjs)-1), R_dist[1:]/(1-R_dist[0]), alpha=0.4, label=r'$n_I$='+str(nI))
plt.xlabel('r (secondary infections)')
plt.ylabel(r'$\mathbb{P}(R = r | r \neq 0)$')
plt.legend()

plt.subplot(1,2,1)
for i, nI in enumerate(nIs):
    plt.bar(i, probs_0R[i], alpha=0.4)
plt.xlabel(r'$n_I$')
plt.ylabel('$\mathbb{P}(R=0)$')
plt.xticks(range(len(nIs)), nI_labels)
plt.savefig('R_dist_infectious_stages_bursting.png', bbox_inches='tight')


'''
###########################################################
functions to calculate the probability of viral extinction
###########################################################
'''

def pgf_nb(s, k, prob):
    #pgf of a negative binomal distrbution with parameters k and prob
    return ((1-prob)/(1-prob*s))**k

def pgf(s, nI, theta, bursting):
    #pgf of reproduction number distribution
    #called \pi(s) in paper
    if bursting == False:
        #budding case
        #Eq. (2.13) in paper with \nu_E=0
        prob=tauI*theta*p/(tauI*theta*p + nI + tauI*nuI)
        first_sum = sum([(nI/(nI+tauI*nuI))**(k-1)*(tauI*nuI/(nI+tauI*nuI))*pgf_nb(s, k, prob) for k in range(1, nI)])
        return first_sum + (nI/(nI+tauI*nuI))**(nI-1)*pgf_nb(s, nI, prob)
    else:
        #bursting case
        #Eq. (2.14) in paper with \nu_E=0
        return 1 - (nI/(nI+tauI*nuI))**(nI)*(1 - ((nI +tauI*nuI)/(nI + tauI*nuI + (1-s)*tauI*theta*p))**nI)
 
def roots(s, nI, theta, bursting):
    #function to obtain s-pi(s)
    #we want to find roots of this function, i,e., the fixed points of the pgf,
    #since the smallest root gives the probability that a population
    #of infected cells starting with 1 infected cell will go extinct
    return s-pgf(s, nI, theta, bursting)

def prob_i(nI, theta, i, bursting):
    #probability of virus extinction given you start with i infected cells
    #i.e., probability that all i will go extinct
    return (fsolve(roots, start, args=(nI, theta, bursting)))**i

'''
#####################################################
plots for probability of extinction in bursting case
#####################################################
'''

nIs = [1, 2, 5, 10, 20]
Is=np.logspace(0,2,100)
theta_vals = np.logspace(-2.7,0,500)

f1, ax1 = plt.subplots(1,2,figsize=(21,7))
plt.subplots_adjust(left=None, bottom=None, right=None, top=1.24, wspace=None, hspace=None)

probs = []
labels = []
for i, nI in enumerate(nIs):
    R_mean = mean_R(nI, True)
    probs.append(prob_i(nI, theta, 1, True))
    labels.append('$n_I$='+str(nI)+'\n'+r' $\bar R=$'+str(round(R_mean,2)))
    
plt.subplot(122)
for i, nI in enumerate(nIs):
    plt.bar(i, probs[i], alpha=0.4)
plt.title(r'$\theta=$'+str(round(theta,2)))
plt.ylabel('Probability of extinction from one infected cell')
plt.xticks(range(len(nIs)), labels)

plt.subplot(121)
for nI in [1,2,5,10,20]:
    B_mean = mean_R(nI, True)/theta
    probs=np.array([fsolve(roots, start, args=(nI, theta, True)) for theta in theta_vals])
    plt.plot(theta_vals, probs,label='$n_I$='+str(nI)+r', $\mathbb{E}[B]=$'+str(round(B_mean,2)))
theta=beta*T0/(beta*T0+c)
plt.plot([theta, theta], [0.6, 1], '--')
plt.gca().set_xscale("log")
plt.legend(loc=3)
plt.xlabel(r'$\theta$')
plt.ylabel('Probability of extinction from one infected cell')
plt.savefig('prob_extinction_bursting.png', bbox_inches='tight')
