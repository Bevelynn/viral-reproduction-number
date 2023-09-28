
from scipy.special import binom as binom_coef
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


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

c = 7

'''
###########################################################################
functions to calculate the burst size and reproduction number distributions
###########################################################################
'''

def burst_size_prob(p, b):
    '''calculates the probability of burst size b'''
    prob = binom_coef(nI+b-1,b)*(p/(p+deltI))**b*deltI**nI/(p+deltI)**nI
    return prob

def pmf_B(p, max_b):
    '''gives the burst size distribution for values of b = 0, ..., max_b'''
    return [burst_size_prob(p, b) for b in range(max_b+1)]



def pmf_R_case1(beta, T0, max_r):
    '''gives the reproduction number distribution for values of r = 0, ..., max_r'''
    '''Case 1 (i.e. assuming a constant number of T0 uninfected target cells)'''
    theta = beta*T0/(beta*T0 + c)
    return pmf_B(theta*p, max_r)



def P(beta, T0, max_b, max_r):
    '''returns a matrix of elements p_r[b], which is the probability of r secondary infections, given burst size b'''
    #the index of the row corresponds to the value of b
    #the index of the column corresponds to the value of r
    #the matrix elements will be zero for r > min(b, T0)
    #T0 is the initial number of target cells in the population
    
    #initialize a matrix with zeros to fill in
    N=np.zeros((max_b+1,max_r+1))
    N[:,0]=[(c/(beta*T0+c))**b for b in range(max_b+1)]
    for b in range(1,max_b+1):
        for r in range(1,max_r+1):
            if r <= min(b,T0):
                N[b,r]=N[b-1,r-1]*beta*(T0+1-r)/(beta*(T0+1-r)+c)+N[b-1,r]*c/(beta*(T0-r)+c)
    return N

def pmf_R_case2(beta, T0, max_b, max_r):
    '''gives the probability mass function of the reproduction number distribution for Case 2'''
    #burst size distribution with nuE = nuI = 0
    burst_size_dist = pmf_B(p, max_b)
    mat = P(beta, T0, max_b, max_r)
    pmf = []
    for r in range(max_r+1):
        probability=sum([mat[b,r]*burst_size_dist[b] for b in range(max_b+1)])
        pmf.append(probability)
    return pmf



def hellinger(P,Q):
     '''calculates the hellinger distance between two probability distributions'''
     sum=0
     for k in range(len(P)):
         sum+=(np.sqrt(P[k])-np.sqrt(Q[k]))**2
     dist=np.sqrt(sum)/np.sqrt(2)
     return dist

'''
#####################################################
examples of differences between case 1 and case 2
#####################################################
'''

max_b=2500

betas=[10**-2, 10**-4, 10**-2, 10**-4]
beta_labels=['$10^{-2}$','$10^{-4}$', '$10^{-2}$', '$10^{-4}$']

Ts=[10, 1000, 1000, 100000]
T_labels=['10','$10^3$', '$10^3$', '$10^5$']

plt.figure(figsize=(12,8))
plt.subplots_adjust(left=None, bottom=None, right=None, top=1.24, wspace=0.4, hspace=0.4)

for i,T0 in enumerate(Ts):
    beta=betas[i]
    
    theta=beta*T0/(beta*T0+c)
    print('theta=', theta)
    
    mean_case1=tauI*theta*p
    print('Rbar=',mean_case1)
    
    max_r=int(2*mean_case1+5)
    js=range(max_r+1)
    
    case1=pmf_R_case1(beta, T0, max_r)
    #print('sum of case 1 distribution =', sum(case1))
    
    case2=pmf_R_case2(beta, T0, max_b, max_r)
    #print('sum of case 2 distribution =', sum(case2))
    
    dist = hellinger(case1,case2)
    print('hellinger distance =', dist)
    
    plt.subplot(2,2,i+1)
    if i==0 or i==1:
        plt.bar(range(max_r+1), case1, alpha=0.5, color='tomato', label='Case 1')
        plt.bar(range(max_r+1)[:T0+1], case2[:T0+1], alpha=0.5, color='forestgreen', label='Case 2')
        plt.ylim(-0.0001,0.3)
    else:
        plt.fill_between(range(max_r+1),[0]*(max_r+1), case1, alpha=0.5, color='tomato', label='Case 1')
        plt.fill_between(range(max_r+1)[:min(max_r,T0)+1],[0]*(min(max_r,T0)+1), case2[:T0+1], alpha=0.5, color='forestgreen', label='Case 2')
        plt.ylim(-0.0001,0.0035)
    plt.xlabel('r (secondary infections)')
    if i==0 or i==2:
        plt.ylabel('$\mathbb{P}(R=r)$')
    plt.title('$T_0$='+T_labels[i]+r', $\beta$='+beta_labels[i])
    plt.legend()
plt.savefig('R_dist_examples_case1_case2.png', bbox_inches='tight')


'''
#######
heatmap
#######
'''

betarange=np.logspace(-5,0,11)
T0range=np.logspace(0,5,11)

#create an empty array to store the Hellinger distances in
dists=np.empty((len(betarange),len(T0range)))

#each column corresponds to a value of T0
#each row corresponds to a value of beta
for i,T0 in enumerate(T0range):
    T0=int(T0)
    for j,beta in enumerate(betarange):
        theta=beta*T0/(beta*T0+c)
        mean_case1=tauI*theta*p
        max_r=int(2*mean_case1+5)
        case1=pmf_R_case1(beta, T0, max_r)
        #print('sum of case 1 dist=',sum(case1))
        case2=pmf_R_case2(beta, T0, max_b, max_r)
        #print('sum of case 2 dist=',sum(case2))
        dist = hellinger(case1,case2)
        dists[j,i]=dist

sns.light_palette("seagreen", as_cmap=True)

x_axis_labels = [0,'',1,'',2,'',3,'',4,'',5] # labels for x-axis
y_axis_labels = [-5,'',-4,'',-3, '', -2, '', -1, '', 0] # labels for y-axis


plt.figure()
ax=sns.heatmap(dists, xticklabels=x_axis_labels, yticklabels=y_axis_labels, cmap=sns.light_palette("seagreen", as_cmap=True), vmin=0, vmax=1)
plt.setp(ax.get_yticklabels(), rotation=360)
ax.invert_yaxis()
ax.collections[0].colorbar.set_label("Hellinger distance")
plt.xlabel('$log_{10}T_0$ (target cells)')
plt.ylabel(r'$log_{10}\beta$ (infection rate)')
plt.savefig('case1_case2_heatmap.png' ,bbox_inches='tight')
