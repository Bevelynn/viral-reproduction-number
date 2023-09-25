## Description of .py files

### "infection_rate.py":  
Calculates the burst size and reproduction number distributions for the model with $\nu_E=\nu_I=0$, under the assumption of a constant number of uninfected target cells (Case 1) and for the case of a depleting target cell population due to infection (Case 2). The reproduction number distributions for Case 1 and Case 2 are plotted for some key values of the infection rate, $\beta$, and the number of target cells, $T_0$, creating Figure 3. The Hellinger distance is calculated between the distributions resulting from the Case 1 and Case 2 methods, for a range of values of $\beta$ and $T_0$. These are presented in a heatmap, creating Figure 4.

### "immune_killing.py":  
The burst size distribution, reproduction number distribution, and probability of viral extinction, are plotted for the budding and bursting models, for different values of the immune clearance rate, $\nu_I$. This reproduces Figures 5-7.

### "infectious_stages.py":  
The burst size and reproduction number distributions are plotted for the budding and bursting models, for different numbers of infectious stages, $n_I$. The probability of viral extinction is also calculated. This reproduces Figures 8 and 9.
