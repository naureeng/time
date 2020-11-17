"""

02.07.2020: first iteration of code


"""
# PLOT_STIMULUS_SWITCHING_TASK.PY plots baseline and change states of a single trial for switching task variant

# step 1: import dependencies
import numpy as np
import matplotlib
matplotlib.use('PDF')


import matplotlib.pyplot as plt
plt.rc('font', family='Arial')
import seaborn as sns

# step 2: create vector of TF changes [Hz]
mus = np.log2([0.25, 0.5, 1, 2, 4]) # mean of change - baseline [TF in Hz]
sigma = 0.25 # standard deviation of stimulus distribution [Hz]
colors = [[153, 153, 153], # color scheme of blues
          [0, 255, 255],
          [0, 204, 255],
          [0, 102, 255],
          [0, 0, 255],
          [0, 0, 153]]

# step 3: create stimulus distribution based on mean [mus] + standard deviation [sigma]
def make_stimulus(change=0):
    nsamples = 100 # 100 is length of output 
    ntrials = 1;
    xvec = np.random.normal(size=(nsamples, ntrials)) * sigma; 
    xvec[-40:,:] += change
    return xvec

# step 4: plot baseline state
t = np.linspace(-2.95, 2.0, num=100) # a vector spaced from -2.95 to 2 with 100 samples
plt.figure(figsize=(3,2)) # xlim = [3,2] to match [-2.95, 2]
plt.plot([t[0], t[-1]], [1,1], color='k', linestyle=':') # plot from t[0] = -2.95 to t[-1] = 2 and scale both up by [1,1] so baseline occurs at 2 Hz

# step 5: plot change states
for mu, c in zip(mus, colors):
    plt.step(t, make_stimulus(change=mu)+1, color=np.array(c)/255, linewidth=0.5)

# step 6: clean up plot 
plt.ylim([-2., 4.])
plt.xlim([-3., 2.])
plt.yticks(ticks=np.arange(-2., 4.), labels=2**np.arange(-2., 4.))
plt.xticks(ticks=np.arange(-3., 3.))
sns.despine(trim=True, offset=10)
plt.ylabel('TF [Hz]')
plt.xlabel('Time to change [s]')
plt.tight_layout()
plt.savefig('figures/stimuli.pdf')
