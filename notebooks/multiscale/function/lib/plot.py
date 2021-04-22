"""
For plotting.
"""
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy.stats import pearsonr

from lib.calcium import get_section

def view_trials(trace, condition, show_lines=True):
    
  plt.figure(figsize=(16,16))

  angle_list = np.unique(condition)
  for i in range(16):

    ax = plt.subplot(4,4,i+1)

    angle = angle_list[i]
    section_list = get_section(condition, angle)
    
    if show_lines:
      for s in section_list:
        n = s[1] - s[0]
        plt.plot(np.arange(0,n*0.0674-0.001,0.0674), trace[s[0]:s[1]], linewidth=0.5)

    else:
      n = 15
      trace_all = np.ones((len(section_list),n))*np.nan
      for j in range(len(section_list)):

        s = section_list[j]

        trace_section = trace[s[0]:s[1]]
        trace_all[j,:trace_section.shape[0]] = trace_section

      mean_trace = np.nanmean(trace_all, axis=0)
      std = np.nanstd(trace_all, axis=0)

      t = np.arange(0,n*0.0674-0.001,0.0674)
      ax.plot(t, mean_trace)
      ax.fill_between(t, mean_trace-std, mean_trace, facecolor='blue', alpha=0.3)
      ax.fill_between(t, mean_trace, mean_trace+std, facecolor='blue', alpha=0.3)
        
    plt.ylim([trace.min(), trace.max()])
    plt.title(str(np.round(angle,2)))
        
    if i == 0: 
      plt.xlabel("Time (s)")
      plt.ylabel("Response (a.u.)")
        
  plt.show()


def plot_linear_fit(xval, yval, xlab="", ylab="", return_result=True):

  Xval = sm.add_constant(xval)
  re = sm.OLS(yval, Xval).fit()

  xrng = np.max(xval) - np.min(xval)
  xlin = np.linspace(np.min(xval)-0.01*xrng,np.max(xval)+0.01*xrng,100)
  Xlin = sm.add_constant(xlin)
  dt = re.get_prediction(Xlin).summary_frame(alpha = 0.2)

  plt.figure()
  plt.plot(xval, yval, 'k.', alpha=0.7)
  plt.plot(xlin, dt["mean"], '-', alpha=0.8)
  plt.fill_between(
          xlin, dt["obs_ci_lower"], dt["obs_ci_upper"],
          alpha=0.2, color="gray")
  plt.xlabel(xlab, fontsize=12)
  plt.ylabel(ylab, fontsize=12)
  plt.show()
    
  if return_result:
    r,p = pearsonr(xval, yval)
    print("r = {}, p = {}".format(r,p))


def plot_spatial_ax(ax, xval, yval, nbins=5, xlab="", ylab="In-conn. density ($\mu m^{-1}$)"):
    
  bins = np.linspace(xval.min(), xval.max()*1.001, nbins)
  
  Xval = sm.add_constant(xval)
  re = sm.OLS(yval, Xval).fit()

  xrng = np.max(xval) - np.min(xval)
  xlin = np.linspace(np.min(xval)-0.05*xrng,np.max(xval)+0.05*xrng,10)
  Xlin = sm.add_constant(xlin)
  dt = re.get_prediction(Xlin).summary_frame(alpha = 0.2)

  ax.plot(xval, yval, 'k.', alpha=0.7)
  ax.fill_between(
      xlin, dt["obs_ci_lower"], dt["obs_ci_upper"],
      alpha=0.2, color="gray")
  ylim = ax.get_ylim()
  ax.vlines(bins, ylim[0], ylim[1], "k", linestyles="--", linewidth=0.5)
  
  ax.set_xlim([xlin.min(), xlin.max()])
  ax.set_ylim(ylim)
  
  ax.set_xlabel(xlab, fontsize=20, fontname="Helvetica")
  ax.set_ylabel(ylab, fontsize=20, fontname="Helvetica")
  ax.set_yticks(np.arange(0,0.009,0.002))