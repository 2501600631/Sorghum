from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, FuncFormatter
from matplotlib.lines import Line2D

HERE=Path(__file__).resolve().parent
TARGETS=['crude_protein','moisture','crude_ash','crude_fiber']
NAMES=['Crude protein','Moisture','Crude ash','Crude fibre']
COLORS=['#176A9C','#D18127','#657382']
METHODS=['nested_strong_whole','nested_ground_assisted','ground_reference']
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
 'font.size':9,'axes.titlesize':10,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,
 'axes.edgecolor':'#303740','axes.labelcolor':'#303740','text.color':'#202A34',
 'xtick.color':'#303740','ytick.color':'#303740','axes.linewidth':0.8,
 'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})

def setup(title,subtitle):
 fig,axs=plt.subplots(2,2,figsize=(7.2,5.5))
 fig.subplots_adjust(left=.09,right=.98,bottom=.115,top=.82,wspace=.31,hspace=.52)
 fig.text(.09,.965,title,ha='left',va='top',fontsize=13,weight='bold')
 fig.text(.09,.917,subtitle,ha='left',va='top',fontsize=8.3,color='#53606D')
 for j,ax in enumerate(axs.flat):
  ax.set_title(NAMES[j],loc='left',fontweight='bold',pad=8)
  ax.text(-.19,1.07,chr(97+j),transform=ax.transAxes,weight='bold',fontsize=12)
  ax.spines[['top','right']].set_visible(False)
  ax.tick_params(direction='out',length=3,width=.7)
  ax.set_axisbelow(True);ax.yaxis.grid(True,color='#E8EBEE',linewidth=.6)
  ax.yaxis.set_major_locator(MaxNLocator(4));ax.xaxis.set_major_locator(MaxNLocator(4))
 return fig,axs

def export(fig,number):
 for ext in ['png','pdf','svg']:
  fig.savefig(HERE/f'Figure_{number:02d}.{ext}',dpi=600)
 plt.close(fig)
 print(f'Figure {number}: workbook checks passed; PNG, PDF and SVG exported.')

wb=HERE/'Figure_03_data.xlsx'
pol=pd.read_excel(wb,sheet_name='Policy_repeats')
summary=pd.read_excel(wb,sheet_name='Plot_data')
assert len(pol)==840 and len(summary)==84
assert not pol.duplicated(['target','policy','requested_fraction','outer_repeat']).any()
for (t,p,f),g in pol.groupby(['target','policy','requested_fraction']):
 s=summary[(summary.target==t)&(summary.policy==p)&np.isclose(summary.requested_fraction,f)].iloc[0]
 assert len(g)==10
 for q,n in [(.25,'q25_rmse'),(.5,'median_rmse'),(.75,'q75_rmse')]:
  assert np.isclose(g.rmse.quantile(q),s[n],rtol=1e-10,atol=1e-12)
fig,axs=setup('Selective grinding burden and prediction error','Retrospective policies | median and IQR across 10 outer repeats | 98 coded samples')
fig.subplots_adjust(bottom=.17)
policies=['uncertainty_guided','random_matched','oracle_upper_bound']
styles=['-','--',':'];markers=['o','s','^']
for j,(t,ax) in enumerate(zip(TARGETS,axs.flat)):
 for k,p in enumerate(policies):
  g=summary[(summary.target==t)&(summary.policy==p)].sort_values('median_grinding_burden')
  x=g.median_grinding_burden.to_numpy()*100
  ax.fill_between(x,g.q25_rmse.to_numpy(),g.q75_rmse.to_numpy(),color=COLORS[k],alpha=.10,lw=0)
  ax.plot(x,g.median_rmse.to_numpy(),styles[k],marker=markers[k],markersize=3.5,lw=1.15,color=COLORS[k],markerfacecolor='white' if k==1 else COLORS[k],markeredgewidth=.8)
 ax.set_xlim(-2,102);ax.set_xticks([0,20,40,60,80,100]);ax.set_xlabel('Samples routed to grinding (%)');ax.set_ylabel('RMSE (percentage points)')
 a=summary[summary.target==t];lo=a.q25_rmse.min();hi=a.q75_rmse.max();pad=(hi-lo)*.12;ax.set_ylim(max(0,lo-pad),hi+pad)
handles=[Line2D([0],[0],color=COLORS[k],ls=styles[k],marker=markers[k],markersize=4,lw=1.2,label=s) for k,s in enumerate(['Disagreement guided','Matched random','Hindsight oracle'])]
fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.53,.04),ncol=3,frameon=False,fontsize=8.5,handlelength=2.4,columnspacing=2)
fig.text(.09,.018,'Bands are descriptive IQRs, not confidence intervals. Oracle requires known outcomes. Panel-specific y scales.',fontsize=7.8,color='#53606D')
export(fig,3)
