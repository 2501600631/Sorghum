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

wb=HERE/'Figure_01_data.xlsx'
raw=pd.read_excel(wb,sheet_name='OOF_predictions')
assert len(raw)==3920 and not raw.duplicated(['target','sample_id','outer_repeat']).any()
assert raw.groupby(['target','sample_id']).y_true.nunique().eq(1).all()
assert raw.groupby(['target','sample_id']).outer_repeat.nunique().eq(10).all()
g=raw.groupby(['target','sample_id'],as_index=False).agg(y_true=('y_true','first'),mean_oof_prediction=('y_pred','mean'),sd_oof_prediction=('y_pred','std'))
saved=pd.read_excel(wb,sheet_name='Plot_data')
c=g.merge(saved,on=['target','sample_id'],validate='one_to_one',suffixes=('_calc','_saved'))
for name in ['y_true','mean_oof_prediction','sd_oof_prediction']:
 assert np.allclose(c[name+'_calc'],c[name+'_saved'],rtol=1e-10,atol=1e-12)
metrics=pd.read_excel(wb,sheet_name='Repeat_metrics')
fig,axs=setup('Whole grain predictions','98 coded samples per panel | mean held-out prediction across 10 repeats')
for j,(t,ax) in enumerate(zip(TARGETS,axs.flat)):
 a=g[g.target==t];m=metrics[metrics.target==t]
 ax.errorbar(a.y_true,a.mean_oof_prediction,yerr=a.sd_oof_prediction,fmt='none',ecolor='#BACAD6',elinewidth=.6,alpha=.8,zorder=1)
 ax.scatter(a.y_true,a.mean_oof_prediction,s=15,c=COLORS[0],alpha=.78,edgecolors='white',linewidths=.35,zorder=3)
 lo=min(a.y_true.min(),(a.mean_oof_prediction-a.sd_oof_prediction).min())
 hi=max(a.y_true.max(),(a.mean_oof_prediction+a.sd_oof_prediction).max())
 pad=(hi-lo)*.10;lim=(lo-pad,hi+pad)
 ax.plot(lim,lim,ls='--',lw=.9,color='#657382',zorder=0)
 ax.set_xlim(lim);ax.set_ylim(lim)
 ax.set_xlabel('Reference content (%)');ax.set_ylabel('Predicted content (%)')
 ax.text(.045,.95,f'Median CV $R^2$ = {m.r2.median():.3f}\nMedian CV RMSE = {m.rmse.median():.3f}',transform=ax.transAxes,ha='left',va='top',fontsize=8,bbox={'facecolor':'white','edgecolor':'none','alpha':.87,'pad':2})
fig.text(.09,.025,'Vertical intervals: ±1 SD across repeats; dashed line: identity. Contents are expressed as percentages.',fontsize=8,color='#53606D')
export(fig,1)
