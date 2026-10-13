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

wb=HERE/'Figure_02_data.xlsx'
raw=pd.read_excel(wb,sheet_name='OOF_predictions')
rep=pd.read_excel(wb,sheet_name='Plot_data')
assert len(raw)==11760 and len(rep)==120
assert not raw.duplicated(['target','method','sample_id','outer_repeat']).any()
for (t,m,r),g in raw.groupby(['target','method','outer_repeat']):
 assert len(g)==98
 rmse=np.sqrt(np.mean((g.y_true-g.y_pred)**2))
 saved=rep[(rep.target==t)&(rep.method==m)&(rep.outer_repeat==r)].rmse.iloc[0]
 assert np.isclose(rmse,saved,rtol=1e-10,atol=1e-12)
fig,axs=setup('Prediction errors across measurement workflows','10 repeated nested cross-validation estimates per workflow | 98 coded samples')
for j,(t,ax) in enumerate(zip(TARGETS,axs.flat)):
 a=rep[rep.target==t].pivot(index='outer_repeat',columns='method',values='rmse')[METHODS]
 vals=[a[m].to_numpy() for m in METHODS]
 box=ax.boxplot(vals,positions=[1,2,3],widths=.46,patch_artist=True,showfliers=False,
  medianprops={'color':'#18232D','linewidth':1.35},whiskerprops={'linewidth':.8,'color':'#77828C'},capprops={'linewidth':.8,'color':'#77828C'},boxprops={'linewidth':.9})
 for b,col in zip(box['boxes'],COLORS):b.set_facecolor(col+'28');b.set_edgecolor(col)
 for k,row in enumerate(a.to_numpy()):
  xx=np.arange(1,4)+(k-4.5)*.018
  ax.plot(xx,row,color='#C3CBD1',lw=.55,alpha=.7,zorder=1)
  for q in range(3):ax.scatter(xx[q],row[q],s=19,marker=['o','s','^'][q],color=COLORS[q],edgecolors='white',linewidths=.4,zorder=4)
 ax.set_xticks([1,2,3],['Whole','Assisted','Ground']);ax.set_xlim(.5,3.5)
 ax.set_ylabel('RMSE (percentage points)')
 low=np.min(a.to_numpy());high=np.max(a.to_numpy());margin=(high-low)*.24
 ax.set_ylim(max(0,low-margin),high+margin)
fig.text(.09,.025,'Boxes: median and IQR; points: individual repeats; connectors: matched repeats. Panel-specific y scales.',fontsize=8,color='#53606D')
export(fig,2)
