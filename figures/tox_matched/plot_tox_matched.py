"""Generate TOX, PDCD1 and HAVCR2 panels matching Figure 1g/h geometry and style."""
from pathlib import Path
import argparse,json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.transforms import Bbox
S=Path(__file__).resolve().parent
parser=argparse.ArgumentParser();parser.add_argument('--data',type=Path,default=S/'TOX_matched_source.csv');parser.add_argument('--output-dir',type=Path,default=S/'outputs');args=parser.parse_args();O=args.output_dir;O.mkdir(parents=True,exist_ok=True)
t=pd.read_csv(args.data,float_precision='round_trip');post=t[t.time_tx=='combo'];assert len(post)==7445
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,'axes.titlesize':10,'xtick.labelsize':7,'ytick.labelsize':7,'axes.linewidth':.65,'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
cmap=LinearSegmentedColormap.from_list('expression',['#E8E9ED','#F1A58E','#B71832'])
params=[]
for gene,condition,filename in [('PDCD1','combo','PDCD1_matched'),('HAVCR2','combo','HAVCR2_matched'),('TOX','combo','TOX_matched'),('TOX','Baseline','TOX_baseline_matched')]:
 # Original Figure 1g subplot location: preserves physical axes and colorbar sizes.
 fig=plt.figure(figsize=(13.2,11.2));gs=fig.add_gridspec(3,12,height_ratios=[2.35,2.05,2.3],left=.065,right=.965,top=.90,bottom=.095,hspace=.66,wspace=1.8)
 slot=gs[2,6:9];pos=slot.get_position(fig);ax=fig.add_subplot(slot)
 d=t[t.time_tx==condition].sample(frac=1,random_state=941);lim=np.quantile(t[gene] if gene=='TOX' else d[gene],[.02,.98]);d=d.sort_values(gene)
 im=ax.scatter(d.UMAP_1,d.UMAP_2,s=3,c=d[gene],lw=0,cmap=cmap,vmin=lim[0],vmax=lim[1],rasterized=True)
 cb=fig.colorbar(im,ax=ax,orientation='horizontal',fraction=.035,pad=.14,shrink=.65);cb.ax.tick_params(labelsize=6);cb.set_label('Normalized expression',fontsize=7)
 ax.set_title(gene+' expression',loc='left',pad=10);ax.set_xticks([]);ax.set_yticks([]);ax.set_xlabel('UMAP 1');ax.set_ylabel('UMAP 2');ax.spines[['top','right']].set_visible(False)
 for c,fn in [('UMAP_1',ax.set_xlim),('UMAP_2',ax.set_ylim)]:
  lo,hi=t[c].min(),t[c].max();fn(lo-.03*(hi-lo),hi+.03*(hi-lo))
 # Fixed crop shared across all three genes; no content-dependent resizing.
 box=Bbox.from_extents(pos.x0*13.2-.38,pos.y0*11.2-.45,pos.x1*13.2+.10,pos.y1*11.2+.34)
 for ext in ['png','pdf','svg']:fig.savefig(O/(filename+'.'+ext),dpi=300,bbox_inches=box,pad_inches=0)
 fig.canvas.draw();params.append(dict(gene=gene,condition=condition,cells=len(d),scatter_area_pt2=3,alpha=1,edge_linewidth=0,vmin=float(lim[0]),vmax=float(lim[1]),axes_width_inches=ax.get_position().width*13.2,axes_height_inches=ax.get_position().height*11.2))
 plt.close(fig)
(O/'Plot_parameters.json').write_text(json.dumps(params,indent=2));print('Matched-format panels saved:',O)
