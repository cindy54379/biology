"""Conventional single-cell and clinical aggregate displays. No synthetic curves."""
from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from PIL import Image

import argparse, shutil
S=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description='Replot figures from supplied CSV source data; no GEO download or private intermediate inputs required.')
parser.add_argument('--output-dir',type=Path,default=S/'cytotoxic_option_figures')
args=parser.parse_args()
R=args.output_dir.resolve();O=R;R.mkdir(parents=True,exist_ok=True)
def read(name):
 path=S/name
 if not path.exists() and name.endswith('.gz'):
  path=S/name[:-3]
 return pd.read_csv(path,float_precision='round_trip')
params=json.loads((S/'Analysis_parameters.json').read_text(encoding='utf-8'))
t=read('Cells_tumor_baseline_week4.csv.gz')
late=read('Cells_exact_Fig3E.csv.gz')
annotations=read('Cluster_annotations.csv')
labels=dict(zip(annotations.cluster_id,annotations.English))
dots=read('Fig1_marker_dotplot.csv')
program_rows=read('Fig1_cluster_program_heatmap.csv')
program_order=program_rows[['program','label']].drop_duplicates()
b=dict(t=t,late=late,cut=pd.Series({'Module_1':params['module1_cutoff'],'Module_8':params['module8_cutoff']}),labels=labels,
 summary=read('Fig2_group_summary.csv'),comp=read('Fig2_quadrant_proportions.csv'),patient=read('Fig2_patient_joint_programs.csv'),
 dots=dots,genes=dots.gene.drop_duplicates().tolist(),program_rows=program_rows,programs=program_order.program.tolist(),program_names=program_order.label.tolist())
assert len(t)==36185 and len(late)==6062
assert t.cb.is_unique and late.cb.is_unique
assert 'TOX' in t and 'novershtern_eff_cd8' in t
for name in ['Sensitivity_cutoffs_and_patient_weighting.csv','Fig3_clinical_best_response.csv','Fig3_clinical_landmark_estimates.csv','Fig3_clinical_hazard_ratios.csv']:
 if (S/name).resolve()!=(O/name).resolve():shutil.copyfile(S/name,O/name)

t=b['t'];late=b['late'];cut=b['cut'];labels=b['labels'];summary=b['summary']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,'axes.titlesize':10,'xtick.labelsize':7,'ytick.labelsize':7,'axes.linewidth':.65,'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
navy='#3D5488';red='#C52E38';salmon='#EDAA91';teal='#64ADAA';grey='#707A89'
cluster_colors=['#7DBEB6','#87A0C4','#C797B8','#A4B284','#F3B982','#6C80AB','#D96D64','#C42D2B','#8B6CA8']
div=LinearSegmentedColormap.from_list('navy_white_red',[navy,'#FAFAFA',red])
exprc=LinearSegmentedColormap.from_list('expression',['#E8E9ED','#F1A58E','#B71832'])
old_effector_palette=LinearSegmentedColormap.from_list('original_effector',['#E6E7E9','#F4BCA8','#BD1F2D'])
audit=read('Effector_comparison_cluster_statistics.csv')
effector_limits=audit.loc[audit.cohort=='Original full tumor atlas',['q02','q98']].iloc[0].to_numpy(float)
groups=['Baseline','rela','nivo','combo'];names={'Baseline':'Baseline','rela':'Rela | week 4','nivo':'Nivo | week 4','combo':'Nivo + Rela | week 4'}
states=['Exhaustion-low / cytotoxic-low','Exhaustion-low / cytotoxic-high','Exhaustion-high / cytotoxic-low','Exhaustion-high / cytotoxic-high']
state_colors=['#DCE2E9','#79B4BD','#8C9DBC','#C84043']

def panel(ax,letter,title):
 ax.set_title(title,loc='left',pad=10)
 ax.text(-.06,1.05,letter,transform=ax.transAxes,fontweight='bold',fontsize=14,va='bottom')
 ax.spines[['top','right']].set_visible(False)

def save(fig,name):
 for ext in ['pdf','svg','png']:fig.savefig(O/(name+'.'+ext),dpi=300)
 fig.savefig(R/(name+'_preview.png'),dpi=110)
 plt.close(fig)

def umap(ax,d,value=None):
 d=d.sample(frac=1,random_state=941)
 if value is None:
  for i in range(1,10):
   q=d[d.cluster_id==i];ax.scatter(q.UMAP_1,q.UMAP_2,s=1.5,c=cluster_colors[i-1],lw=0,alpha=.65,rasterized=True)
  for i,q in d.groupby('cluster_id'):
   ax.text(q.UMAP_1.median(),q.UMAP_2.median(),f'C{i}',ha='center',fontsize=7,bbox=dict(facecolor='white',alpha=.8,edgecolor='none',pad=.5))
 else:
  # Display clipping is declared and never applied to exported source scores.
  if value=='novershtern_eff_cd8':
   # Restore the older display convention without changing the week-4 cohort.
   lim=effector_limits;d=d.iloc[np.random.default_rng(813).permutation(len(d))]
   im=ax.scatter(d.UMAP_1,d.UMAP_2,s=1.8,c=d[value],lw=0,cmap=old_effector_palette,vmin=lim[0],vmax=lim[1],rasterized=True)
  else:
   lim=np.quantile(d[value],[.02,.98]);d=d.sort_values(value)
   im=ax.scatter(d.UMAP_1,d.UMAP_2,s=3,c=d[value],lw=0,cmap=exprc,vmin=lim[0],vmax=lim[1],rasterized=True)
  cb=ax.figure.colorbar(im,ax=ax,orientation='horizontal',fraction=.035,pad=.14,shrink=.65)
  cb.ax.tick_params(labelsize=6);cb.set_label('Cytotoxicity module score' if value=='Module_8' else ('Deposited score' if value.startswith('Module') or value in b['programs'] else 'Deposited normalized expression'),fontsize=7)
 ax.set_xticks([]);ax.set_yticks([]);ax.set_xlabel('UMAP 1');ax.set_ylabel('UMAP 2')
 ax.spines[['top','right']].set_visible(False)

# FIGURE 1: a complete transcriptional atlas, not presence of a labeled cluster alone.
fig=plt.figure(figsize=(13.2,11.2))
gs=fig.add_gridspec(3,12,height_ratios=[2.35,2.05,2.3],left=.065,right=.965,top=.90,bottom=.095,hspace=.66,wspace=1.8)
fig.text(.065,.968,'Figure 1 | CD8 T-cell identities and exhaustion-associated transcriptional programs',fontsize=14,fontweight='bold')
fig.text(.065,.942,f'GSE225063 | Baseline + week-4 tumor atlas: {len(t):,} cells | Author embeddings, cluster IDs and scores retained',fontsize=9,color=grey)
ax=fig.add_subplot(gs[0,:4]);panel(ax,'a','Tumor CD8 atlas');umap(ax,t)
ax=fig.add_subplot(gs[0,4:8]);ax.axis('off');panel(ax,'b','Cluster identities')
for i in range(1,10):
 y=.96-(i-1)*.107;ax.scatter(.0,y,s=28,color=cluster_colors[i-1],transform=ax.transAxes,clip_on=False)
 text=labels[i]
 if i==7:text='Terminal Tex /\nTCR signaling'
 if i==8:text='Progenitor-1 Tex / terminal Tex'
 ax.text(.055,y,f'C{i}  '+text,va='center',transform=ax.transAxes,fontsize=7.2,linespacing=1.0)
ax=fig.add_subplot(gs[0,8:]);panel(ax,'c','Programs across clusters')
p=b['program_rows'].pivot(index='label',columns='cluster_id',values='program_zscore').reindex(b['program_names'])
im=ax.imshow(p,cmap=div,vmin=-2,vmax=2,aspect='auto');ax.set_yticks(range(len(p)),p.index,fontsize=6.7);ax.set_xticks(range(9),[f'C{i}' for i in range(1,10)]);ax.set_xlabel('Original cluster')
cb=fig.colorbar(im,ax=ax,fraction=.03,pad=.03);cb.ax.set_title('z-score',fontsize=6);cb.ax.tick_params(labelsize=6)
ax=fig.add_subplot(gs[1,:]);panel(ax,'d','Marker expression supports cell-state identities')
d=b['dots'];genes=b['genes']
for i in range(1,10):
 q=d[d.cluster_id==i].set_index('gene').loc[genes]
 im=ax.scatter(range(len(genes)),np.repeat(i-1,len(genes)),s=2+125*q.fraction_detected,c=q.gene_zscore,cmap=div,vmin=-2,vmax=2,lw=0)
ax.set_xticks(range(len(genes)),genes,rotation=45,ha='right',fontstyle='italic');ax.set_yticks(range(9),[f'C{i}' for i in range(1,10)]);ax.set_xlim(-.6,len(genes)-.4);ax.set_ylim(8.6,-.6)
handles=[ax.scatter([],[],s=2+125*f,c='#858E9B',label=f'{100*f:.0f}%') for f in [.25,.5,1]]
ax.legend(handles=handles,title='Cells with expression > 0',ncol=3,frameon=False,loc='upper right',bbox_to_anchor=(.99,1.32),fontsize=6,title_fontsize=6)
cb=fig.colorbar(im,ax=ax,fraction=.015,pad=.018);cb.ax.set_title('z-score',fontsize=6)
post=t[t.time_tx=='combo'];lims={c:(t[c].min(),t[c].max()) for c in ['UMAP_1','UMAP_2']}
for i,(value,title) in enumerate([('Module_1','Terminal exhaustion module'),('Module_8','Cytotoxic module'),('PDCD1','PDCD1 expression'),('HAVCR2','HAVCR2 expression')]):
 ax=fig.add_subplot(gs[2,i*3:(i+1)*3]);panel(ax,'efgh'[i],title);umap(ax,post,value)
 for c,fn in [('UMAP_1',ax.set_xlim),('UMAP_2',ax.set_ylim)]:
  lo,hi=lims[c];fn(lo-.03*(hi-lo),hi+.03*(hi-lo))
fig.text(.065,.038,f'e–h: Nivo + Rela, week-4 tumor CD8 cells (n={len(post):,}; {post.patient_id.nunique()} patients). e–h: within-panel 2nd–98th percentile color limits; high scores/expression plotted last.\nDot size: detection fraction; heatmap/dot color: program-wise/gene-wise z-score across nine clusters (display ±2).',fontsize=7,color=grey)
save(fig,'Figure_1_CD8_atlas_and_programs')

