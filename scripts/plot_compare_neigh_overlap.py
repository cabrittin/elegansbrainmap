"""
compare_neigh_overlap.py

Plots distributions of Jaccard distances for overlapping ipsilateral 
neighborhoods (blue) and homologous contralateral neighborhoods (red) 
in the adult and L4.

crated: Christopher Brittin
data: 01 November 2018

"""
import os
import sys
from configparser import ConfigParser,ExtendedInterpolation
import argparse
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import mannwhitneyu
import pandas as pd
import seaborn as sns

mpl.rcParams['xtick.labelsize'] = 5
mpl.rcParams['ytick.labelsize'] = 5

SOURCE = "data/neighborhood_similarity.csv"

def get_stat_sig(df,comparison,col,group1,group2):
    d1 = df[(df['Comparison'] == comparison) & (df[col] == group1)]['Jaccard Distance'].to_numpy()
    d2 = df[(df['Comparison'] == comparison) & (df[col] == group2)]['Jaccard Distance'].to_numpy()
    tval,pval = mannwhitneyu(d1,d2)
    return pval


def run(_cfg,fout=None,source_data=None):
    cfg = ConfigParser(interpolation=ExtendedInterpolation())
    cfg.read(_cfg)

    if os.path.isfile(SOURCE):
        print(f"Found source data file: {SOURCE}")
    else:
        print(f"Source data file not found: {SOURCE}")
        print("Generating data source file now")
        print("To speed up plotting, generate source file with:")
        print("python analysis2/compare_neigh_overlap.py")
        import compare_neigh_overlap
        compare_neigh_overlap.run(_cfg,source_data=SOURCE)


    sns.set() 
    sns.set_theme(style="whitegrid")
    df = pd.read_csv(SOURCE)

    df_all = df.loc[(df['Network'] == 'all') & (df['Edge threshold'] == 'all')]
    df_edge = df.loc[(df['Network'] == 'all') & 
            (df['Edge threshold'].isin(['high','mid','low']))& 
            (df['Measure'] == 'homologous')]
    #df_lateral = df.loc[(df['Network'].isin(['ipsilateral','contralateral'])) &
    #        (df['Edge threshold'] == 'all') & (df['Measure'] == 'homologous')]


    fig,ax = plt.subplots(2,1,figsize=(2.1,3.2))
    flierprops = dict(markersize=1,marker='d',markerfacecolor='k')
    medianprops = dict(linestyle='-',linewidth=0.5,color='k')
    whiskerprops = dict(linestyle='-',linewidth=0.3,color='k')
    capprops = dict(linewidth=0.3)
    sns.boxplot(x="Comparison",y="Jaccard Distance",hue="Measure",
            palette=["#0b38c7","#ffe126"],data=df_all,width=0.3,ax=ax[0],linewidth=0.3,
            flierprops=flierprops,medianprops=medianprops,capprops=capprops)
    sns.boxplot(x="Comparison",y="Jaccard Distance",hue="Edge threshold",
            palette=["#42fc30","#A4A4A4","#7c3aff"],data=df_edge,width=0.3,linewidth=0.3,
            ax=ax[1],hue_order=["low","mid","high"],
            flierprops=flierprops,medianprops=medianprops,capprops=capprops)
    #sns.boxplot(x="Comparison",y="Jaccard Distance",hue="Network",
    #        palette=["m","g"],data=df_lateral,width=0.3,ax=ax[2])
    ax[0].set_ylim([0,1])
    ax[0].set_xlabel("")
    ax[0].set_ylabel("Jaccard Index",fontsize=7)
    ax[1].set_ylabel("Jaccard Index",fontsize=7)
    ax[1].set_xlabel("")
    ax[0].legend(loc="lower center",fontsize=6)
    ax[1].legend(loc="upper right",fontsize=6)
    for _ax in ax:
        for tick in _ax.xaxis.get_major_ticks(): tick.label.set_fontsize(7) 
        for tick in _ax.yaxis.get_major_ticks(): tick.label.set_fontsize(7)
        _ax.set_yticks([0,0.2,0.4,0.6,0.8,1.0])
        _ax.set_yticklabels(['0','0.2','0.4','0.6','0.8','1.0'],fontsize=5)

    plt.tight_layout()
    if fout: plt.savefig(fout)
    if source_data:
        sd = source_data.replace('.csv','_all.csv')
        df_all.to_csv(sd,index=False)
        sd = source_data.replace('.csv','_edge.csv')
        df_edge.to_csv(sd,index=False)

    plt.show()

if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-c','--config',
                        dest = 'config',
                        action = 'store',
                        default = CONFIG,
                        required = False,
                        help = 'Config file')
    
    parser.add_argument('-o','--fout',
                        action = 'store',
                        dest = 'fout',
                        default = None,
                        required = False,
                        help = 'Output svg file')
    
    parser.add_argument('-sd','--source_data',
                        action = 'store',
                        dest = 'source_data',
                        default = None,
                        required = False,
                        help = 'Output source data file')
    
    params = parser.parse_args()
    run(params.config,fout=params.fout,source_data = params.source_data)
