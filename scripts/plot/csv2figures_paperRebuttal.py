import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from sklearn import metrics
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from organelle_measure.data import read_results

# Global Variables
plt.rcParams["figure.autolayout"]=True
plt.rcParams['font.size'] = '26'
px_x,px_y,px_z = 0.41,0.41,0.20

organelles = [
    "peroxisome",
    "vacuole",
    "ER",
    "golgi",
    "mitochondria",
    "LD"
]
subfolders  = ["paperRebuttal"]
experiments = {"rebuttal": "paperRebuttal"}
exp_names   = list(experiments.keys())
exp_folder  = [experiments[i] for i in exp_names]
extremes    = {"paperRebuttal": [0.,    100.]}
list_colors = {"rebuttal": [1,2,3,4,0,5]}

df_bycell = read_results(Path("./data/results"),subfolders,(px_x,px_y,px_z))

# DATAFRAME FOR CORRELATION COEFFICIENT
pv_bycell = df_bycell.set_index(['folder','condition','field','idx_cell'])
df_corrcoef = pd.DataFrame(index=pv_bycell.loc[pv_bycell["organelle"].eq("ER")].index)
df_corrcoef.loc[:,'cell length'] = np.sqrt(pv_bycell.loc[pv_bycell["organelle"].eq("ER"),'cell-area']/np.pi)
df_corrcoef.loc[:,'cell area'] = pv_bycell.loc[pv_bycell["organelle"].eq("ER"),'cell-area']
df_corrcoef.loc[:,'cell volume'] = pv_bycell.loc[pv_bycell["organelle"].eq("ER"),'cell-volume']
properties = []
rename_dict = {
    'mean'          : "average volume",
    'count'         : "number",
    'total'         : "total volume",
    'total-fraction': "volume fraction"
}
for orga in [*organelles,"non-organelle"]:
    for prop in ['mean','count','total','total-fraction']:
        if (orga in ["ER","vacuole","non-organelle"]) and (prop in ["count","mean"]):
            continue
        prop_new = f"{orga} {rename_dict[prop]}"
        properties.append(prop_new)
        df_corrcoef.loc[:,prop_new] = pv_bycell.loc[pv_bycell["organelle"].eq(orga),prop]
df_corrcoef.reset_index(inplace=True)

columns = ['cell length','cell area','cell volume',*properties]

# Correlation coefficient
for folder in subfolders:
    np_corrcoef = df_corrcoef.loc[df_corrcoef["folder"].eq(folder),columns].to_numpy()
    corrcoef = np.corrcoef(np_corrcoef,rowvar=False)
    fig = px.imshow(
            corrcoef,
            x=columns,y=columns,
            color_continuous_scale = "RdBu_r",range_color=[-1,1]
        )
    fig.update_layout(
        # font_family="Arial Black",
        font_size=14
    )
    fig.write_html(f'{Path("./data/REBUTTAL_8HOURS/correlation")}/corrcoef-nocond_{folder}.html')
    
    # trivial visualization change:
    # only keep lower half
    rows,cols = corrcoef.shape
    corrcoef_triangle_low = corrcoef.copy()
    for r in range(rows-1):
        for c in range(r+1,cols):
            corrcoef_triangle_low[r,c] = 0
    fig = px.imshow(
            corrcoef_triangle_low,
            x=columns,y=columns,
            color_continuous_scale = "RdBu_r",range_color=[-1,1]
        )
    fig.update_layout(
        # font_family="Arial Black",
        font_size=14
    )
    fig.write_html(f'{Path("./data/REBUTTAL_8HOURS/correlation")}/corrcoef-nocond_{folder}_lower.html')
    # only keep upper half
    corrcoef_triangle_upp = corrcoef.copy()
    for c in range(cols-1):
        for r in range(c+1,rows):
            corrcoef_triangle_upp[r,c] = 0
    fig = px.imshow(
            corrcoef_triangle_upp,
            x=columns,y=columns,
            color_continuous_scale = "RdBu_r",range_color=[-1,1]
        )
    fig.update_layout(
        # font_family="Arial Black",
        font_size=14
    )
    fig.update_xaxes(side="top")
    fig.update_yaxes(side="right")
    fig.write_html(f'{Path("./data/REBUTTAL_8HOURS/correlation")}/corrcoef-nocond_{folder}_upper.html')


    for condi in df_corrcoef.loc[df_corrcoef["folder"].eq(folder),"condition"].unique():
        np_corrcoef = df_corrcoef.loc[(df_corrcoef["folder"].eq(folder) & df_corrcoef['condition'].eq(condi)),columns].to_numpy()
        corrcoef = np.corrcoef(np_corrcoef,rowvar=False)
        fig = px.imshow(
                corrcoef,
                x=columns,y=columns,
                color_continuous_scale="RdBu_r",range_color=[-1,1]
            )
        fig.write_html(f'{Path("./data/REBUTTAL_8HOURS/correlation")}/groupby_conditions/corrcoef-nocond_{folder}_{str(condi).replace(".","-")}.html')


# Plot power law between log(cell volume) and log(cytoplasm)
dfs = df_bycell.loc[df_bycell["folder"].eq("paperRebuttal") & df_bycell["organelle"].eq("ER"),["condition","cell-volume"]]
dfs = dfs.reset_index().drop("index",axis=1)
dfs["cell-volume"] = np.log(dfs["cell-volume"])
for orga in [*organelles,"non-organelle"]:
    dfs[orga] = np.log(df_bycell.loc[df_bycell["folder"].eq("paperRebuttal") & df_bycell["organelle"].eq(orga),"total"].reset_index().drop("index",axis=1))
dfs.replace([np.inf, -np.inf], np.nan, inplace=True)
dfs.dropna(axis=0,inplace=True)
dfs = dfs[dfs['non-organelle']>2]
df_scambled = dfs[dfs["condition"].eq(100)]
df_scambled["cell-volume"] = np.random.permutation(df_scambled["cell-volume"]) 

plt.rcParams['font.size'] = '26'
fig,ax = plt.subplots(figsize=(12,9))
for i,condi in enumerate(np.sort(dfs["condition"].unique())):
    print(i,condi,list_colors["rebuttal"][i])
    # ax.axis("equal")

    # ax.set_xlim(2,7)
    # ax.set_ylim(3,7)
    ax.scatter(
        x=dfs.loc[dfs["condition"].eq(condi),"non-organelle"],
        y=dfs.loc[dfs["condition"].eq(condi),"cell-volume"],
        label=f"{condi/100*2}% glucose",
        color=sns.color_palette("tab10")[list_colors["rebuttal"][i]],alpha=0.5,edgecolors='w'
    )
    ax.set_adjustable('datalim')
inlet = fig.add_axes([0.62,0.2,0.25,0.25])
inlet.axis("equal")
# inlet.set_xbound(1.5,7.5)
# inlet.set_ybound(3,8)
inlet.scatter(
    x=df_scambled["non-organelle"],
    y=df_scambled["cell-volume"],
    label="scambled",
    color="grey",alpha=0.5,edgecolors='w'
)
ax.set_adjustable('datalim')
x_min,x_max = dfs["non-organelle"].min(),dfs["non-organelle"].max()
y_min,y_max = dfs["cell-volume"].min(),dfs["cell-volume"].max()
ax.plot(
    [x_min,x_max],[y_min,y_min+(x_max-x_min)],"k--"
)
ax.plot(
    [x_min,x_max],[y_min+0.5,y_min+0.5+(2/3)*(x_max-x_min)],"r--"
)
ax.set_xlabel(r"$log(V_{cyto})$")
ax.set_ylabel(r"$log(V_{cell})$")
# ax.legend()
plt.savefig("data/REBUTTAL_8HOURS/power_law/power_law_rectangular.png",dpi=600)
plt.close()

# generalize to other pairs of volumes.
def plot_loglog_vol(experiment,output=""):
    df_loglogs = df_bycell.loc[df_bycell["folder"].eq(experiments[experiment]) & df_bycell["organelle"].eq("ER"),["condition","cell-volume"]]
    df_loglogs = df_loglogs.reset_index().drop("index",axis=1)
    df_loglogs["cell-volume"] = np.log(df_loglogs["cell-volume"])
    for orga in [*organelles,"non-organelle"]:
        df_loglogs[orga] = np.log(df_bycell.loc[df_bycell["folder"].eq(experiments[experiment]) & df_bycell["organelle"].eq(orga),"total"].reset_index().drop("index",axis=1))
    df_loglogs.replace([np.inf, -np.inf], np.nan, inplace=True)
    df_loglogs.dropna(axis=0,inplace=True)

    fitted = np.zeros((8,8))
    plt.rcParams['font.size'] = '26'
    fig_pair,pairplots = plt.subplots(
        nrows=8,ncols=8,
        figsize=(50,40),
        sharex=False,sharey=False
    )

    uniq_conds = np.sort(df_loglogs["condition"].unique())
    for i,prop1 in enumerate(["cell-volume",*organelles,"non-organelle"]): # row
        for j,prop2 in enumerate(["cell-volume",*organelles,"non-organelle"]): # col
            dfs2fit = []
            dfs2go  = []
            pct10x = np.percentile(df_loglogs[prop2],10)
            pct10y = np.percentile(df_loglogs[prop1],10)
            print(i,j,pct10y,pct10x)
            for condi in uniq_conds:
                dfs2fit.append(
                    df_loglogs[
                        df_loglogs["condition"].eq(condi) & 
                        df_loglogs[prop2].ge(pct10x) &
                        df_loglogs[prop1].ge(pct10y)
                ])
                dfs2go.append(
                    df_loglogs[
                        df_loglogs["condition"].eq(condi) & 
                        (df_loglogs[prop2].le(pct10x) |
                        df_loglogs[prop1].le(pct10y))
                ])
            dfs2fit = pd.concat(dfs2fit,ignore_index=True)
            dfs2go  = pd.concat(dfs2go, ignore_index=True)
            fitted[i,j] = LinearRegression().fit(dfs2fit[prop2].to_numpy().reshape(-1,1),dfs2fit[prop1]).coef_[0]
            # fitted[i,j] = scipy.odr.ODR(
            #         scipy.odr.Data(
            #             dfs2fit[prop2].to_numpy(),
            #             dfs2fit[prop1].to_numpy()
            #         ),
            #         scipy.odr.unilinear
            #     ).run().beta[0]
            for k,condi in enumerate(uniq_conds):
                pairplots[i,j].scatter(
                    dfs2fit.loc[dfs2fit["condition"].eq(condi),prop2],
                    dfs2fit.loc[dfs2fit["condition"].eq(condi),prop1],
                    color=sns.color_palette("tab10")[list_colors[experiment][k]],
                    label=f"{experiment} {condi}",
                    alpha=0.5,edgecolors='w'
                )
                pairplots[i,j].scatter(
                    dfs2go.loc[dfs2go["condition"].eq(condi),prop2],
                    dfs2go.loc[dfs2go["condition"].eq(condi),prop1],
                    color='k',alpha=0.3*(k+1)/len(uniq_conds),edgecolors='w'
                )
            xmin,xmax,xmean,ymean,k = dfs2go[prop2].min(),dfs2fit[prop2].max(),dfs2fit[prop2].mean(),dfs2fit[prop1].mean(),fitted[i,j]
            pairplots[i,j].plot(
                [xmin,xmax],[k*(xmin-xmean)+ymean,k*(xmax-xmean)+ymean],
                'k--'
            )
            # if i==7:
            #     pairplots[i,j].set_xlabel(f"log[V({prop2})]")
            # if j==0:
            #     pairplots[i,j].set_ylabel(f"log[V({prop1})]")
    pairplots[7,7].legend()
    if not Path(f"{output}/power_law/").exists():
        Path(f"{output}/power_law/").mkdir()
    fig_pair.savefig(f"{output}/power_law/pairwise_{experiment}.png",dpi=600)
    np.savetxt(f"{output}/power_law/pairwise_{experiment}.txt",fitted)
    print(f"Pairwise log-log regression: {experiment}")
    return None

for exps in exp_names:
    plot_loglog_vol(exps,output="data/REBUTTAL_8HOURS")


# PCA 
def make_pca_plots(experiment,property,groups=None,has_volume=False,is_normalized=False,non_organelle=False,saveto="./data/"):
    for pca_subfolder in [ "pca_data/",
                           "pca_compare/",
                           "pca_projection_extremes/",
                           "pca_projection_all_plt/",
                           "pca_projection_all/"
                         ]:
        if not (Path(saveto)/pca_subfolder).exists():
            (Path(saveto)/pca_subfolder).mkdir()
        continue
    folder = experiments[experiment]
    name = f"{'all-conditions' if groups is None else 'extremes'}_{'has-cytoplasm' if non_organelle else 'no-cytoplasm'}_{'cell-volume' if has_volume else 'organelle-only'}_{property}_{'norm-mean-std' if is_normalized else 'raw'}"
    print("PCA Anaysis: ",folder,name)

    df_orga_perfolder = df_bycell[df_bycell["folder"].eq(folder)]
    df_orga_perfolder.set_index(["condition","field","idx_cell"],inplace=True)
    idx = df_orga_perfolder.groupby(["condition","field","idx_cell"]).count().index
    
    columns = [*organelles,"non-organelle"] if non_organelle else organelles
    df_pca = pd.DataFrame(index=idx,columns=columns)
    num_pc = 7 if non_organelle else 6
    
    for orga in columns:
        df_pca[orga] = df_orga_perfolder.loc[df_orga_perfolder["organelle"].eq(orga),property]
    
    if has_volume:
        df_pca["cell-volume"] = df_orga_perfolder.loc[df_orga_perfolder["organelle"].eq("ER"),"cell volume"]
        columns = ["cell volume",*columns]
        num_pc += 1
    
    if is_normalized:
        for col in columns:
            df_pca[col] = (df_pca[col]-df_pca[col].mean())/df_pca[col].std()
    
    df_pca.reset_index(inplace=True)

    # Find the the direction of the condition change:
    df_centroid = df_pca.groupby("condition")[columns].mean()
    fitter_centroid = LinearRegression(fit_intercept=False)
    np_centroid = df_centroid.to_numpy()
    fitter_centroid.fit(np_centroid,np.ones(np_centroid.shape[0]))
    
    vec_centroid_start = df_centroid.loc[groups[0],:].to_numpy()
    vec_centroid_start[-1] = (1 - np.dot(fitter_centroid.coef_[:-1],vec_centroid_start[:-1]))/fitter_centroid.coef_[-1]

    vec_centroid_ended = df_centroid.loc[groups[-1],:].to_numpy()
    vec_centroid_ended[-1] = (1 - np.dot(fitter_centroid.coef_[:-1],vec_centroid_ended[:-1]))/fitter_centroid.coef_[-1]

    vec_centroid = vec_centroid_ended - vec_centroid_start
    vec_centroid = vec_centroid/np.linalg.norm(vec_centroid)
    np.savetxt(str(Path(saveto)/"pca_data/condition-vector_{folder}_{name}.txt"),vec_centroid)

    # Get Principal Components (PCs)
    np_pca = df_pca[columns].to_numpy()
    pca = PCA(n_components=num_pc)
    pca.fit(np_pca)
    pca_components = pca.components_
    pca_var_ratios = pca.explained_variance_ratio_

    # Calculate cos<condition,PCs>, and realign the PCs
    cosine_pca = np.dot(pca_components,vec_centroid)
    for c in range(len(cosine_pca)):
        if cosine_pca[c] < 0:
            pca_components[c] = -pca_components[c]
    cosine_pca = np.abs(cosine_pca)
    # save and plot PCs without sorting.
    np.savetxt(f'{saveto}/pca_data/cosine_{folder}_{name}.txt',cosine_pca)
    np.savetxt(f'{saveto}/pca_data/pca-components_{folder}_{name}.txt',pca_components)
    np.savetxt(f'{saveto}/pca_data/pca-explained-ratio_{folder}_{name}.txt',pca_var_ratios)
    fig_components = px.imshow(
        pca_components,
        x=columns, y=[f"PC{i}" for i in range(num_pc)],
        color_continuous_scale="RdBu_r", color_continuous_midpoint=0
    )
    fig_components.write_html(f'{saveto}/pca_data/pca_components_{folder}_{name}.html')

    # Sort PCs according to the cosine
    arg_cosine = np.argsort(cosine_pca)[::-1]
    pca_components_sorted = pca_components[arg_cosine]
    # save and plot the PCs with sorting
    np.savetxt(f'{saveto}/pca_compare/condition-sorted-index_{folder}_{name}.txt',arg_cosine)
    np.savetxt(f'{saveto}/pca_compare/condition-sorted-cosine_{folder}_{name}.txt',cosine_pca[arg_cosine])
    np.savetxt(f'{saveto}/pca_compare/condition-sorted-pca-components_{folder}_{name}.txt',pca_components_sorted)
    fig_components_sorted = px.imshow(
        pca_components_sorted,
        x=columns, y=[f"PC{i}" for i in arg_cosine],
        color_continuous_scale="RdBu_r", color_continuous_midpoint=0
    )
    fig_components_sorted.write_html(f'{saveto}/pca_compare/condition-sorted-pca_components_{folder}_{name}.html')
    # plot the cosine 
    plt.figure(figsize=(15,12))
    plt.barh(np.arange(num_pc),cosine_pca[arg_cosine[::-1]],align='center')
    plt.yticks(np.arange(num_pc),[f"PC{i}" for i in arg_cosine[::-1]])
    plt.xlabel(r"$cos\left<condition\ vector,PC\right>$")
    plt.title(f"{folder}")
    plt.savefig(f'{saveto}/pca_compare/condition-sorted-cosine_{folder}_{name}.png')
    plt.close()


    # Draw projections onto the PCs
    for i_pc in range(len(pca_components)):
        base = pca_components[i_pc]
        df_pca[f"proj{i_pc}"] = df_pca.apply(lambda x:np.dot(base,x.loc[columns]),axis=1)
    pc2proj = arg_cosine[:3]
    df_pca_extremes = df_pca.loc[df_pca["condition"].isin(groups)]

    # 3d projection
    figproj = plt.figure(figsize=(15,12))
    ax = figproj.add_subplot(projection="3d")
    for condi in groups[::-1]:
        pc_x = df_pca_extremes.loc[df_pca_extremes["condition"].eq(condi),f"proj{pc2proj[0]}"],
        pc_y = df_pca_extremes.loc[df_pca_extremes["condition"].eq(condi),f"proj{pc2proj[1]}"],
        pc_z = df_pca_extremes.loc[df_pca_extremes["condition"].eq(condi),f"proj{pc2proj[2]}"],
        ax.scatter(
            pc_x, pc_y, pc_z,
            s=55,alpha=0.3,label=f"{condi}"
        )
    ax.set_xlabel(f"proj {pc2proj[0]}")
    ax.set_ylabel(f"proj {pc2proj[1]}")
    ax.set_zlabel(f"proj {pc2proj[2]}")
    ax.xaxis.pane.set_edgecolor('black')
    ax.yaxis.pane.set_edgecolor('black')
    ax.zaxis.pane.set_edgecolor('black')
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.set_xlim(*(np.percentile(df_pca_extremes[f"proj{pc2proj[0]}"].to_numpy(),[1,99])+np.array([-0.1,0.1])))
    ax.set_ylim(*(np.percentile(df_pca_extremes[f"proj{pc2proj[1]}"].to_numpy(),[1,99])+np.array([-0.1,0.1])))
    ax.set_zlim(*(np.percentile(df_pca_extremes[f"proj{pc2proj[2]}"].to_numpy(),[1,99])+np.array([-0.1,0.1])))
    ax.legend(loc=(1.04,1.0))
    figproj.savefig(f'{saveto}/pca_projection_extremes/pca_projection3d_{folder}_{name}_pc{"".join([str(p) for p in pc2proj])}.png')
    plt.close(figproj)

    # 3d projections, all conditions
    for d,condi in enumerate(np.sort(df_pca["condition"].unique())):
        if condi == groups[1]:
            continue
        figproj = plt.figure(figsize=(15,12))
        ax = figproj.add_subplot(projection="3d")
        
        pc_x = df_pca.loc[df_pca["condition"].eq(groups[-1]),f"proj{pc2proj[0]}"],
        pc_y = df_pca.loc[df_pca["condition"].eq(groups[-1]),f"proj{pc2proj[1]}"],
        pc_z = df_pca.loc[df_pca["condition"].eq(groups[-1]),f"proj{pc2proj[2]}"],
        ax.scatter(
            pc_x, pc_y, pc_z,
            edgecolor='white',facecolor=sns.color_palette('tab10')[0],
            s=55,alpha=0.3,label=f"{groups[-1]}"
        )
        
        pc_x = df_pca.loc[df_pca["condition"].eq(condi),f"proj{pc2proj[0]}"],
        pc_y = df_pca.loc[df_pca["condition"].eq(condi),f"proj{pc2proj[1]}"],
        pc_z = df_pca.loc[df_pca["condition"].eq(condi),f"proj{pc2proj[2]}"],
        ax.scatter(
            pc_x, pc_y, pc_z,
            edgecolor='white',facecolor=sns.color_palette('tab10')[list_colors[experiment][d]],
            s=55,alpha=0.3,label=f"{condi}"
        )
        ax.set_xlabel(f"proj {pc2proj[0]}")
        ax.set_ylabel(f"proj {pc2proj[1]}")
        ax.set_zlabel(f"proj {pc2proj[2]}")
        ax.xaxis.pane.set_edgecolor('black')
        ax.yaxis.pane.set_edgecolor('black')
        ax.zaxis.pane.set_edgecolor('black')
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False
        ax.set_xlim(*(np.percentile(df_pca_extremes[f"proj{pc2proj[0]}"].to_numpy(),[1,99])+np.array([-0.1,0.1])))
        ax.set_ylim(*(np.percentile(df_pca_extremes[f"proj{pc2proj[1]}"].to_numpy(),[1,99])+np.array([-0.1,0.1])))
        ax.set_zlim(*(np.percentile(df_pca_extremes[f"proj{pc2proj[2]}"].to_numpy(),[1,99])+np.array([-0.1,0.1])))
        ax.legend(loc=(1.04,0.5))
        figproj.savefig(f'{saveto}/pca_projection_all_plt/pca_projection3d_{folder}_{name}_condi-{str(condi).replace(".","-")}_pc{"".join([str(p) for p in pc2proj])}.png')
        plt.close(figproj)
 
    # 2d projections
    sns.set_style("whitegrid")
    for first,second in ((0,1),(0,2),(1,2)):
        plt.figure(figsize=(15,12))
        sns_plot = sns.scatterplot(
            data=df_pca_extremes[df_pca_extremes["condition"].eq(groups[0])],
            x=f"proj{pc2proj[first]}",y=f"proj{pc2proj[second]}",
            color=sns.color_palette("tab10")[1],s=100,alpha=0.75
        )
        sns_plot = sns.scatterplot(
            data=df_pca_extremes[df_pca_extremes["condition"].eq(groups[1])],
            x=f"proj{pc2proj[first]}",y=f"proj{pc2proj[second]}",
            color=sns.color_palette("tab10")[0],s=100,alpha=0.75
        )
        sns_plot.figure.savefig(f'{saveto}/pca_projection_extremes/pca_projection2d_{folder}_{name}_pc{pc2proj[first]}{pc2proj[second]}.png')
        plt.close()

    # 2d projections, all conditions
    for d,condi in enumerate(np.sort(df_pca["condition"].unique())):
        if condi == groups[1]:
            continue
        for first,second in ((0,1),(0,2),(1,2)):
            plt.figure(figsize=(15,12))
            sns_plot = sns.scatterplot(
                data=df_pca[df_pca["condition"].eq(groups[1])],
                x=f"proj{pc2proj[first]}",y=f"proj{pc2proj[second]}",
                color=sns.color_palette("tab10")[0],s=100,alpha=0.75
            )
            sns_plot = sns.scatterplot(
                data=df_pca[df_pca["condition"].eq(condi)],
                x=f"proj{pc2proj[first]}",y=f"proj{pc2proj[second]}",
                color=sns.color_palette("tab10")[list_colors[experiment][d]],s=100,alpha=0.75
            )
            plt.xlim(*(np.percentile(df_pca_extremes[f"proj{pc2proj[first]}"].to_numpy(),[1,99])+np.array([-0.1,0.1])))
            plt.ylim(*(np.percentile(df_pca_extremes[f"proj{pc2proj[second]}"].to_numpy(),[1,99])+np.array([-0.1,0.1])))
            sns_plot.figure.savefig(f'{saveto}/pca_projection_all_plt/pca_projection2d_{folder}_{name}_condi-{str(condi).replace(".","-")}_pc{pc2proj[first]}{pc2proj[second]}.png')
            plt.close()

    # draw interactive HTML with Plotly:
    figproj = go.Figure()
    for condi in pd.unique(df_pca["condition"]):
        figproj.add_trace(
            go.Scatter3d(
                x=df_pca.loc[df_pca["condition"].eq(condi),f"proj{pc2proj[0]}"],
                y=df_pca.loc[df_pca["condition"].eq(condi),f"proj{pc2proj[1]}"],
                z=df_pca.loc[df_pca["condition"].eq(condi),f"proj{pc2proj[2]}"],
                name=condi,
                mode="markers",
                marker=dict(
                            size=4,
                            # color=figcolors[j],
                            opacity=0.8
                       )   
            )
        )
    figproj.update_layout(
        scene=dict(
            xaxis=dict(
                title=f"proj{pc2proj[0]}",
                backgroundcolor='rgba(0,0,0,0)',
                gridcolor='grey',
                showline = True,
                zeroline=True,
                zerolinecolor='black'
            ),
            yaxis=dict(
                title=f"proj{pc2proj[1]}",
                backgroundcolor='rgba(0,0,0,0)',
                gridcolor='grey',
                showline=True,
                zeroline=True,
                zerolinecolor='black'
            ),
            zaxis=dict(
                title=f"proj{pc2proj[2]}",
                backgroundcolor='rgba(0,0,0,0)',
                gridcolor='grey',
                showline = True,
                zeroline=True,
                zerolinecolor='black'
            )
        )
    )
    figproj.write_html(f'{Path(f"{saveto}/pca_projection_all")}/pca_projection3d_{folder}_{name}_pc{"".join([str(p) for p in pc2proj])}.html')
    return None

make_pca_plots(
    "rebuttal","total-fraction",
    groups=extremes[experiments["rebuttal"]],
    has_volume=False,is_normalized=True,non_organelle=False,
    saveto="data/REBUTTAL_8HOURS"
)
