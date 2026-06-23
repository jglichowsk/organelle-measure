# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import os
import sys
sys.path.append(r'C:\Users\jglic\Documents\School\WashU\Mukherji Lab\organelle-measure\organelle_measure')
from pathing_variables import expmt_path
from tools import batch_apply
# %%
## Separate cells within dataframe into two groups: those that start in G1 and those in S/G2/M. For each, the next major cell cycle checkpoint
    # serves as the reference point. For example, cells that start in G1 will have the G1/S checkpoint as their reference point for binning
    # organelle measurements in time, and as an alignment point if/when applicable. 

## So split into two groups, then within each determine frames between which cell cycle checkpoint occurs, use this to determine waiting times,
    # and then use those for plots & analysis.

## Find and replace g1_ and s_ with G1_ and S_ respectively after done typing, if still want to. 

#fractional/percent error in organelle volume fraction measurements for ey rbow data
#############These error bars are currently meaningless; update ###############
org_err={
"er":0.01,
"px":0.01,
"vo":0.01,
"mt":0.01,
"gl":0.01,
"ld":0.01
}
# orgs=['er','px','vo','mt','gl','ld']
orgs=['ld','gl','vo']
xbound=30
# %% 3c-ey2795 pathing
acdc_paths=os.listdir(expmt_path+'/cell_measure')
# %% pathing for debugging convenience
# pacdc=r'C:\Users\jglic\Downloads\11-5-25 myo1-mLemon\Run_2\BF-timelapse_acdc_output_cpsam-d10-ms8.csv'
# pmyo1=r'C:\Users\jglic\Downloads\11-5-25 myo1-mLemon\Run_2\eyrbow-yellowconfocal_zstack-avg-proj_myo1_fov2_upscaled.csv'
# pbend=r'C:\Users\jglic\Downloads\11-5-25 myo1-mLemon\Run_2\BF-timelapse_manual-annotations.csv'
# palign=r'C:\Users\jglic\Downloads\11-5-25 myo1-mLemon\Run_2\manual_alignment_key_upscaled.csv'
# pmask=r'C:\Users\jglic\Downloads\11-5-25 myo1-mLemon\Run_2\BF-timelapse_cpsam-d10-ms8.tif'

# %% Parse file name metadata
def parse_meta_organelle(name: str):
    """
    Args: Name is the stem of the ORGANELLE measure csv file.
    Outptuts: Dictionary containing experiment metadata
    """
    #Unpack experiment metadata according to file naming convention. Field is fov #, and time is pre- or post-BF capture.
    organelle,date,strain,condition,field_time=name.split('_')
    field,time=field_time.split('-')
    return {
        "organelle":  organelle,
        "date":       date,
        "strain":     strain,
        "condition":  condition,
        "field":      field,
        "time":       time
    } #########don't actually need time???

# %% Linear regression
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
def linreg(dataset: np.ndarray, n: int=1000, train_frac=0.5, plot_graph=False):
    """
    Args: Array of data to perform regression upon, number of runs, fraction of data to use for training, optional graph.
    Outputs: List of regression coeff scores.
    """
    scores =[]
    for run in range(n):
        train,test=train_test_split(dataset,test_size=0.5)
        reg=LinearRegression(fit_intercept=False).fit(train[:,0].reshape(-1,1),train[:,1].reshape(-1,1)) #reshaping to revert back to vertical array
        score=reg.score(test[:,0].reshape(-1,1),test[:,1].reshape(-1,1))
        scores.append(score)
    return scores

# fig,axes=plt.subplots(nrows=6,ncols=2,figsize=(15,15))
# for i in range(len(orgs)):
#     org=orgs[i]
#     n0,b0,p0=axes[i,0].hist(sdict[org][0],bins='doane', alpha=0.75, edgecolor='black', color='b', label=org+' G1 cells')
#     axes[i,0].legend(fontsize='x-large')
#     n1,b1,p1=axes[i,1].hist(sdict[org][1], bins='doane', alpha=0.75, edgecolor='black', color='y', label=org+' S cells')
#     axes[i,1].legend(fontsize='x-large')    
#     # axes[i,0].set_title('G1 '+org+' (n='f"{len(sdict[org][0])})")
#     # axes[i,0].set_xlabel('Coeff of Determination')
#     # axes[i,1].set_title('S '+org+' (n='f"{len(sdict[org][1])})")
#     # axes[i,1].set_xlabel('Coeff of Determination')
# fig.supxlabel('Coeff of Determination')
# fig.tight_layout()

# %% Preprocess and sort input data by cc
def pp_acdc_output(df: pd.DataFrame) -> pd.DataFrame:
    """
    Args: Dataframe containing Cell-ACDC output metrics and analysis information.
    Outputs: Dataframe containing the same information sans dead or manually-excluded cells.
    """
    # dead_indices=df[df['is_cell_dead']==1].index
    if len(df['is_cell_excluded'])>0:
        excl_indices=df[df['is_cell_excluded']==1].index
        pp_df=df.drop(excl_indices)
        pp_df.columns=pp_df.columns.str.strip() #remove leading and trailing spaces
    return pp_df

def cc_sort(cell_path: str) -> (pd.DataFrame, pd.DataFrame):
    """
    Args: Path to csv containing cell cycle, cell size, and organelle metric information. 
    Outputs: Two cell metric dataframes (mothers AND daughters) distinguished by initial cell cycle position.
    """
    raw_df=pd.read_csv(cell_path)
    cell_df=pp_acdc_output(raw_df)
    g1_cells=[]
    s_cells=[]
    nb_mothers=0
    nd_mothers=0
    for cell_id in np.unique(cell_df['Cell_ID'].values):
        cell_rows=cell_df[cell_df['Cell_ID']==cell_id]
        if cell_rows['cell_cycle_stage'].values[0]=='G1':
            if len(cell_rows.loc[cell_rows.cell_cycle_stage=='S'])>0:
                g1_cells.append(cell_rows)
            # else: #exclude mothers that will not bud
            #     nb_mothers+=1
        elif cell_rows['relationship'].values[0]=='mother': 
            if len(cell_rows.loc[cell_rows.cell_cycle_stage=='G1'])>0:
                s_cells.append(cell_rows)
            # else: #exclude first-frame mothers that will not divide
            #     nd_mothers+=1
        elif cell_rows['relationship'].values[0]=='bud':
            if cell_rows['frame_i'].values[0]==0 and len(cell_rows.loc[cell_rows.cell_cycle_stage=='G1', 'frame_i'])>0:
                s_cells.append(cell_rows)
            elif cell_rows['frame_i'].values[0]>0:
                g1_cells.append(cell_rows)
            else:
                continue
        
        else:
            print('Error: Cell state not recognized for cell ' + str(cell_id))

    g1_output=pd.concat(g1_cells, ignore_index=True)
    s_output=pd.concat(s_cells, ignore_index=True)

    return g1_output, s_output
    
# for path in acdc_paths:
#     gg,ss=cc_sort(expmt_path+'/cell_measure'+'/'+path)
#     print(len(np.unique(gg.loc[gg.relationship=='mother','Cell_ID'].values)))
#     print(len(np.unique(ss.loc[ss.relationship=='mother','Cell_ID'].values)))
# %% Calc cc phase lengths
def extract_streaks(plist: list):
    """
    Args: List of cc phase annotations for a given cell.
    Outputs: Three lists together containing complete cc phase lengths, corresponding cc phase, and # of complete cc phases.
    """
    ind=1
    lens=[]
    phases=[]
    for i in range(len(plist)): #extract lengths of continuous streaks
        if 0<i<(len(plist)-1):
            if plist[i]==plist[i-1]:
                ind+=1
            else:
                lens.append(ind)
                phases.append(plist[i-1])
                ind=1
        if i==(len(plist)-1):
            lens.append(ind)
            phases.append(plist[i])
    if len(lens)<=2:
        # return 0,phases,0 #####################
        return [],[],[]
    else:
        del lens[0] #drop any streaks contacting the edges 
        del lens[-1]
        del phases[0]
        del phases[-1]
        return lens, phases, len(lens)

def cc_lengths(cell_path: str):
    """
    Args: Path to cell segm output csv
    Outputs: Arrays of frame lengths of G1, S/G2/M phases and lists of complete cc-phase # distribution and composition.
    """
    g1_df,s_df=cc_sort(cell_path)
    g1_lens, s_lens, lens_dist =[],[],[]
    ind=0
    for cell_id in np.unique(g1_df.loc[g1_df.relationship=='mother','Cell_ID'].values): #for each G1 mother
        # ind+=1
        cell_rows=g1_df.loc[g1_df.Cell_ID==cell_id].copy()
        plist=cell_rows['cell_cycle_stage'].values
        lens,phases,lens_len = extract_streaks(plist)
        # if len(lens)==0: 
        #     print('gabbagool')
        if len(phases)>0:
            lens_dist.append(lens_len)
            for i in range(len(lens)):
                if phases[i]=='G1':
                    g1_lens.append(lens[i])
                else:
                    s_lens.append(lens[i])
        # if cell_rows.loc[cell_rows.cell_cycle_stage=='S', 'frame_i'].values[-1]<cell_rows['frame_i'].values[-1]: #if G2/M transition happens
        #     s_lens.append(len(cell_rows.loc[cell_rows.cell_cycle_stage=='S'].values))
    
    for cell_id in np.unique(s_df.loc[s_df.relationship=='mother','Cell_ID'].values): #for each S mother
        # ind+=1
        cell_rows=s_df.loc[s_df.Cell_ID==cell_id].copy()
        plist=cell_rows['cell_cycle_stage'].values
        lens,phases,lens_len = extract_streaks(plist)
        if len(phases)>0:
            lens_dist.append(lens_len)
            for i in range(len(lens)):
                if phases[i]=='G1':
                    g1_lens.append(lens[i])
                else:
                    s_lens.append(lens[i])
        # if cell_rows.loc[cell_rows.cell_cycle_stage=='G1', 'frame_i'].values[-1]<cell_rows['frame_i'].values[-1]: #if G1/S transition happens
        #     g1_lens.append(len(cell_rows.loc[cell_rows.cell_cycle_stage=='G1'].values))
    # print(ind)
    return g1_lens, s_lens, lens_dist

# g_lens,s_lens,lens_l=[],[],[]
# for path in acdc_paths:
#     gl,sl,lens_dist=cc_lengths(expmt_path+'/cell_measure'+'/'+path)
#     g_lens.append(gl)
#     s_lens.append(sl)
#     lens_l.append(lens_dist)
# g_lens=[element for innerList in g_lens for element in innerList]
# s_lens=[element for innerList in s_lens for element in innerList]
# lens_l=[element for innerList in lens_l for element in innerList]
# plt.figure()
# plt.xlabel('# frames (5 min interval)')
# plt.title('CC Phase Annotation Lengths')
# n0,b0,p0=plt.hist(s_lens,bins='doane', alpha=0.75, edgecolor='black', color='b', label='S/G2/M Phase n=('+f'{len(s_lens)})')
# n1,b1,p1=plt.hist(g_lens, bins=b0, alpha=0.75, edgecolor='black', color='y', label='G1 Phase n=('+f'{len(g_lens)})')
# n2,b2,p2=plt.hist(lens_l, bins=np.arange(-0.5,5.5,1), alpha=0.75, edgecolor='black', color='r', label='# Complete cc phases n=('+f'{len(lens_l)})')
# plt.legend()
# %% find cc transitions
def find_cc_transitions(g1_df: pd.DataFrame, s_df: pd.DataFrame):
    """
    Args: cc_sort output Dataframes containing cells in G1 and S/G2/M phases respectively at frame 0.
    Outputs: Same dataframes with new index column relative to cc checkpoint progression.
                Currently dropping g1 but not s buds###############
    """
    # ind=0
    # for cell_id in np.unique(g1_df.loc[g1_df.relationship=='bud','Cell_ID'].values):
    #     ind+=1
    # print(ind)
    g1_dfs=[]
    for cell_id in np.unique(g1_df.loc[g1_df.relationship=='mother','Cell_ID'].values): #for each G1 mother
        cell_rows=g1_df.loc[g1_df.Cell_ID==cell_id].copy()
        start_frame=cell_rows.loc[cell_rows.cell_cycle_stage=='S', 'frame_i'].values[0]
        start_index=[(cell_rows['frame_i'].values[i]-start_frame) for i in range(len(cell_rows))]
        cell_rows.loc[cell_rows.Cell_ID==cell_id, 'Start_Index']=start_index
        plist=cell_rows['cell_cycle_stage'].values
        lens,phases,lens_len = extract_streaks(plist)

        if len(lens)==0: #those cells with ONE annotated cc transition
            cell_rows.loc[cell_rows.Cell_ID==cell_id, 'Transitions']= 1
            g1_dfs.append(cell_rows)

        elif len(lens)==1: #those cells with TWO annotated cc transitions
            div_frame=cell_rows.loc[cell_rows.cell_cycle_stage=='S', 'frame_i'].values[-1]+1
            div_index=[(cell_rows['frame_i'].values[i]-div_frame) for i in range(len(cell_rows))]
            cell_rows.loc[cell_rows.Cell_ID==cell_id, 'Div_Index']=div_index
            cell_rows.loc[cell_rows.Cell_ID==cell_id, 'Transitions']= 2
            g1_dfs.append(cell_rows)

        elif len(lens)==2: #those cells with THREE annotated cc transitions
            # print('deux')
            continue
        else: #those cells with >THREE annotated cc transitions
            # print('trois')
            continue
    g1_mothers=pd.concat(g1_dfs)

    smom_dfs=[]
    sbud_dfs=[] ################### propagate mult cc transitions code
    for cell_id in np.unique(s_df['Cell_ID'].values): #for each S/G2/M cell
        cell_rows=s_df.loc[s_df.Cell_ID==cell_id].copy()
        div_frame=cell_rows.loc[cell_rows.cell_cycle_stage=='G1', 'frame_i'].values[0]
        div_index=[(cell_rows['frame_i'].values[i]-div_frame) for i in range(len(cell_rows))]
        cell_rows.loc[cell_rows.Cell_ID==cell_id, 'Div_Index']=div_index
        plist=cell_rows['cell_cycle_stage'].values
        lens,phases,lens_len = extract_streaks(plist)

        if len(lens)==0: #those cells with ONE annotated cc transition
            cell_rows.loc[cell_rows.Cell_ID==cell_id, 'Transitions']= 1
            if cell_rows['relationship'].values[0]=='mother':
                smom_dfs.append(cell_rows)
            else:
                sbud_dfs.append(cell_rows)

        elif len(lens)==1: #those cells with TWO annotated cc transitions
            start_frame=cell_rows.loc[cell_rows.cell_cycle_stage=='G1', 'frame_i'].values[-1]+1
            start_index=[(cell_rows['frame_i'].values[i]-start_frame) for i in range(len(cell_rows))]
            cell_rows.loc[cell_rows.Cell_ID==cell_id, 'Start_Index']=start_index
            cell_rows.loc[cell_rows.Cell_ID==cell_id, 'Transitions']= 2
            if cell_rows['relationship'].values[0]=='mother':
                smom_dfs.append(cell_rows)
            else:
                sbud_dfs.append(cell_rows)

        elif len(lens)==2: #those cells with THREE annotated cc transitions
            # print('deux')
            continue
        else: #those cells with >THREE annotated cc transitions
            # print('trois')
            continue

    s_mothers=pd.concat(smom_dfs)
    s_buds=pd.concat(sbud_dfs)

    return g1_mothers, s_mothers, s_buds

# %% cc size analysis
# def plot_cc_scatter(df: pd.DataFrame, axis: int, org_label: str):
#     fig,axes=plt.subplots(nrows=1,ncols=2)

#     axes[axis].scatter(offset,vol_frac,c='m',marker='o') #scatter of organelle volume fraction versus relative frame index
#     axes[axis].errorbar(offset, vol_frac, yerr=org_err[org_label]*vol_frac, c='m',marker='o')

def cc_size(cell_dfpath: str, norm_size: bool=False, xbound: int=30):
    """
    Args: Path to cell metric csv file.
    Outputs: Graph displaying cell size profiles vs cell cycle position.
    """
    g1_df, s_df=cc_sort(cell_dfpath)
    g1_moms,s_moms,s_buds=find_cc_transitions(g1_df,s_df)
    g1_buds=g1_df.loc[g1_df['relationship']=='bud']

    fig,axes=plt.subplots(nrows=1,ncols=2)
    for frame in np.unique(g1_moms['Relative_Index'].values):
        axes[0].errorbar(frame, np.mean(g1_moms.loc[g1_moms.Relative_Index==frame,'cell_area_pxl'].values), yerr= np.std(g1_moms.loc[g1_moms.Relative_Index==frame,'cell_area_pxl'].values), c='m',marker='o')
    
    # for cell_id in np.unique(g1_df['Cell_ID'].values):
    #     if g1_df.loc[g1_df.Cell_ID==cell_id, 'relationship'].values[0]=='mother':
    #         axes[0].plot(g1_moms.loc[g1_moms.Cell_ID==cell_id, 'Relative_Index'].values, g1_moms.loc[g1_moms.Cell_ID==cell_id,'cell_area_pxl'].values,'k-',label='mother')
    #     else:
    #         axes[0].plot(g1_buds.loc[g1_buds.Cell_ID==cell_id, 'frame_i'].values, g1_buds.loc[g1_buds.Cell_ID==cell_id,'cell_area_pxl'].values,'b--',label='bud') 
    for cell_id in np.unique(s_df['Cell_ID'].values):
        if s_df.loc[s_df.Cell_ID==cell_id, 'relationship'].values[0]=='mother':
            axes[1].plot(s_moms.loc[s_moms.Cell_ID==cell_id, 'Div_Index'].values, s_moms.loc[s_moms.Cell_ID==cell_id,'cell_area_pxl'].values, 'k-',label='mother') 
        # else:
        #     axes[1].plot(s_buds.loc[s_buds.Cell_ID==cell_id, 'frame_i'].values, s_buds.loc[s_buds.Cell_ID==cell_id,'cell_area_pxl'].values,'b--',label='bud') 
    

    axes[0].set_title('Cell size vs CC position')
    axes[0].set_xlabel('Frames relative to Start (5 min interval)')
    axes[0].set_ylabel('cell area (px)')
    # axes[0].legend()
    axes[0].set_xlim(-xbound,xbound)
    axes[1].set_title('Cell size vs CC position')
    axes[1].set_xlabel('Frames relative to Division (5 min interval)')
    axes[1].set_ylabel('cell area (px)')
    # axes[1].legend()
    axes[1].set_xlim(-xbound,xbound)
    if norm_size==True:
        axes[0].set_ylabel('Normalized cell area (px)')
        axes[0].set_ylim(-.1,1.1)
        axes[1].set_ylim(-.1,1.1)
        axes[1].set_ylabel('Normalized cell area (px)')    
    fig.tight_layout()

    return 

#%% cc organelle analysis
def cc_org_analysis(cell_dfpaths: list, pre_org_dfpaths: list, vol_frac=False) -> (pd.DataFrame,pd.DataFrame):
    """
    Args: Lists of paths to csv files containing extracted cell and organelle information respectively.
    Outputs: Two dataframes containing selected cell and organelle metrics, and relative cc transition frames (G1/S 
            for G1 df, M/G1 for S/G2/M df). Optional plots
    """
    g1_dfs, s_dfs, s_b = [],[],[]
    for i in range(len(cell_dfpaths)):
        g1_df,s_df=cc_sort(cell_dfpaths[i]) #prelim sort and process
        g1_mothers, s_mothers, s_buds = find_cc_transitions(g1_df, s_df) #find annotated cc transitions
        path_parts=cell_dfpaths[i].stem.split("_")
        date=path_parts[1]
        fov=path_parts[4]
        fov_org_csvs=[p for p in pre_org_dfpaths if p.stem.split("_")[1]==date and p.stem.split("_")[4][:4]==fov] #parse out org csvs for this fov

        for cell_id in np.unique(g1_mothers['Cell_ID'].values): #for each cell...
            transitions=g1_mothers.loc[g1_mothers.Cell_ID==cell_id, 'Transitions'].values[0]
            pre_start_offset=g1_mothers.loc[g1_mothers.Cell_ID==cell_id, 'Start_Index'].values[0]
            if transitions==1:
                post_label="post_start_offset"
                post_offset=g1_mothers.loc[g1_mothers.Cell_ID==cell_id, 'Start_Index'].values[-1]
            elif transitions==2:
                post_label="post_div_offset"
                post_offset=g1_mothers.loc[g1_mothers.Cell_ID==cell_id, 'Div_Index'].values[-1]
            else:
                continue
            cell_metrics={
                "idx_cell"          : [cell_id],
                "date"              : [date], 
                "fov"               : [fov],
                "relationship"      : ["mother"],
                "transitions"       : [transitions],
                "pre_start_offset"  : [pre_start_offset],
                post_label          : [post_offset]
            }
            result=cell_metrics

            for org_path in fov_org_csvs: #for each organelle...
                pre_org_df=pd.read_csv(org_path) #read in organelle dfs
                post_org_df=pd.read_csv(Path(expmt_path+'/org_measure')/f"{org_path.stem[:-3]}post.csv")
                meta=parse_meta_organelle(org_path.stem) 
                org_label=meta['organelle']
                pre_org_vol=pre_org_df.loc[pre_org_df.idx_cell==cell_id, 'volume-pixel']
                if len(pre_org_vol)>0: #handle non-existent pre image case
                    if len(pre_org_vol)>1: #handle fragmented organelle case
                        pre_org_vol=np.sum(pre_org_vol)
                    else:
                        pre_org_vol=pre_org_vol.values[0]
                    pre_cell_vol=g1_mothers.loc[g1_mothers.Cell_ID==cell_id, 'cell_vol_vox'].values[0]
                    if vol_frac==True:
                        pre_org_vol=pre_org_vol/pre_cell_vol
                else:
                    pre_org_vol=None
                
                post_org_vol=post_org_df.loc[post_org_df.idx_cell==cell_id, 'volume-pixel']
                if len(post_org_vol)>0: #handle non-existent post image case
                    if len(post_org_vol)>1: #handle fragmented organelle case
                        post_org_vol=np.sum(post_org_vol)
                    else:
                        post_org_vol=post_org_vol.values[0]
                    post_cell_vol=g1_mothers.loc[g1_mothers.Cell_ID==cell_id, 'cell_vol_vox'].values[-1]
                    if vol_frac==True:
                        post_org_vol=post_org_vol/post_cell_vol
                else:
                    post_org_vol=None
                if type(post_org_vol)==np.float64 and type(pre_org_vol)==np.float64:
                    vol_difference=post_org_vol-pre_org_vol
                else:
                    vol_difference=None
                org_metrics={
                    org_label+"_vol_pre"   : [pre_org_vol],
                    org_label+"_vol_post"  : [post_org_vol]#,
                    # org_label+"_vol_diff": [vol_difference]
                }
            # result=meta | cell_metrics | org_metrics
                result=result | org_metrics
            g1_dfs.append(pd.DataFrame(result))

        for cell_id in np.unique(s_mothers['Cell_ID'].values):
            pre_div_offset=s_mothers.loc[s_mothers.Cell_ID==cell_id, 'Div_Index'].values[0]
            transitions=s_mothers.loc[s_mothers.Cell_ID==cell_id, 'Transitions'].values[0]
            if transitions==1:
                post_label="post_div_offset"
                post_offset=s_mothers.loc[s_mothers.Cell_ID==cell_id, 'Div_Index'].values[-1]
            elif transitions==2:
                post_label="post_start_offset"
                post_offset=s_mothers.loc[s_mothers.Cell_ID==cell_id, 'Start_Index'].values[-1]
            else:
                continue
            cell_metrics={
                "idx_cell"          : [cell_id],
                "date"              : [date], 
                "fov"               : [fov],
                "relationship"      : ["mother"],
                "transitions"       : [transitions],
                "pre_div_offset"    : [pre_div_offset],
                post_label          : [post_offset]
            }
            result=cell_metrics

            for org_path in fov_org_csvs: #for each organelle...
                pre_org_df=pd.read_csv(org_path) #read in organelle dfs
                post_org_df=pd.read_csv(Path(expmt_path+'/org_measure')/f"{org_path.stem[:-3]}post.csv")
                meta=parse_meta_organelle(org_path.stem) 
                org_label=meta['organelle']
                pre_org_vol=pre_org_df.loc[pre_org_df.idx_cell==cell_id, 'volume-pixel']
                if len(pre_org_vol)>0: #handle non-existent pre image case
                    if len(pre_org_vol)>1: #handle fragmented organelle case
                        pre_org_vol=np.sum(pre_org_vol)
                    else:
                        pre_org_vol=pre_org_vol.values[0]
                    pre_cell_vol=s_mothers.loc[s_mothers.Cell_ID==cell_id, 'cell_vol_vox'].values[0]
                    if vol_frac==True:
                        pre_org_vol=pre_org_vol/pre_cell_vol
                else:
                    pre_org_vol=None
                    
                post_org_vol=post_org_df.loc[post_org_df.idx_cell==cell_id, 'volume-pixel']
                if len(post_org_vol)>0: #handle non-existent post image case
                    if len(post_org_vol)>1: #handle fragmented organelle case
                        post_org_vol=np.sum(post_org_vol)
                    else:
                        post_org_vol=post_org_vol.values[0]
                    post_cell_vol=s_mothers.loc[s_mothers.Cell_ID==cell_id, 'cell_vol_vox'].values[0]
                    if vol_frac==True:
                        post_org_vol=post_org_vol/post_cell_vol
                else:
                    post_org_vol=None
                if type(post_org_vol)==np.float64 and type(pre_org_vol)==np.float64:
                    vol_difference=post_org_vol-pre_org_vol
                else:
                    vol_difference=None
                org_metrics={
                    org_label+"_vol_pre"   : [pre_org_vol],
                    org_label+"_vol_post"  : [post_org_vol]#,
                    # org_label+"_vol_diff": [vol_difference]
                }
                result=result | org_metrics
            s_dfs.append(pd.DataFrame(result))

        # for cell_id in np.unique(s_buds['Cell_ID'].values):
        #     pre_div_offset=s_buds.loc[s_buds.Cell_ID==cell_id, 'Relative_Index'].values[0]
        #     post_div_offset=s_buds.loc[s_buds.Cell_ID==cell_id, 'Relative_Index'].values[-1]+1
            
        #     pre_org_vol=pre_org_df.loc[pre_org_df.idx_cell==cell_id, 'volume-pixel']
        #     if len(pre_org_vol)>0: #handle non-existent pre image case
        #         if len(pre_org_vol)>1: #handle fragmented organelle case
        #             pre_org_vol=np.sum(pre_org_vol)
        #         else:
        #             pre_org_vol=pre_org_vol.values[0]
        #         pre_cell_vol=s_buds.loc[s_buds.Cell_ID==cell_id, 'cell_vol_vox'].values[0]
        #         if vol_frac==True:
        #             pre_org_vol=pre_org_vol/pre_cell_vol
        #     else:
        #         pre_org_vol=None
                
        #     post_org_vol=post_org_df.loc[post_org_df.idx_cell==cell_id, 'volume-pixel']
        #     if len(post_org_vol)>0: #handle non-existent post image case
        #         if len(post_org_vol)>1: #handle fragmented organelle case
        #             post_org_vol=np.sum(post_org_vol)
        #         else:
        #             post_org_vol=post_org_vol.values[0]
        #         post_cell_vol=s_buds.loc[s_buds.Cell_ID==cell_id, 'cell_vol_vox'].values[0]
        #         if vol_frac==True:
        #             post_org_vol=post_org_vol/post_cell_vol
        #     else:
        #         post_org_vol=None
        #     if type(post_org_vol)==np.float64 and type(pre_org_vol)==np.float64:
        #         vol_difference=post_org_vol-pre_org_vol
        #     else:
        #         vol_difference=None
        #     metrics={
        #         "idx_cell": [cell_id],
        #         "relationship": ['bud'],
        #         "pre_org_vol": [pre_org_vol],
        #         "post_org_vol": [post_org_vol],
        #         "pre_offset": [pre_div_offset],
        #         "post_offset": [post_div_offset],
        #         "vol_diff": [vol_difference]
        #     }
        #     result=meta | metrics
        #     s_dfs.append(pd.DataFrame(result))

    g1_org_df=pd.concat(g1_dfs, ignore_index=True)
    s_org_df=pd.concat(s_dfs, ignore_index=True)
    
    return g1_org_df, s_org_df

def plot_ccorg(g1_org_df: pd.DataFrame, s_org_df: pd.DataFrame, save_fig=False, xbound: int=30,ms=3):
    """
    Args: Two dataframes for either G1 or S/G2/M. Bool to save figure, x-axis bound param, set markersize.
    Outputs: Plots of org metric vs cc positions for either df. Option to save plots.
    """ 
    s_moms=s_org_df.loc[s_org_df.relationship=='mother']
    s_buds=s_org_df.loc[s_org_df.relationship=='bud']
    #Plot desired organelle metric as function of cc position wrt reference cc checkpoint for both groups.
    fig,axes=plt.subplots(nrows=len(orgs),ncols=2)
    for i in range(len(orgs)): #should i keep the four offset columns or distinguish using transitions count in here?
        org_label=orgs[i]
        #all the pre's
        axes[i, 0].errorbar(g1_org_df['pre_start_offset'].values, g1_org_df[org_label+"_vol_pre"].values, yerr=org_err[org_label]*g1_org_df[org_label+"_vol_pre"].values, ls='none', c='m',marker='o',ms=ms)
        axes[i, 1].errorbar(s_moms['pre_div_offset'].values, s_moms[org_label+"_vol_pre"].values, yerr=org_err[org_label]*s_moms[org_label+"_vol_pre"].values,ls='none', c='m',marker='o',ms=ms) 
        #then then the 1-transitions
        axes[i, 0].errorbar(g1_org_df.loc[g1_org_df.transitions==1,'post_start_offset'].values,g1_org_df.loc[g1_org_df.transitions==1,org_label+"_vol_post"].values, yerr=org_err[org_label]*g1_org_df.loc[g1_org_df.transitions==1,org_label+"_vol_post"].values,ls='none', c='m',marker='o',ms=ms)
        axes[i, 1].errorbar(s_moms.loc[s_moms.transitions==1,'post_div_offset'].values, s_moms.loc[s_moms.transitions==1,org_label+"_vol_post"].values, yerr=org_err[org_label]*s_moms.loc[s_moms.transitions==1,org_label+"_vol_post"].values,ls='none', c='m',marker='o',ms=ms)
        #then the 2-transitions
        axes[i, 1].errorbar(g1_org_df.loc[g1_org_df.transitions==2,'post_div_offset'].values,g1_org_df.loc[g1_org_df.transitions==2,org_label+"_vol_post"].values, yerr=org_err[org_label]*g1_org_df.loc[g1_org_df.transitions==2,org_label+"_vol_post"].values,ls='none', c='m',marker='o',ms=ms)
        axes[i, 0].errorbar(s_moms.loc[s_moms.transitions==2,'post_start_offset'].values, s_moms.loc[s_moms.transitions==2,org_label+"_vol_post"].values, yerr=org_err[org_label]*s_moms.loc[s_moms.transitions==2,org_label+"_vol_post"].values,ls='none', c='m',marker='o',ms=ms)
       
        # axes[i, 1].errorbar(s_buds['pre_offset'].values, s_buds[org_label+"_vol_pre"].values, yerr=org_err[org_label]*s_buds[org_label+"_vol_pre"].values,ls='none', c='r',marker='s') 
        # axes[i, 1].errorbar(s_buds['post_offset'].values, s_buds[org_label+"_vol_post"].values, yerr=org_err[org_label]*s_buds[org_label+"_vol_post"].values,ls='none', c='r',marker='s')

        axes[i, 0].set_title('G1 mother '+org_label+' vs CC position')
        axes[i, 0].set_xlabel('Frames relative to Start (5 min interval)')
        # axes[i, 0].set_ylabel(org_label+' volume (voxels)')
        axes[i, 0].set_ylabel(org_label+' volume fraction')
        axes[i, 0].set_xlim(-xbound,xbound)

        axes[i, 1].set_title('S/G2/M '+org_label+' vs CC position')
        axes[i, 1].set_xlabel('Frames relative to Division (5 min interval)')
        # axes[i, 1].set_ylabel(org_label+' volume (voxels)')
        axes[i, 1].set_ylabel(org_label+' volume fraction')    
        axes[i, 1].set_xlim(-xbound,xbound)

        fig.tight_layout()

    if save_fig==True:
        if vol_frac==True:
            plt.savefig(Path(expmt_path+'/cc_measure')/f"{org_label}_cc-analysis_normalized.png")
        else:
            plt.savefig(Path(expmt_path+'/cc_measure')/f"{org_label}_cc-analysis.png")
    return None

# %% Bin organelle measurements by cc postion
def binned_org(g1_org_df: pd.DataFrame, s_org_df: pd.DataFrame,plot_graph: bool=True):
    """
    Args: Dataframes of extracted cell and organelle metrics for each cc group. Option to plot.
    Outputs: Graphs displaying org vs cc trends, binned along the cc position axis.
    """
    all_cells=pd.concat([g1_org_df,s_org_df], ignore_index=True)
    if plot_graph==True:
        plt_ind=0
        fig,axes=plt.subplots(nrows=len(orgs),ncols=2,sharex=True)

    scores={}
    for org in orgs:
        start_out=[]
        for time in np.unique(all_cells['pre_start_offset'].values):
            bin_vals=all_cells.loc[all_cells.pre_start_offset==time,org+"_vol_pre"].values
            bin_vals=[entry for entry in bin_vals if ~np.isnan(entry)]
            if len(bin_vals)>0:
                bin_avg=np.mean(bin_vals)
                bin_std=np.std(bin_vals)
                start_out.append([time,bin_avg,bin_std])
        for time in np.unique(all_cells['post_start_offset'].values):
            bin_vals=all_cells.loc[all_cells.post_start_offset==time,org+"_vol_post"].values
            bin_vals=[entry for entry in bin_vals if ~np.isnan(entry)]
            if len(bin_vals)>0:
                bin_avg=np.mean(bin_vals)
                bin_std=np.std(bin_vals)
                start_out.append([time,bin_avg,bin_std])

        div_out=[]
        for time in np.unique(all_cells['pre_div_offset'].values):
            bin_vals=all_cells.loc[all_cells.pre_start_offset==time,org+"_vol_pre"].values
            bin_vals=[entry for entry in bin_vals if ~np.isnan(entry)]
            if len(bin_vals)>0:
                bin_avg=np.mean(bin_vals)
                bin_std=np.std(bin_vals)
                div_out.append([time,bin_avg,bin_std])
        for time in np.unique(all_cells['post_div_offset'].values):
            bin_vals=all_cells.loc[all_cells.post_start_offset==time,org+"_vol_post"].values
            bin_vals=[entry for entry in bin_vals if ~np.isnan(entry)]
            if len(bin_vals)>0:
                bin_avg=np.mean(bin_vals)
                bin_std=np.std(bin_vals)
                div_out.append([time,bin_avg,bin_std])

        # g1_out=[]
        # for time in np.unique(g1_org_df['pre_offset'].values):
        #     bin_vals=g1_org_df.loc[g1_org_df.pre_offset==time,org+"_vol_pre"].values
        #     bin_vals=[entry for entry in bin_vals if ~np.isnan(entry)]
        #     if len(bin_vals)>0:
        #         bin_avg=np.mean(bin_vals)
        #         bin_std=np.std(bin_vals)
        #         g1_out.append([time,bin_avg,bin_std])
        # for time in np.unique(g1_org_df['post_offset'].values):
        #     bin_vals=g1_org_df.loc[g1_org_df.post_offset==time,org+"_vol_post"].values
        #     bin_vals=[entry for entry in bin_vals if ~np.isnan(entry)]
        #     if len(bin_vals)>0:
        #         bin_avg=np.mean(bin_vals)
        #         bin_std=np.std(bin_vals)
        #         g1_out.append([time,bin_avg,bin_std])

        # s_out=[]
        # for time in np.unique(s_org_df['pre_offset'].values):
        #     bin_vals=s_org_df.loc[s_org_df.pre_offset==time,org+"_vol_pre"].values
        #     bin_vals=[entry for entry in bin_vals if ~np.isnan(entry)]
        #     if len(bin_vals)>0:
        #         bin_avg=np.mean(bin_vals)
        #         bin_std=np.std(bin_vals)
        #         s_out.append([time,bin_avg,bin_std])
        # for time in np.unique(s_org_df['post_offset'].values):
        #     bin_vals=s_org_df.loc[s_org_df.post_offset==time,org+"_vol_post"].values
        #     bin_vals=[entry for entry in bin_vals if ~np.isnan(entry)]
        #     if len(bin_vals)>0:
        #         bin_avg=np.mean(bin_vals)
        #         bin_std=np.std(bin_vals)
        #         s_out.append([time,bin_avg,bin_std])

        
        start_out,div_out=np.array(start_out),np.array(div_out)
        # g1_out,s_out=np.array(g1_out),np.array(s_out)

        #regression analysis here ###############
        # g1_scores=linreg(g1_out[:,0:2])
        # s_scores=linreg(s_out[:,0:2])
        # score={
        #     org:[g1_scores,s_scores]
        # }
        # scores=scores|score

        if plot_graph==True:
            axes[plt_ind,0].errorbar(start_out[:,0],start_out[:,1],yerr=start_out[:,2],ls='none', c='m',marker='o')
            axes[plt_ind,1].errorbar(div_out[:,0],div_out[:,1],yerr=div_out[:,2],ls='none', c='m',marker='o')
            # axes[plt_ind,0].errorbar(g1_out[:,0],g1_out[:,1],yerr=g1_out[:,2],ls='none', c='m',marker='o')
            # axes[plt_ind,1].errorbar(s_out[:,0],s_out[:,1],yerr=s_out[:,2],ls='none', c='m',marker='o')
            
            axes[plt_ind,0].set_title('G1 mother '+org+' vs CC position')
            # axes[plt_ind,0].set_xlabel('Frames relative to Start (5 min interval)')
            # axes[plt_ind,0].set_ylabel('Binned'+org+' volume')
            axes[plt_ind,0].set_xlim(-xbound,xbound)
            # axes[plt_ind,0].set_ylim(-.05,1)
            axes[plt_ind,1].set_title('S/G2/M mother '+org+' vs CC position')
            # axes[plt_ind,1].set_xlabel('Frames relative to Division (5 min interval)')
            # axes[plt_ind,1].set_ylabel('Binned'+org+' volume')
            axes[plt_ind,1].set_xlim(-xbound,xbound)
            # axes[plt_ind,1].set_ylim(-.05,1)
            fig.tight_layout()
            plt_ind+=1
        # axes[plt_ind-1,0].set_xlabel('Frames relative to Start (5 min interval)')        
        # axes[plt_ind-1,1].set_xlabel('Frames relative to Division (5 min interval)')
    if plot_graph==True:
        fig.supxlabel('Frames relative to Annotated Bud Emergence/Division (5 min interval')
        fig.supylabel('Binned organelle volume')
        fig.tight_layout()
    # return scores
    # return g1_out,s_out
    return start_out,div_out
# %% Myo1 dynamics analysis
def ticks(img_df: pd.DataFrame, thresh: float=0.5):
    """
    Args: Dataframe containing extracted roi info to be analyzed.
    Outputs: Dataframe containing up and/or down tick frames for each unique roi. 
    """
    dfs=[]
    for roi_id in np.unique(img_df['ROI_ID'].values):
        vals=img_df.loc[img_df.ROI_ID==roi_id, 'sum'].values
        vals=vals/np.max(vals)
        yi=img_df.loc[img_df.ROI_ID==roi_id, 'y'].values[0]
        xi=img_df.loc[img_df.ROI_ID==roi_id, 'x'].values[0]
        frames=img_df.loc[img_df.ROI_ID==roi_id, 'frame'].values
        on_indices=np.where(vals>thresh)[0]
        off_indices=np.where(vals<thresh)[0]
        start_frame=None
        div_frame=None
        # if vals[0]<thresh:
        pre_on=[x for x in off_indices if x<on_indices[0]]
        # print(len(pre_on))
        if len(pre_on)>0:
            start_frame=frames[pre_on[-1]]+1  
            # print(start_frame)
        # elif vals[-1]<thresh:
        post_on=[x for x in off_indices if x>on_indices[-1]]
        if len(post_on)>0:
            div_frame=frames[post_on[0]]
        # else:
        #     continue

        metrics={
            "ROI_ID"      : [roi_id],
            "y"           : [yi],
            "x"           : [xi],
            "Start_Frame" : [start_frame],
            "Div_Frame"   : [div_frame]
        }

        roi_df=pd.DataFrame(metrics)
        dfs.append(roi_df)
    output_df=pd.concat(dfs, ignore_index=True)
    return output_df

def diff_dist(img_df: pd.DataFrame):
    """
    Args: Dataframe containing myo1 roi intensity information.
    Outputs: 1D array containing distribution of changes in roi intensity.
    """
    diffs=[]
    for roi_id in np.unique(img_df['ROI_ID'].values):
        vals=img_df.loc[img_df.ROI_ID==roi_id, 'sum'].values
        vals=vals/np.max(vals)
        diffs.append(np.diff(vals))
    return np.concatenate(diffs)

def sum_dist(img_df: pd.DataFrame):
    """
    Args: Dataframe containing myo1 roi intensity information.
    Outputs: 1D array containing distribution of roi intensities.
    """
    sums=[]
    for roi_id in np.unique(img_df['ROI_ID'].values):
        vals=img_df.loc[img_df.ROI_ID==roi_id, 'sum'].values
        vals=vals/np.max(vals)
        sums.append(vals)
    return np.concatenate(sums)
# %% Autoassign myo1 roi and cellacdc IDs
from skimage import io, draw
from scipy import stats
def autoassign(myo1df_path: str, cellmask_path: str):
    """
    Args: Dataframe containing detected myo1 ROI coords; cell segmentation mask file path.
    Outputs: Dataframe containing ROI-Cell_ID assignments.
    """
    mask=io.imread(cellmask_path)
    imgdf=pd.read_csv(myo1df_path)
    rowsout=[]
    missed=0
    for roi_id in np.unique(imgdf['ROI_ID'].values):
        # print(roi_id)
        roi_rows=imgdf.loc[imgdf.ROI_ID==roi_id].copy()
        frames=roi_rows['frame'].values
        for frame_i in frames:
            yi=roi_rows.loc[roi_rows.frame==frame_i,'y'].values[0]
            xi=roi_rows.loc[roi_rows.frame==frame_i,'x'].values[0]
            rr, cc = draw.disk((yi, xi), 8, shape=(mask.shape[1],mask.shape[2]))
            unique, counts = np.unique(mask[frame_i,rr,cc][mask[frame_i,rr,cc]>0], return_counts=True) 
            # cell_id = stats.mode(mask[frame_i,rr, cc])
            if len(counts)>0:
                cell_id=unique[np.argmax(counts)]
                # print(yi,xi,frame_i, cell_id)
                # cell_id=mask[frame_i, int(yi), int(xi)] #f,y,x order in ImageJ and io.imread(). y=0 at top of image.
                roi_rows.loc[roi_rows.frame==frame_i, 'Cell_ID']=cell_id
            else:
                missed+=1
                continue
        rowsout.append(roi_rows)
    return pd.concat(rowsout)
# %% cc size and myo1 overlay
# def overlay_myo1size(cell_path: str, myo1_path: str):
#     """
#     Args: Paths to dataframes containing cell size, myo1 intensity, and cc phase info
#     Outputs:  
#     """
#     return
# %% finding cc differences
import math
def find_ccdiffs(cell_path: str, myo1roi_path: str, bend_path: str):
    """
    Args: Paths to cell segm output csv, myo1 img csv, manual bend annotation csv, and ROI/Cell_ID alignment csv
    Outputs: Tuple of arrays containing the frame differences in cc progression between myo1 signal and manual annotations
    """
    # acdcdf=pd.read_csv(cell_path)
    myo1df=pd.read_csv(myo1roi_path)
    benddf=pd.read_csv(bend_path)
    ticksdf=ticks(myo1df)
    g1_df,s_df=cc_sort(cell_path)
    g1_mothers, s_mothers, s_buds = find_cc_transitions(g1_df, s_df) #find some way to not leave S to G1 to budding mothers on the table
    # print(len(np.unique(g1_mothers['Cell_ID'].values)),len(np.unique(s_mothers['Cell_ID'].values)))
    myo1start, myo1div, myo1bend, divbend =[],[],[],[]
    count=0
    nov0525_ignore=[490,409,472,192,323,304,434,468,468,482,524,453,465]
    for cell_id in np.unique(g1_mothers['Cell_ID'].values):
        if cell_id not in nov0525_ignore:
            if cell_id in myo1df['Cell_ID'].values:
                if (roi_id:=myo1df.loc[myo1df.Cell_ID==cell_id, 'ROI_ID'].values[0]) in ticksdf['ROI_ID'].values:
                    count+=1
                    cell_rows=g1_mothers.loc[g1_mothers.Cell_ID==cell_id].copy()
                    annot_start=cell_rows.loc[cell_rows.Relative_Index==0, 'frame_i'].values[0]
                    myo1_start=ticksdf.loc[ticksdf.ROI_ID==roi_id, 'Start_Frame'].values[0]
                    if myo1_start!=None:
                        myo1start.append(annot_start-myo1_start)
                        # if annot_start-myo1_start>4:
                        #     print(roi_id, cell_id)
                        #     print(myo1_start, annot_start)
                    # else:
                    #     print(ticksdf.loc[ticksdf.ROI_ID==roi_id, 'Div_Frame'].values[0])
                    #     continue

    for cell_id in np.unique(s_mothers['Cell_ID'].values):
        if cell_id not in nov0525_ignore:
            if math.isnan(s_mothers.loc[s_mothers.Cell_ID==cell_id, 'Start_Index'].values[0])==False:
                cell_rows=s_mothers.loc[s_mothers.Cell_ID==cell_id].copy()
                annot_start=cell_rows.loc[cell_rows.Start_Index==0, 'frame_i'].values[0]
                if cell_id in myo1df['Cell_ID'].values and (roi_id:=myo1df.loc[myo1df.Cell_ID==cell_id, 'ROI_ID'].values[0]) in ticksdf['ROI_ID'].values:
                    # count+=1
                    roi_id=myo1df.loc[myo1df.Cell_ID==cell_id, 'ROI_ID'].values[0]
                    myo1_start=ticksdf.loc[ticksdf.ROI_ID==roi_id, 'Start_Frame'].values[0]
                    if myo1_start!=None:
                        myo1start.append(annot_start-myo1_start)
                        # if annot_start-myo1_start>4:
                        #     print(roi_id, cell_id)
                        #     print(myo1_start, annot_start)

            if math.isnan(s_mothers.loc[s_mothers.Cell_ID==cell_id, 'Div_Index'].values[0])==False:
                cell_rows=s_mothers.loc[s_mothers.Cell_ID==cell_id].copy()
                annot_div=cell_rows.loc[cell_rows.Div_Index==0, 'frame_i'].values[0]
                if cell_id in myo1df['Cell_ID'].values and (roi_id:=myo1df.loc[myo1df.Cell_ID==cell_id, 'ROI_ID'].values[0]) in ticksdf['ROI_ID'].values:
                    # count+=1
                    roi_id=myo1df.loc[myo1df.Cell_ID==cell_id, 'ROI_ID'].values[0]
                    myo1_div=ticksdf.loc[ticksdf.ROI_ID==roi_id, 'Div_Frame'].values[0]
                    if myo1_div!=None:
                        myo1div.append(annot_div-myo1_div)
                        # if annot_div-myo1_div==1:
                            # print(roi_id, cell_id)
                            # print(myo1_div, annot_div)
                        if cell_id in np.unique(benddf.loc[benddf.relationship=='mother', 'Cell_ID'].values):
                            bend_frame=benddf.loc[benddf.Cell_ID==cell_id, 'frame_b'].values[0]-1
                            myo1bend.append(bend_frame-myo1_div)
                    # else:
                    #     continue

                if cell_id in np.unique(benddf.loc[benddf.relationship=='mother', 'Cell_ID'].values):
                    bend_frame=benddf.loc[benddf.Cell_ID==cell_id, 'frame_b'].values[0]-1
                    divbend.append(annot_div-bend_frame)

    # print(len(myo1start), len(myo1div), len(myo1bend), len(divbend))
    # print(count)
    return myo1start, myo1div, myo1bend, divbend
    # return g1_mothers

# def plot_distr():
#     plt.figure()
#     plt.xlabel('Frame # (5 min interval)')
#     plt.title('Relative difference between myo1 on/off and division annotations')
#     n0,b0,p0=plt.hist(md,bins='doane', alpha=0.75, edgecolor='black', color='b', label='Divisions n=('+f'{len(md)})')
#     n1,b1,p1=plt.hist(mb, bins=b0, alpha=0.75, edgecolor='black', color='y', label='Bends n=('+f'{len(mb)})')
#     plt.legend()
    # plt.figure()
    # plt.xlabel('# Frames after myo1 off (5 min interval)')
    # plt.title('Relative difference between myo1 off and division annotations')
    # n0,b0,p0=plt.hist(md,bins=np.arange(-.5,7.5,1), alpha=0.75, edgecolor='black', color='b', label='Division (n='+f'{len(ms)})')
    # n1,b1,p1=plt.hist(mb, bins=b0, alpha=0.75, edgecolor='black', color='y', label='Bends (n='+f'{len(mb)})')
    # plt.legend()
#     return

# %%
list_cell = [] #Initialize lists for cell and org csvs
list_in   = [] 

if not os.path.exists(newpath:=Path(expmt_path+'/cc_measure')):
    print('Creating folder ',str(expmt_path+'/cc_measure'))
    os.makedirs(newpath)

for path_in in Path(expmt_path+'/org_measure').glob('*pre.csv'):
    path_parts=path_in.stem.split("_")
    fov=path_parts[4][:4] ########HARD-CODED BELOW = BAD
    # path_cell=Path(expmt_path+'/cell_measure')/f"BF-timelapse_03042026_eyrbow_glucose-2.0_{fov}_acdc_output_cpsam_d10.csv"

    # path_end="_".join(path_parts[:-1])
    cell_parts=path_in.stem.split('-')
    cell_end="-".join(cell_parts[:3])[3:]
    path_acdc=Path(expmt_path+'/cell_measure')/f"BF-timelapse_{cell_end}_acdc_output_cpsam-d25.csv"

    list_in.append(path_in)
    list_cell.append(path_acdc)
# list_in=np.unique(list_in)
# list_cell=np.unique(list_cell)

# args = pd.DataFrame({
#     "path_in":   list_in,
#     "path_cell": list_cell
# })

# %%
