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
###TREATING AS 10% ACROSS THE BOARD FOR NOW UNTIL UPDATE FROM SIMON'S PAPER
org_err={
"er":0.1,
"px":0.1,
"vo":0.1,
"mt":0.1,
"gl":.1,
"ld":.1
}
xbound=30
# %%
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
# %% find cc transitions
def find_cc_transitions(g1_df: pd.DataFrame, s_df: pd.DataFrame):
    """
    Args: cc_sort output Dataframes containing cells in G1 and S/G2/M phases respectively at frame 0.
    Outputs: Same dataframes with new index column relative to cc checkpoint progression.
                Currently dropping g1 but not s buds###############
    """
    g1_dfs=[]
    for cell_id in np.unique(g1_df.loc[g1_df.relationship=='mother','Cell_ID'].values): #for each G1 mother
        cell_rows=g1_df.loc[g1_df.Cell_ID==cell_id].copy()
        # if len(cell_rows.loc[cell_rows.cell_cycle_stage=='S', 'frame_i'].values>0): #if G1/S transition happens...
        alignment_frame=cell_rows.loc[cell_rows.cell_cycle_stage=='S', 'frame_i'].values[0]
        relative_index=[(cell_rows['frame_i'].values[i]-alignment_frame) for i in range(len(cell_rows))]
        cell_rows.loc[cell_rows.Cell_ID==cell_id, 'Relative_Index']=relative_index
        g1_dfs.append(cell_rows)
    g1_mothers=pd.concat(g1_dfs)

    smom_dfs=[]
    sbud_dfs=[] ################### propagate mult cc transitions code
    for cell_id in np.unique(s_df['Cell_ID'].values): #for each S/G2/M cell
        cell_rows=s_df.loc[s_df.Cell_ID==cell_id].copy()
        # if len(cell_rows.loc[cell_rows.cell_cycle_stage=='G1', 'frame_i'])>0: #if M/G1 transition happens...
        div_frame=cell_rows.loc[cell_rows.cell_cycle_stage=='G1', 'frame_i'].values[0]
        div_index=[(cell_rows['frame_i'].values[i]-div_frame) for i in range(len(cell_rows))]
        cell_rows.loc[cell_rows.Cell_ID==cell_id, 'Div_Index']=div_index
        if cell_rows.loc[cell_rows.cell_cycle_stage=='S', 'frame_i'].values[-1]>div_frame:
            start_frame=cell_rows.loc[cell_rows.cell_cycle_stage=='G1', 'frame_i'].values[-1]+1
            start_index=[(cell_rows['frame_i'].values[i]-start_frame) for i in range(len(cell_rows))]
            cell_rows.loc[cell_rows.Cell_ID==cell_id, 'Start_Index']=start_index

        if cell_rows['relationship'].values[0]=='mother':
            smom_dfs.append(cell_rows)
        else:
            sbud_dfs.append(cell_rows)
    s_mothers=pd.concat(smom_dfs)
    s_buds=pd.concat(sbud_dfs)

    return g1_mothers, s_mothers, s_buds

# %% cc size analysis
# def plot_cc_scatter(df: pd.DataFrame, axis: int, org_label: str):
#     fig,axes=plt.subplots(nrows=1,ncols=2)

#     axes[axis].scatter(offset,vol_frac,c='m',marker='o') #scatter of organelle volume fraction versus relative frame index
#     axes[axis].errorbar(offset, vol_frac, yerr=org_err[org_label]*vol_frac, c='m',marker='o')

def cc_size_analysis(cell_dfpath: str, norm_size: bool=False, specify_range: bool=False, xbound: int=30):
    """
    Args: Path to cell metric csv file.
    Outputs: Graph displaying cell size profiles vs cell cycle position.
    """
    g1_df, s_df=cc_sort(cell_dfpath)
    g1_mothers,s_mothers,s_buds=find_cc_transitions(g1_df,s_df)
    g1_buds=g1_df.loc[g1_df['relationship']=='bud']

    # fig,axes=plt.subplots(nrows=1,ncols=2)
    for cell_id in np.unique(g1_df['Cell_ID'].values):
        if g1_df.loc[g1_df.Cell_ID==cell_id, 'relationship'].values[0]=='mother':
            axes[0].plot(g1_mothers.loc[g1_mothers.Cell_ID==cell_id, 'Relative Index'].values, g1_mothers.loc[g1_mothers.Cell_ID==cell_id,'cell_area_pxl'].values,'k-',label='mother')
        else:
            axes[0].plot(g1_buds.loc[g1_buds.Cell_ID==cell_id, 'frame_i'].values, g1_buds.loc[g1_buds.Cell_ID==cell_id,'cell_area_pxl'].values,'b--',label='bud') 
    for cell_id in np.unique(s_df['Cell_ID'].values):
        if s_df.loc[s_df.Cell_ID==cell_id, 'relationship'].values[0]=='mother':
            axes[1].plot(s_mothers.loc[s_mothers.Cell_ID==cell_id, 'Relative Index'].values, s_mothers.loc[s_mothers.Cell_ID==cell_id,'cell_area_pxl'].values, 'k-',label='mother') 
        else:
            axes[1].plot(s_buds.loc[s_buds.Cell_ID==cell_id, 'Relative Index'].values, s_buds.loc[s_buds.Cell_ID==cell_id,'cell_area_pxl'].values,'b--',label='bud') 
    
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

    return g1_df, g1_mothers ####################

#%% cc organelle analysis
def cc_org_analysis(cell_dfpaths: list, pre_org_dfpaths: list, plot_graph: bool=False, xbound: int=30, save_fig: bool=False, vol_frac=False) -> (pd.DataFrame,pd.DataFrame):
    """
    Args: Lists of paths to csv files containing extracted cell and organelle information respectively.
    Outputs: Two dataframes containing selected cell and organelle metrics, and relative cc transition frames (G1/S 
            for G1 df, M/G1 for S/G2/M df). Optional plot(s)
    """
    g1_dfs, s_dfs, s_b = [],[],[]
    for i in range(len(cell_dfpaths)):
        g1_df,s_df=cc_sort(cell_dfpaths[i])
        pre_org_df=pd.read_csv(pre_org_dfpaths[i]) 
        post_org_df=pd.read_csv(Path(expmt_path+'/org_measure')/f"{pre_org_dfpaths[i].stem[:-3]}post.csv")
        meta=parse_meta_organelle(pre_org_dfpaths[i].stem)
        g1_mothers, s_mothers, s_buds = find_cc_transitions(g1_df, s_df)
        org_label=meta['organelle']
      
        for cell_id in np.unique(g1_mothers['Cell_ID'].values):
            pre_start_offset=g1_mothers.loc[g1_mothers.Cell_ID==cell_id, 'Relative Index'].values[0]
            post_start_offset=g1_mothers.loc[g1_mothers.Cell_ID==cell_id, 'Relative Index'].values[-1]+1
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
            metrics={
                "idx_cell": [cell_id],
                "relationship": ["mother"],
                "pre_org_vol": [pre_org_vol],
                "post_org_vol": [post_org_vol],
                "pre_offset": [pre_start_offset],
                "post_offset": [post_start_offset],
                "vol_diff": [vol_difference]
            }
            result=meta | metrics
            g1_dfs.append(pd.DataFrame(result))

        for cell_id in np.unique(s_mothers['Cell_ID'].values):
            pre_div_offset=s_mothers.loc[s_mothers.Cell_ID==cell_id, 'Relative Index'].values[0]
            post_div_offset=s_mothers.loc[s_mothers.Cell_ID==cell_id, 'Relative Index'].values[-1]+1
            
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
            metrics={
                "idx_cell": [cell_id],
                "relationship": ['mother'],
                "pre_org_vol": [pre_org_vol],
                "post_org_vol": [post_org_vol],
                "pre_offset": [pre_div_offset],
                "post_offset": [post_div_offset],
                "vol_diff": [vol_difference]
            }
            result=meta | metrics
            s_dfs.append(pd.DataFrame(result))

        for cell_id in np.unique(s_buds['Cell_ID'].values):
            pre_div_offset=s_buds.loc[s_buds.Cell_ID==cell_id, 'Relative Index'].values[0]
            post_div_offset=s_buds.loc[s_buds.Cell_ID==cell_id, 'Relative Index'].values[-1]+1
            
            pre_org_vol=pre_org_df.loc[pre_org_df.idx_cell==cell_id, 'volume-pixel']
            if len(pre_org_vol)>0: #handle non-existent pre image case
                if len(pre_org_vol)>1: #handle fragmented organelle case
                    pre_org_vol=np.sum(pre_org_vol)
                else:
                    pre_org_vol=pre_org_vol.values[0]
                pre_cell_vol=s_buds.loc[s_buds.Cell_ID==cell_id, 'cell_vol_vox'].values[0]
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
                post_cell_vol=s_buds.loc[s_buds.Cell_ID==cell_id, 'cell_vol_vox'].values[0]
                if vol_frac==True:
                    post_org_vol=post_org_vol/post_cell_vol
            else:
                post_org_vol=None
            if type(post_org_vol)==np.float64 and type(pre_org_vol)==np.float64:
                vol_difference=post_org_vol-pre_org_vol
            else:
                vol_difference=None
            metrics={
                "idx_cell": [cell_id],
                "relationship": ['bud'],
                "pre_org_vol": [pre_org_vol],
                "post_org_vol": [post_org_vol],
                "pre_offset": [pre_div_offset],
                "post_offset": [post_div_offset],
                "vol_diff": [vol_difference]
            }
            result=meta | metrics
            s_dfs.append(pd.DataFrame(result))

    g1_org_df=pd.concat(g1_dfs, ignore_index=True) 
    s_org_df=pd.concat(s_dfs, ignore_index=True)
    #Plot desired organelle metric as function of cc position wrt reference cc checkpoint for both groups.
    if plot_graph==True:
        # fig,axes=plt.subplots(nrows=1,ncols=2)
        s_moms=s_org_df.loc[s_org_df.relationship=='mother']
        s_buds=s_org_df.loc[s_org_df.relationship=='bud']
        axes[0].errorbar(g1_org_df['pre_offset'].values, g1_org_df['pre_org_vol'].values, yerr=org_err[org_label]*g1_org_df['pre_org_vol'].values, ls='none', c='m',marker='o')
        axes[0].errorbar(g1_org_df['post_offset'].values,g1_org_df['post_org_vol'].values, yerr=org_err[org_label]*g1_org_df['post_org_vol'].values,ls='none', c='m',marker='o') #scatter of organelle volume fraction versus relative frame index
        axes[1].errorbar(s_moms['pre_offset'].values, s_moms['pre_org_vol'].values, yerr=org_err[org_label]*s_moms['pre_org_vol'].values,ls='none', c='m',marker='o') 
        axes[1].errorbar(s_moms['post_offset'].values, s_moms['post_org_vol'].values, yerr=org_err[org_label]*s_moms['post_org_vol'].values,ls='none', c='m',marker='o')
        axes[1].errorbar(s_buds['pre_offset'].values, s_buds['pre_org_vol'].values, yerr=org_err[org_label]*s_buds['pre_org_vol'].values,ls='none', c='r',marker='s') 
        axes[1].errorbar(s_buds['post_offset'].values, s_buds['post_org_vol'].values, yerr=org_err[org_label]*s_buds['post_org_vol'].values,ls='none', c='r',marker='s')

        axes[0].set_title('G1 mother '+org_label+' vs CC position')
        axes[0].set_xlabel('Frames relative to Start (5 min interval)')
        axes[0].set_ylabel(org_label+' volume (voxels)')
        axes[0].set_xlim(-xbound,xbound)
    
        axes[1].set_title('S/G2/M '+org_label+' vs CC position')
        axes[1].set_xlabel('Frames relative to Division (5 min interval)')
        axes[1].set_ylabel(org_label+' volume (voxels)')
        axes[1].set_xlim(-xbound,xbound)
        if vol_frac==True:
            axes[0].set_ylabel(org_label+' volume fraction')
            axes[0].set_ylim(-.1,1.1)
            axes[1].set_ylim(-.1,1.1)
            axes[1].set_ylabel(org_label+' volume fraction')    
        fig.tight_layout()

    if save_fig==True:
        if vol_frac==True:
            plt.savefig(Path(expmt_path+'/cc_measure')/f"{org_label}_cc-analysis_normalized.png")
        else:
            plt.savefig(Path(expmt_path+'/cc_measure')/f"{org_label}_cc-analysis.png")

    return g1_org_df, s_org_df

# %% Bin organelle measurements by cc postion
def binned_org(g1_ls: list,s_ls: list,plot_graph: bool=True):
    """
    Args: Lists containing Dataframes of extracted cell and organelle metrics for each cc group.
    Outputs: Graphs displaying org vs cc trends, binned along the cc position axis.
    """
    # g1_df,s_df=cc_org_analysis(cell_path,pre_org_path,plot_graph=False,save_fig=False,vol_frac=True)
    g1_df=pd.concat(g1_ls, ignore_index=True)
    s_df=pd.concat(s_ls, ignore_index=True)
    org_label=g1_df['organelle'][0]

    g1_out=[]
    for time in np.unique(g1_df['pre_offset'].values):
        bin_vals=g1_df.loc[g1_df.pre_offset==time,'pre_org_vol'].values
        bin_sum=np.sum(bin_vals)
        bin_avg=bin_sum/len(bin_vals)
        bin_std=np.std(bin_vals)
        g1_out.append([time,bin_avg,bin_std])
    for time in np.unique(g1_df['post_offset'].values):
        bin_vals=g1_df.loc[g1_df.post_offset==time,'post_org_vol'].values
        bin_sum=np.sum(bin_vals)
        bin_avg=bin_sum/len(bin_vals)
        bin_std=np.std(bin_vals)
        g1_out.append([time,bin_avg,bin_std])

    s_out=[]
    for time in np.unique(s_df['pre_offset'].values):
        bin_vals=s_df.loc[s_df.pre_offset==time,'pre_org_vol'].values
        bin_sum=np.sum(bin_vals)
        bin_avg=bin_sum/len(bin_vals)
        bin_std=np.std(bin_vals)
        s_out.append([time,bin_avg,bin_std])
    for time in np.unique(s_df['post_offset'].values):
        bin_vals=s_df.loc[s_df.post_offset==time,'post_org_vol'].values
        bin_sum=np.sum(bin_vals)
        bin_avg=bin_sum/len(bin_vals)
        bin_std=np.std(bin_vals)
        s_out.append([time,bin_avg,bin_std])

    g1_out,s_out=np.array(g1_out),np.array(s_out)

    if plot_graph==True:
        axes[0].errorbar(g1_out[:,0],g1_out[:,1],yerr=g1_out[:,2],ls='none', c='m',marker='o')
        axes[1].errorbar(s_out[:,0],s_out[:,1],yerr=s_out[:,2],ls='none', c='m',marker='o')
        
        axes[0].set_title('G1 mother '+org_label+' vs CC position')
        axes[0].set_xlabel('Frames relative to Start (5 min interval)')
        axes[0].set_ylabel('Binned'+org_label+' volume')
        axes[0].set_xlim(-xbound,xbound)
    
        axes[1].set_title('S/G2/M mother '+org_label+' vs CC position')
        axes[1].set_xlabel('Frames relative to Division (5 min interval)')
        axes[1].set_ylabel('Binned'+org_label+' volume')
        axes[1].set_xlim(-xbound,xbound)
        fig.tight_layout()  
    return g1_out,s_out

# %% Myo1 globular analysis
def ticks(img_df: pd.DataFrame, thresh: float=0.7):
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

        if vals[0]<thresh:
            start_indices=np.where(vals<thresh)[0]
            start_frame=frames[start_indices[-1]]+1  
            div_frame=None
        elif vals[-1]<thresh:
            end_indices=np.where(vals>thresh)[0]
            div_frame=frames[end_indices[-1]]+1
            start_frame=None
        else:
            continue

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

import skimage.io as io
def autoassign(myo1df_path: str, cellmask_path: str):
    """
    Args: Dataframe containing detected myo1 ROI coords; cell segmentation mask file path.
    Outputs: Dataframe containing ROI-Cell_ID assignments.
    """
    mask=io.imread(cellmask_path)
    imgdf=pd.read_csv(myo1df_path)
    rowsout=[]
    for roi_id in np.unique(imgdf['ROI_ID'].values):
        roi_rows=imgdf.loc[imgdf.ROI_ID==roi_id].copy()
        yi=roi_rows['y'].values[0]
        xi=roi_rows['x'].values[0]
        frames=roi_rows['frame'].values
        for frame_i in frames:
            #draw circle about roi centroid, take mode of mask pixel values excl zero.
            cell_id=mask[frame_i, int(yi), int(xi)] #f,y,x order in ImageJ and io.imread(). y=0 at top of image.
            roi_rows.loc[roi_rows.frame==frame_i, 'Cell_ID']=cell_id
        rowsout.append(roi_rows)
    return pd.concat(rowsout)

# %% finding cc differences
import math
def find_ccdiffs(acdc_path: str, myo1roi_path: str, bend_path: str, alignids_path: str):
    """
    Args: Paths to acdc output csv, myo1 img csv, manual bend annotation csv, and ROI/Cell_ID alignment csv
    Outputs: Tuple of arrays containing the frame differences in cc progression between myo1 signal and manual annotations
    """
    # acdcdf=pd.read_csv(acdc_path)
    myo1df=pd.read_csv(myo1roi_path)
    benddf=pd.read_csv(bend_path)
    aligndf=pd.read_csv(alignids_path)

    ticksdf=ticks(myo1df)

    g1_df,s_df=cc_sort(acdc_path)
    g1_mothers, s_mothers, s_buds = find_cc_transitions(g1_df, s_df) #find some way to not leave S to G1 to budding mothers on the table
    myo1start, myo1div, myo1bend, divbend =[],[],[],[]
    count=0
    for cell_id in np.unique(g1_mothers['Cell_ID'].values):
        if cell_id in aligndf['mother_ID'].values and (roi_id:=aligndf.loc[aligndf.mother_ID==cell_id, 'ROI_ID'].values[0]) in ticksdf['ROI_ID'].values:
            count+=1
            cell_rows=g1_mothers.loc[g1_mothers.Cell_ID==cell_id].copy()
            annot_start=cell_rows.loc[cell_rows.Relative_Index==0, 'frame_i'].values[0]
            # roi_id=aligndf.loc[aligndf.mother_ID==cell_id, 'ROI_ID'].values[0]
            myo1_start=ticksdf.loc[ticksdf.ROI_ID==roi_id, 'Start_Frame'].values[0]
            if myo1_start!=None:
                myo1start.append(myo1_start-annot_start)
                print(roi_id, cell_id)
                print(myo1_start, annot_start)
            else:
                continue

    for cell_id in np.unique(s_mothers['Cell_ID'].values):
        if math.isnan(s_mothers.loc[s_mothers.Cell_ID==cell_id, 'Start_Index'].values[0])==False:
            cell_rows=s_mothers.loc[s_mothers.Cell_ID==cell_id].copy()
            annot_start=cell_rows.loc[cell_rows.Start_Index==0, 'frame_i'].values[0]
            if cell_id in aligndf['mother_ID'].values and (roi_id:=aligndf.loc[aligndf.mother_ID==cell_id, 'ROI_ID'].values[0]) in ticksdf['ROI_ID'].values:
                # count+=1
                roi_id=aligndf.loc[aligndf.mother_ID==cell_id, 'ROI_ID'].values[0]
                myo1_start=ticksdf.loc[ticksdf.ROI_ID==roi_id, 'Start_Frame'].values[0]
                if myo1_start!=None:
                    myo1start.append(myo1_start-annot_start)

        if math.isnan(s_mothers.loc[s_mothers.Cell_ID==cell_id, 'Div_Index'].values[0])==False:
            cell_rows=s_mothers.loc[s_mothers.Cell_ID==cell_id].copy()
            annot_div=cell_rows.loc[cell_rows.Div_Index==0, 'frame_i'].values[0]
            if cell_id in aligndf['mother_ID'].values and (roi_id:=aligndf.loc[aligndf.mother_ID==cell_id, 'ROI_ID'].values[0]) in ticksdf['ROI_ID'].values:
                # count+=1
                roi_id=aligndf.loc[aligndf.mother_ID==cell_id, 'ROI_ID'].values[0]
                myo1_div=ticksdf.loc[ticksdf.ROI_ID==roi_id, 'Div_Frame'].values[0]
                if myo1_div!=None:
                    myo1div.append(myo1_div-annot_div)
                    if cell_id in np.unique(benddf.loc[benddf.relationship=='mother', 'Cell_ID'].values):
                        bend_frame=benddf.loc[benddf.Cell_ID==cell_id, 'frame_b'].values[0]-1
                        myo1bend.append(myo1_div-bend_frame)
                else:
                    continue

            if cell_id in np.unique(benddf.loc[benddf.relationship=='mother', 'Cell_ID'].values):
                bend_frame=benddf.loc[benddf.Cell_ID==cell_id, 'frame_b'].values[0]-1
                divbend.append(annot_div-bend_frame)

    # print(len(myo1start), len(myo1div), len(myo1bend), len(divbend))
    # print(count, len(aligndf['mother_ID'].values))
    return myo1start, myo1div, myo1bend, divbend
    # return g1_mothers
# %%
# def extract_trend():
#     """
#     Args:
#     Outputs:
#     """

#     return trend_array 
# %%
list_in   = [] #Initialize lists for org and cell csvs
list_cell = []

if not os.path.exists(newpath:=Path(expmt_path+'/cc_measure')):
    print('Creating folder ',str(expmt_path+'/cc_measure'))
    os.makedirs(newpath)

for path_in in Path(expmt_path+'/org_measure').glob('*pre.csv'):
    path_parts=path_in.stem.split("_")
    fov=path_parts[4][:4]

    path_cell=Path(expmt_path+'/cell_measure')/f"BF-timelapse_03042026_eyrbow_glucose-2.0_{fov}_acdc_output_cpsam_d10.csv"

    list_in.append(path_in)
    list_cell.append(path_cell)

args = pd.DataFrame({
    "path_in":   list_in,
    "path_cell": list_cell
})

# %% pathing for debugging convenience
# pacdc=r'C:\Users\jglic\Downloads\11-5-25 myo1-mLemon\Run_2\BF-timelapse_acdc_output_cpsam-d10-ms8.csv'
# pmyo1=r'C:\Users\jglic\Downloads\11-5-25 myo1-mLemon\Run_2\eyrbow-yellowconfocal_zstack-avg-proj_myo1_fov2.csv'
# pbend=r'C:\Users\jglic\Downloads\11-5-25 myo1-mLemon\Run_2\BF-timelapse_manual-annotations.csv'
# palign=r'C:\Users\jglic\Downloads\11-5-25 myo1-mLemon\Run_2\manual_alignment_key.csv'

# %%
