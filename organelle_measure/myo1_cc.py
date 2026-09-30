# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import os
import sys
sys.path.append(r'C:\Users\jglic\Documents\School\WashU\Mukherji Lab\organelle-measure\organelle_measure')
from pathing_variables import expmt_path
from org_cc import cc_sort, find_all_transitions
# from tools import batch_apply
# %%
julcpath=r'C:\Users\jglic\Downloads\07172026 myo1 high res\cell_measure\BF-timelapse_acdc_output_07172026_cpsam-d25-ms20.csv'
julmaskpath=r'C:\Users\jglic\Downloads\07172026 myo1 high res\BF-timelapse_07172026_myo1-mlemon_glucose-2.0_fov1_segm-afftransf-expand.tif'
julmyo1path=r'C:\Users\jglic\Downloads\07172026 myo1 high res\confyellow-jg_07172026_myo1-mlemon_glucose-2.0_fov1_despeckle_r=4.csv'
julaapath=r'C:\Users\jglic\Downloads\07172026 myo1 high res\autoassigntable_r=4_norm.csv'

juncpath=r'C:\Users\jglic\Downloads\06102026 myo1 high res\BF-timelapse_acdc_output_cpsam-d25.csv'
junmaskpath=r'C:\Users\jglic\Downloads\06102026 myo1 high res\BF-timelapse_cpsam-d25_afftransf-TEST2.tif'
junmyo1path=r'C:\Users\jglic\Downloads\06102026 myo1 high res\SUM-confyellow-jg_06102026_ey2795-myo1_glucose-2.0_fov1_despeckle_r=4.csv'
junaapath=r'C:\Users\jglic\Downloads\06102026 myo1 high res\autoassigntable_r=4_norm.csv'

# cpath=r'C:\Users\jglic\Downloads\11052025 myo1\BF-timelapse_acdc_output_cpsam-d10-ms8.csv'
# maskpath=r'C:\Users\jglic\Downloads\11052025 myo1\BF-timelapse_cpsam-d10-ms8_afftransf_no-border.tif'
# myo1path=r'C:\Users\jglic\Downloads\11052025 myo1\eyrbow-yellowconfocal_zstack-avg-proj_myo1_fov2.csv'
# aapath=r'C:\Users\jglic\Downloads\11052025 myo1\autoassign-table_r=3.csv'
# %% Myo1 dynamics analysis
def ticks(myo1_df: pd.DataFrame, thresh: float=0.4):
    """
    Args: Dataframe containing extracted roi info to be analyzed.
    Outputs: Dataframe containing up and/or down tick frames for each unique roi. 
    """
    dfs=[]
    for roi_id in np.unique(myo1_df['ROI_ID'].values):
        vals=myo1_df.loc[myo1_df.ROI_ID==roi_id, 'sum'].values
        vals=vals/np.max(vals)
        yi=myo1_df.loc[myo1_df.ROI_ID==roi_id, 'y'].values[0]
        xi=myo1_df.loc[myo1_df.ROI_ID==roi_id, 'x'].values[0]
        frames=myo1_df.loc[myo1_df.ROI_ID==roi_id, 'frame'].values

        # on_thresh=min([50, 0.5*np.max(vals)])
        on_thresh=thresh
        on_indices=np.where(vals>on_thresh)[0]
        # off_thresh=max([50, 0.4*np.max(vals)])
        off_thresh=thresh
        off_indices=np.where(vals<off_thresh)[0]
        start_frame=None
        div_frame=None

        pre_on=[x for x in off_indices if x<on_indices[0]]
        if len(pre_on)>0:
            start_frame=frames[pre_on[-1]]+1  
        post_on=[x for x in off_indices if x>on_indices[-1]]
        if len(post_on)>0:
            div_frame=frames[post_on[0]]

        if len(pre_on)>0 or len(post_on)>0: #if there is an up/downtick, save info
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
def autoassign(myo1df_path: str, cellmask_path: str, cellacdc_path: str, camsamprate: int, consamprate: int, disk_radius: int=4):
    """
    Args: Dataframe containing detected myo1 ROI coords, cell segmentation mask file path, and cell-acdc output csv.
    Outputs: Dataframe containing ROI-Cell_ID assignments.
    """
    imgdf=pd.read_csv(myo1df_path)
    mask=io.imread(cellmask_path)
    acdcdf=pd.read_csv(cellacdc_path)
    rowsout=[]
    missed=0
    for roi_id in np.unique(imgdf['ROI_ID'].values):
        roi_rows=imgdf.loc[imgdf.ROI_ID==roi_id].copy()
        frames=roi_rows['frame'].values
        for frame_i in frames:
            real_time=frame_i*consamprate
            mask_frame=int(real_time/camsamprate)

            yi=roi_rows.loc[roi_rows.frame==frame_i,'y'].values[0]
            xi=roi_rows.loc[roi_rows.frame==frame_i,'x'].values[0]
            rr, cc = draw.disk((yi, xi), disk_radius, shape=(mask.shape[1],mask.shape[2])) #draw disk around the glob
            unique, counts = np.unique(mask[mask_frame,rr,cc][mask[mask_frame,rr,cc]>0], return_counts=True) #extract which cell masks fall within the disk and how much
            unique, counts = list(unique), list(counts)
            leave=False
            while leave==False: #loop here to only assign rois to mother cells
                if len(counts)>0:
                    max_index=np.argmax(counts)
                    cell_id=unique[max_index]
                    cell_rows=acdcdf.loc[acdcdf.Cell_ID==cell_id]
                    # print(cell_id, real_time, frame_i, mask_frame)
                    if cell_rows.empty==True or cell_rows.loc[cell_rows.frame_i==mask_frame].empty==True:
                        del unique[max_index]
                        del counts[max_index]
                    elif cell_rows.loc[cell_rows.frame_i==mask_frame, 'relationship'].values[0]=='mother': 
                        roi_rows.loc[roi_rows.frame==mask_frame, 'Cell_ID']=cell_id
                        leave=True
                    else:
                        del unique[max_index]
                        del counts[max_index]
                else:
                    leave=True
                    # print((yi,xi))
                    missed+=1
                    continue

        rowsout.append(roi_rows)
    print("Missed cell frames: "+str(missed))
    aadf=pd.concat(rowsout, ignore_index=True)
    output_name='autoassigntable_r='+str(disk_radius)+'.csv'
    output_path=Path(myo1df_path).parent/f"{output_name}"
    aadf.to_csv(output_path, index=False)
    return aadf


# %% cc size and myo1 overlay
def analyze_myo1size(cell_path: str, aa_path: str, thresh: float=0.45, cell_metric: str='cell_area_norm', myo1_metric: str='norm_int'):
    """
    Args: Paths to dataframes containing cell segm & annot outputs and myo1 intensity df with auto-matched cell ids.
    Outputs: Df containing myo1 and cell metrics. By default, deals with normalized values. For raw use: ['cell_area_pxl', 'sum']
    """
    celldf=pd.read_csv(cell_path)
    g1_df,s_df=cc_sort(cell_path)
    g1df,sdf=find_all_transitions(g1_df,s_df)
    # g1df,sdf=find_full_transitions(g1_df,s_df)
    g1_moms=g1df.loc[g1df.cc_group=='gmom'].copy()
    s_moms=sdf.loc[sdf.cc_group=='smom'].copy()
    g1_buds=g1df.loc[g1df.cc_group=='gbud'].copy()
    s_buds=sdf.loc[sdf.cc_group=='sbud'].copy()

    myo1df=pd.read_csv(aa_path)
    ticksdf=ticks(myo1df, thresh=thresh)
    
    for roi_id in np.unique(myo1df['ROI_ID'].values):
        if roi_id in ticksdf['ROI_ID'].values:
            roi_rows=myo1df.loc[myo1df.ROI_ID==roi_id]
            cell_id, count = stats.mode(roi_rows['Cell_ID'].values, nan_policy='omit')

            if cell_id in np.unique(g1_moms["Cell_ID"].values):
                cell_rows=g1_moms.loc[g1_moms.Cell_ID==cell_id].copy()

                myo1_start=ticksdf.loc[ticksdf.ROI_ID==roi_id, 'Start_Frame'].values
                if myo1_start!=None:
                    start_index=[(roi_rows['frame'].values[i]-myo1_start[0]) for i in range(len(roi_rows))]
                    myo1df.loc[myo1df.ROI_ID==roi_id, 'myo1_Start_Index']=start_index
                    g1_moms.loc[g1_moms.Cell_ID==cell_id, 'myo1_Start_Index']=[(cell_rows['frame_i'].values[i]-myo1_start[0]) for i in range(len(cell_rows))]
                    
                    annot_start=cell_rows.loc[cell_rows.Start_Index_1==0, 'frame_i']
                    if annot_start.empty==False:
                        if myo1_start[0]<annot_start.values[0]:
                            bud_id=cell_rows.loc[cell_rows.Start_Index_1==0, 'relative_ID'].values[0]
                            bud_rows=g1_buds.loc[g1_buds.Cell_ID==bud_id].copy()
                            g1_buds.loc[g1_buds.Cell_ID==bud_id, 'myo1_Start_Index']=[(bud_rows['frame_i'].values[i]-myo1_start[0]) for i in range(len(bud_rows))]
                
                myo1_div=ticksdf.loc[ticksdf.ROI_ID==roi_id, 'Div_Frame'].values
                if myo1_div!=None:
                    div_index=[(roi_rows['frame'].values[i]-myo1_div[0]) for i in range(len(roi_rows))]
                    myo1df.loc[myo1df.ROI_ID==roi_id, 'myo1_Div_Index']=div_index
                    g1_moms.loc[g1_moms.Cell_ID==cell_id, 'myo1_Div_Index']=[(cell_rows['frame_i'].values[i]-myo1_div[0]) for i in range(len(cell_rows))]
                    
                    annot_div=cell_rows.loc[cell_rows.Div_Index_1==0, 'frame_i']
                    if annot_div.empty==False:
                        if myo1_div[0]<annot_div.values[0]:
                            bud_id=cell_rows.loc[cell_rows.Div_Index_1==-1, 'relative_ID'].values[0]
                            bud_rows=g1_buds.loc[g1_buds.Cell_ID==bud_id].copy()
                            g1_buds.loc[g1_buds.Cell_ID==bud_id, 'myo1_Div_Index']=[(bud_rows['frame_i'].values[i]-myo1_div[0]) for i in range(len(bud_rows))]
                

            if cell_id in np.unique(s_moms["Cell_ID"].values):
                cell_rows=s_moms.loc[s_moms.Cell_ID==cell_id].copy()

                myo1_div=ticksdf.loc[ticksdf.ROI_ID==roi_id, 'Div_Frame'].values
                if myo1_div!=None:
                    div_index=[(roi_rows['frame'].values[i]-myo1_div[0]) for i in range(len(roi_rows))]
                    myo1df.loc[myo1df.ROI_ID==roi_id, 'myo1_Div_Index']=div_index
                    s_moms.loc[s_moms.Cell_ID==cell_id, 'myo1_Div_Index']=[(cell_rows['frame_i'].values[i]-myo1_div[0]) for i in range(len(cell_rows))]
                    
                    annot_div=cell_rows.loc[cell_rows.Div_Index_1==0, 'frame_i']
                    if annot_div.empty==False:
                        if myo1_div[0]<annot_div.values[0]:
                            bud_id=cell_rows.loc[cell_rows.Div_Index_1==-1, 'relative_ID'].values[0]
                            bud_rows=s_buds.loc[s_buds.Cell_ID==bud_id].copy()
                            s_buds.loc[s_buds.Cell_ID==bud_id, 'myo1_Div_Index']=[(bud_rows['frame_i'].values[i]-myo1_div[0]) for i in range(len(bud_rows))]
                
                myo1_start=ticksdf.loc[ticksdf.ROI_ID==roi_id, 'Start_Frame'].values
                if myo1_start!=None:
                    start_index=[(roi_rows['frame'].values[i]-myo1_start[0]) for i in range(len(roi_rows))]
                    myo1df.loc[myo1df.ROI_ID==roi_id, 'myo1_Start_Index']=start_index
                    s_moms.loc[s_moms.Cell_ID==cell_id, 'myo1_Start_Index']=[(cell_rows['frame_i'].values[i]-myo1_start[0]) for i in range(len(cell_rows))]
                    
                    annot_start=cell_rows.loc[cell_rows.Start_Index_1==0, 'frame_i']
                    if annot_start.empty==False:
                        if myo1_start[0]<annot_start.values[0]:
                            bud_id=cell_rows.loc[cell_rows.Start_Index_1==0, 'relative_ID'].values[0]
                            bud_rows=s_buds.loc[s_buds.Cell_ID==bud_id].copy()
                            s_buds.loc[s_buds.Cell_ID==bud_id, 'myo1_Start_Index']=[(bud_rows['frame_i'].values[i]-myo1_start[0]) for i in range(len(bud_rows))]

    dfs=[]

    myo1_minframe=-20
    myo1_maxframe=20                 
    for rel_frame in np.arange(myo1_minframe,myo1_maxframe,1):
        myo1_start_vals=myo1df.loc[myo1df.myo1_Start_Index==rel_frame, myo1_metric].values
        if len(myo1_start_vals)>0:
            myo1_start_metric=np.mean(myo1_start_vals)
            myo1_start_std=np.std(myo1_start_vals)
            metrics={
                "cat"           : "myo1-sta",
                "avg"           : [myo1_start_metric],
                "std"           : [myo1_start_std],
                "rel_frame"     : [rel_frame],
                "vals"          : [myo1_start_vals]
            }
            df=pd.DataFrame(metrics)
            dfs.append(df)

    for rel_frame in np.arange(myo1_minframe,myo1_maxframe,1):
        myo1_div_vals=myo1df.loc[myo1df.myo1_Div_Index==rel_frame, myo1_metric].values
        if len(myo1_div_vals)>0:
            myo1_div_metric=np.mean(myo1_div_vals)
            myo1_div_std=np.std(myo1_div_vals)
            metrics={
                "cat"           : "myo1-div",
                "avg"           : [myo1_div_metric],
                "std"           : [myo1_div_std],
                "rel_frame"     : [rel_frame],
                "vals"          : [myo1_div_vals]
            }
            df=pd.DataFrame(metrics)
            dfs.append(df)

    for rel_frame in np.arange(-50,50,1):
        cell_start_vals=np.concatenate([g1_moms.loc[g1_moms.myo1_Start_Index==rel_frame, cell_metric].values, s_moms.loc[s_moms.myo1_Start_Index==rel_frame, cell_metric].values])
        if len(cell_start_vals)>0:
            cell_start_metric=np.mean(cell_start_vals)
            cell_start_std=np.std(cell_start_vals)
            metrics={
                "cat"           : "mom-sta",
                "avg"           : [cell_start_metric],
                "std"           : [cell_start_std],
                "rel_frame"     : [rel_frame],
                "vals"          : [cell_start_vals]
            }
            df=pd.DataFrame(metrics)
            dfs.append(df)

        cell_div_vals=np.concatenate([s_moms.loc[s_moms.myo1_Div_Index==rel_frame, cell_metric].values, g1_moms.loc[g1_moms.myo1_Div_Index==rel_frame, cell_metric].values])
        if len(cell_div_vals)>0:
            cell_div_metric=np.mean(cell_div_vals)
            cell_div_std=np.std(cell_div_vals)
            metrics={
                "cat"           : "mom-div",
                "avg"           : [cell_div_metric],
                "std"           : [cell_div_std],
                "rel_frame"     : [rel_frame],
                "vals"          : [cell_div_vals]
            }
            df=pd.DataFrame(metrics)
            dfs.append(df)

        cell_div_vals=np.concatenate([s_buds.loc[s_buds.myo1_Div_Index==rel_frame, cell_metric].values, g1_buds.loc[g1_buds.myo1_Div_Index==rel_frame, cell_metric].values])
        if len(cell_div_vals)>0:
            cell_div_metric=np.mean(cell_div_vals)
            cell_div_std=np.std(cell_div_vals)
            metrics={
                "cat"           : "bud-div",
                "avg"           : [cell_div_metric],
                "std"           : [cell_div_std],
                "rel_frame"     : [rel_frame],
                "vals"          : [cell_div_vals]
            }
            df=pd.DataFrame(metrics)
            dfs.append(df)

    for rel_frame in np.arange(0,50,1):
        cell_start_vals=np.concatenate([g1_buds.loc[g1_buds.myo1_Start_Index==rel_frame, cell_metric].values, s_buds.loc[s_buds.myo1_Start_Index==rel_frame, cell_metric].values])
        if len(cell_start_vals)>0:
            cell_start_metric=np.mean(cell_start_vals)
            cell_start_std=np.std(cell_start_vals)
            metrics={
                "cat"           : "bud-sta",
                "avg"           : [cell_start_metric],
                "std"           : [cell_start_std],
                "rel_frame"     : [rel_frame],
                "vals"          : [cell_start_vals]
            }
            df=pd.DataFrame(metrics)
            dfs.append(df)

    metric_df=pd.concat(dfs, ignore_index=True)
    return metric_df

# junmdf=analyze_myo1size(juncpath, junaapath)
# julmdf=analyze_myo1size(julcpath, julaapath)
# %% Plot myo1 intensity profiles
def plot_myo1(mdf: pd.DataFrame, consamprate: int=3, ms: int=5):
    """
    Args: df containing processed cell and myo1 data; misc params for plotting
    Outputs: Plot displaying myo1 intensity as fnc of time.
    """
    fig,axes=plt.subplots(nrows=1,ncols=2,sharex=True)
    hline=0.45

    myo1sta=mdf.loc[mdf.cat=='myo1-sta']
    for time in myo1sta['rel_frame'].values: 
        y=myo1sta.loc[myo1sta.rel_frame==time, 'vals'].values[0]
        x=np.array([int(time)*consamprate for i in range(len(y))])
        axes[0].scatter(x, y,s=5,c='k',marker='o')
        axes[0].axhline(y=hline,xmin=0,xmax=1,ls='--', label='y= '+str(hline))
        # axes[0].set_ylim(-5,5)
    
    myo1div=mdf.loc[mdf.cat=='myo1-div']
    for time in myo1div['rel_frame'].values:      
        y=myo1div.loc[myo1div.rel_frame==time, 'vals'].values[0]
        x=np.array([int(time)*consamprate for i in range(len(y))])
        axes[1].scatter(x, y,s=5,c='k',marker='o')
        axes[1].axhline(y=hline,xmin=0,xmax=1,ls='--', label='y= '+str(hline))

    # axes[0].legend()
    # axes[1].legend()
    axes[0].set_title('Myo1 Intensity vs CC Position \n about Bud Neck Formation')
    axes[0].set_xlabel('Time Relative to Bud Neck Formation (min)')
    axes[0].set_ylabel('Normalized Myo1 Intensity')
    axes[1].set_title('Myo1 Intensity vs CC Position \n about Cell Division')
    axes[1].set_xlabel('Time Relative to Division (min)')
    plt.tight_layout()

# plot_myo1(junmdf, consamprate=5)
# plot_myo1(julmdf)
# %% Plot mother/daughter sizes relative to myo1 on/off
def plot_cc_size(mdf: pd.DataFrame, cam_samprate: int=2, con_samprate: int=5, ms: int=5):
    """
    Args: df containing processed cell and myo1 data; misc params for plotting
    Outputs: Plot displaying cell size overlay as function of myo1 dynamics.
    """
    fig,axes=plt.subplots(nrows=1,ncols=2,sharex=True)

    momsta=mdf.loc[mdf.cat=='mom-sta']
    for time in momsta['rel_frame'].values: 
        y=momsta.loc[momsta.rel_frame==time, 'vals'].values[0]
        x=np.array([int(time)*cam_samprate for i in range(len(y))])
        axes[0].scatter(x, y,s=5,c='k',marker='o')
        # axes[0].set_ylim(-5,5)
    
    momdiv=mdf.loc[mdf.cat=='mom-div']
    for time in momdiv['rel_frame'].values:      
        y=momdiv.loc[momdiv.rel_frame==time, 'vals'].values[0]
        x=np.array([int(time)*cam_samprate for i in range(len(y))])
        axes[1].scatter(x, y,s=5,c='k',marker='o')

    budsta=mdf.loc[mdf.cat=='bud-sta']
    for time in budsta['rel_frame'].values:     
        y=budsta.loc[budsta.rel_frame==time, 'vals'].values[0]
        x=np.array([int(time)*cam_samprate for i in range(len(y))])
        axes[0].scatter(x, y,s=5,c='r',marker='s')

    buddiv=mdf.loc[mdf.cat=='bud-div']
    for time in buddiv['rel_frame'].values:      
        y=buddiv.loc[buddiv.rel_frame==time, 'vals'].values[0]
        x=np.array([int(time)*cam_samprate for i in range(len(y))])
        axes[1].scatter(x, y,s=5,c='r',marker='s')

    axes[0].set_title('Cell Volume vs CC Position \n about Bud Neck Formation')
    axes[0].set_xlabel('Time Relative to Bud Neck Formation (min)')
    axes[0].set_ylabel('Estimated Cell Volume (um^3)')
    axes[1].set_title('Cell Volume vs CC Position \n about Cell Division')
    axes[1].set_xlabel('Time Relative to Division (min)')
    plt.tight_layout()

# %% Plot myo1 intensity - cell size overlay
def plot_myo1_overlay(mdf: pd.DataFrame, cam_samprate: int=1, con_samprate: int=3, bin_ss: bool=True, ms: int=5, cell: bool=True):
    """
    Args: df containing processed cell and myo1 data; misc params for plotting
    Outputs: Plot displaying binned myo1 intensity and cell size as function of myo1 dynamics.
    """
    fig,axes=plt.subplots(nrows=1,ncols=2,sharex=True)
    axes[0].errorbar(mdf.loc[mdf.cat=='myo1-sta', 'rel_frame'].values*con_samprate, mdf.loc[mdf.cat=='myo1-sta', 'avg'].values, yerr=mdf.loc[mdf.cat=='myo1-sta', 'std'].values,ls='none', c='m',marker='o', ms=ms)
    axes[1].errorbar(mdf.loc[mdf.cat=='myo1-div', 'rel_frame'].values*con_samprate, mdf.loc[mdf.cat=='myo1-div', 'avg'].values, yerr=mdf.loc[mdf.cat=='myo1-div', 'std'].values,ls='none', c='m',marker='o', ms=ms)
    if cell==True:
        axes[0].errorbar(mdf.loc[mdf.cat=='mom-sta', 'rel_frame'].values*cam_samprate, mdf.loc[mdf.cat=='mom-sta', 'avg'].values, yerr=mdf.loc[mdf.cat=='mom-sta', 'std'].values,ls='none', c='k',marker='o', ms=ms)
        axes[1].errorbar(mdf.loc[mdf.cat=='mom-div', 'rel_frame'].values*cam_samprate, mdf.loc[mdf.cat=='mom-div', 'avg'].values, yerr=mdf.loc[mdf.cat=='mom-div', 'std'].values,ls='none', c='k',marker='o', ms=ms)
        axes[0].errorbar(mdf.loc[mdf.cat=='bud-sta', 'rel_frame'].values*cam_samprate, mdf.loc[mdf.cat=='bud-sta', 'avg'].values, yerr=mdf.loc[mdf.cat=='bud-sta', 'std'].values,ls='none', c='r',marker='s', ms=ms)
        axes[1].errorbar(mdf.loc[mdf.cat=='bud-div', 'rel_frame'].values*cam_samprate, mdf.loc[mdf.cat=='bud-div', 'avg'].values, yerr=mdf.loc[mdf.cat=='bud-div', 'std'].values,ls='none', c='r',marker='s', ms=ms)

    if bin_ss==True:
        for i in range(len(mdf.loc[mdf.cat=='myo1-sta', 'rel_frame'].values)):
            axes[0].text(mdf.loc[mdf.cat=='myo1-sta', 'rel_frame'].values[i]*con_samprate,  mdf.loc[mdf.cat=='myo1-sta', 'avg'].values[i]+0.025, int(len(mdf.loc[mdf.cat=='myo1-sta', 'vals'].values[i])), ha="center", fontsize="small")
            # axes[1].text(myo1stats_div[:,0][i]*con_samprate, myo1stats_div[:,1][i]+0.025, int(myo1stats_div[:,3][i]), ha="center", fontsize="small")
        # for i in range(len(momstats_sta[:,0])):
        #     axes[0].text(momstats_sta[:,0][i]*cam_samprate, momstats_sta[:,1][i]+0.025, int(momstats_sta[:,3][i]), ha="center", fontsize="small")
        #     axes[1].text(momstats_div[:,0][i]*cam_
        # samprate, momstats_div[:,1][i]+0.025, int(momstats_div[:,3][i]), ha="center", fontsize="small")
        # for i in range(len(budstats_div[:,0])): 
        #     axes[0].text(budstats_sta[:,0][i]*cam_samprate, budstats_sta[:,1][i]+0.025, int(budstats_sta[:,3][i]), ha="center", fontsize="small")
        #     axes[1].text(budstats_div[:,0][i]*cam_samprate, budstats_div[:,1][i]+0.025, int(budstats_div[:,3][i]), ha="center", fontsize="small")

    # axes[0].set_title('Cell Metrics vs CC Position \n about Bud Neck Formation')
    axes[0].set_xlabel('Time Relative to Bud Neck Formation (min)')
    axes[0].set_ylabel('Average Estimated Cell Volume Normalized (vox)')
    axes[0].set_ylim(0,1)
    # axes[1].set_title('Cell Metrics vs CC Position \n about Division')
    axes[1].set_xlabel('Time Relative to Division (min)')
    axes[1].set_ylim(0,1)
    ax2y=axes[1].twinx()
    ax2y.set_ylabel('Average Normalized Myo1 Intensity', color='m')
    plt.tight_layout()

# %% Interpret custom cc annotations
def find_custom_starts(celldf_path: str):
    """
    Args: path to Cell ACDC dataframe that includes custom Start annotations.
    Outputs: Same dataframe excluding cells without custom annotations and includes index column.
    """
    celldf=pd.read_csv(celldf_path)

    for cell_id in np.unique(celldf['Cell_ID'].values):
        cell_rows=celldf.loc[celldf.Cell_ID==cell_id]
        if max(cell_rows['Start'].values)==1:
            start_frame_1=cell_rows.loc[cell_rows.Start==1, 'frame_i'].values[0]
            celldf.loc[celldf.Cell_ID==cell_id, 'custom_Start_Index_1']=[(cell_rows['frame_i'].values[j]-start_frame_1) for j in range(len(cell_rows))]
            if len(cell_rows.loc[cell_rows.Start==1])>1:
                start_frame_2=cell_rows.loc[cell_rows.Start==1, 'frame_i'].values[1]
                celldf.loc[celldf.Cell_ID==cell_id, 'custom_Start_Index_2']=[(cell_rows['frame_i'].values[j]-start_frame_2) for j in range(len(cell_rows))]
        elif max(cell_rows['Start'].values)==0:
            celldf=celldf.drop(cell_rows.index)
        else:
            print('What does this mean?') 

    return celldf

def find_customdiffs(cell_path, aamyo1_path: str, thresh: float=0.5, camsamprate: int=1, consamprate: int=3):
    """
    Args: Paths to cell segm output csv, myo1 img csv. Param for tick function and frame to min unit conversions.
    Outputs: Tuple of arrays containing the frame differences in cc progression between myo1 signal and CUSTOM manual annotations
    """
    myo1df=pd.read_csv(aamyo1_path)
    ticksdf=ticks(myo1df,thresh=thresh)
    all_cells=find_custom_starts(cell_path)
    myo1start, myo1div =[],[]
    dfs=[]

    for roi_id in np.unique(myo1df['ROI_ID'].values):
        if roi_id in ticksdf['ROI_ID'].values:
            roi_rows=myo1df.loc[myo1df.ROI_ID==roi_id]
            cell_id, count = stats.mode(roi_rows['Cell_ID'].values, nan_policy='omit')
            if cell_id in all_cells['Cell_ID'].values:
                cell_rows=all_cells.loc[all_cells.Cell_ID==cell_id].copy()
                # cc_group=cell_rows['cc_group'].values[0]
                custom_start_1=cell_rows.loc[cell_rows.Start==1, 'frame_i'].values[0]
                myo1_start=ticksdf.loc[ticksdf.ROI_ID==roi_id, 'Start_Frame'].values
                if myo1_start!=None:
                    # print('Cell_ID: '+str(cell_id)+', S1: '+str(camsamprate*custom_start_1)+', Myo1: '+str(consamprate*myo1_start[0]))
                    start_diff_1=(camsamprate*custom_start_1-consamprate*myo1_start[0])
                    # if np.abs(start_diff_1:=(camsamprate*custom_start_1-consamprate*myo1_start[0]))<100:
                    myo1start.append(start_diff_1)

                    metrics={
                                "Cell_ID"      : [cell_id], 
                                "ROI_ID"       : [roi_id],
                                # "cc_group"     : [cc_group[1:]],
                                "start_diff"   : [start_diff_1]
                            }
                    metric_df=pd.DataFrame(metrics)
                    dfs.append(metric_df)
                    # else:
                    #     custom_start_2=cell_rows.loc[cell_rows.Start==1, 'frame_i'].values[1]
                    #     start_diff_2=(camsamprate*custom_start_2-consamprate*myo1_start[0])
                    #     myo1start.append(start_diff_2)
                    #     metrics={
                    #                 "Cell_ID"      : [cell_id], 
                    #                 "ROI_ID"       : [roi_id],
                    #                 # "cc_group"     : [cc_group[1:]],
                    #                 "start_diff"   : [start_diff_2],
                    #             }
                    #     metric_df=pd.DataFrame(metrics)
                    #     dfs.append(metric_df)

    print(len(myo1start))
    diff_df=pd.concat(dfs, ignore_index=True)
    return diff_df

# diffdf=find_customdiffs(julcpath,julaapath)
# %% finding cc differences
import math
from scipy import stats
def find_ccdiffs(cell_path: str, aamyo1_path: str, thresh: float=0.4, camsamprate: int=1, consamprate: int=3):
    """
    Args: Paths to cell segm output csv, myo1 img csv. Param for tick function and frame to min unit conversions.
    Outputs: Tuple of arrays containing the frame differences in cc progression between myo1 signal and manual annotations
    """
    # nov0525_ignore=[490,409,472,192,323,304,434,468,468,482,524,453,465]

    myo1df=pd.read_csv(aamyo1_path)
    ticksdf=ticks(myo1df,thresh=thresh)
    g1_df,s_df=cc_sort(cell_path)
    g1df, sdf = find_all_transitions(g1_df, s_df) 
    all_cells=pd.concat([g1df, sdf], ignore_index=True)
    myo1start, myo1div =[],[]
    dfs=[]

    for roi_id in np.unique(myo1df['ROI_ID'].values):
        if roi_id in ticksdf['ROI_ID'].values:
            roi_rows=myo1df.loc[myo1df.ROI_ID==roi_id]
            cell_id, count = stats.mode(roi_rows['Cell_ID'].values, nan_policy='omit')
            if cell_id in all_cells['Cell_ID'].values:
                # if cell_id not in ignore_list:
                cell_rows=all_cells.loc[all_cells.Cell_ID==cell_id].copy()
                cc_group=cell_rows['cc_group'].values[0]

                annot_start_1=cell_rows.loc[cell_rows.Start_Index_1==0, 'frame_i']
                annot_start_2=cell_rows.loc[cell_rows.Start_Index_2==0, 'frame_i']

                if annot_start_1.empty==False:
                    myo1_start=ticksdf.loc[ticksdf.ROI_ID==roi_id, 'Start_Frame'].values
                    if myo1_start!=None:
                        yi=roi_rows.loc[roi_rows.frame==myo1_start[0],'y'].values[0]
                        xi=roi_rows.loc[roi_rows.frame==myo1_start[0],'x'].values[0]
                        if (start_diff_1:=(camsamprate*annot_start_1.values[0]-consamprate*myo1_start[0]))>0:
                            myo1start.append(start_diff_1)
                            metrics={
                                        "Cell_ID"      : [cell_id], 
                                        "ROI_ID"       : [roi_id],
                                        "cc_group"     : [cc_group[1:]],
                                        "start_diff"   : [start_diff_1],
                                    }
                            metric_df=pd.DataFrame(metrics)
                            dfs.append(metric_df)
                            # if start_diff_1<10:
                                # print("ROI: "+str(roi_id)+', '+str(start_diff_1)+', ('+str(xi)+','+str(yi)+')')

                        elif annot_start_2.empty==False:
                            if (start_diff_2:=(camsamprate*annot_start_2.values[0]-consamprate*myo1_start[0]))>0:
                                myo1start.append(start_diff_2)
                                metrics={
                                        "Cell_ID"      : [cell_id], 
                                        "ROI_ID"       : [roi_id],
                                        "cc_group"     : [cc_group[1:]],
                                        "start_diff"   : [start_diff_2],
                                        }
                                metric_df=pd.DataFrame(metrics)
                                dfs.append(metric_df)
                                # if start_diff_2<10:
                                    # print("ROI: "+str(roi_id)+', '+str(start_diff_2)+', ('+str(xi)+','+str(yi)+')')

                annot_div_1=cell_rows.loc[cell_rows.Div_Index_1==0, 'frame_i']
                annot_div_2=cell_rows.loc[cell_rows.Div_Index_2==0, 'frame_i']
                if annot_div_1.empty==False:
                    myo1_div=ticksdf.loc[ticksdf.ROI_ID==roi_id, 'Div_Frame'].values
                    if myo1_div!=None:
                        if (div_diff_1:=(camsamprate*annot_div_1.values[0]-consamprate*myo1_div[0]))>0:
                            myo1div.append(div_diff_1)
                            metrics={
                                        "Cell_ID"      : [cell_id], 
                                        "ROI_ID"       : [roi_id],
                                        "cc_group"     : [cc_group[1:]],
                                        "div_diff"     : [div_diff_1],
                                    }
                            metric_df=pd.DataFrame(metrics)
                            dfs.append(metric_df)
                        elif annot_div_2.empty==False:
                            if (div_diff_2:=(camsamprate*annot_div_2.values[0]-consamprate*myo1_div[0]))>0:
                                myo1div.append(div_diff_2)
                                metrics={
                                            "Cell_ID"      : [cell_id], 
                                            "ROI_ID"       : [roi_id],
                                            "cc_group"     : [cc_group[1:]],
                                            "div_diff"     : [div_diff_2],
                                        }
                                metric_df=pd.DataFrame(metrics)
                                dfs.append(metric_df)

    print(len(myo1start), len(myo1div))
    diff_df=pd.concat(dfs, ignore_index=True)
    return myo1start, myo1div
    # return diff_df

# %% Gaussian for curve fitting of waiting time distributions
def gaussian(x, u, amp, std):
    """
    Args: Standard parameters for gaussian function
    Outputs: Single value from evaluated fnc
    """
    return amp * np.exp( - (x - u)**2 / (2*std ** 2))

import scipy.special
from scipy.optimize import curve_fit
def poisson(x, lam):
    """
    Args: Standard parameters for poisson distribution function
    Outputs: Single value from evaluated fnc
    """
    return ((lam)**x * np.exp(-lam))/(scipy.special.factorial(x))
    
# %% Plot cc uncertainty histograms
from scipy import stats
def plot_distr(ms: np.ndarray, md: np.ndarray):
    """
    Args: Dataframe containing difference measurements between annotations and myo1 dynamics.
    Outputs: Two plots illustrating the distribution of said measurements.
    """
    # ms=diff_df.loc[diff_df.cc_group=='mom', 'start_diff'].dropna().values
    # bs=diff_df.loc[diff_df.cc_group=='bud', 'start_diff'].dropna().values
    # md=diff_df.loc[diff_df.cc_group=='mom', 'div_diff'].dropna().values
    # bd=diff_df.loc[diff_df.cc_group=='bud', 'div_diff'].dropna().values
    
    #Trim outlier data points; estimated bounds by manual empirical analysis.
    ms = [i for i in ms if 45>i]
    md = [i for i in md if 30>i]
    
    # Fit each distribution
    # plt.figure()
    bud_dist=stats.norm
    bud_res = stats.fit(bud_dist, ms, bounds=[(0,40), (1,10)])
    bud_mean=bud_dist.mean(*bud_res.params)
    bud_std=bud_dist.std(*bud_res.params)

    # plt.figure()
    # div_dist=stats.poisson
    div_dist=stats.skewnorm
    div_res = stats.fit(div_dist, md, bounds=[(0,5), (0,40), (1,10)])
    div_mean=div_dist.mean(*div_res.params)
    div_std=div_dist.std(*div_res.params)

    # custom_diffs=[val for val in diffdf['start_diff'].values if -20<val<30]
    # cus_dist=stats.skewnorm
    # cus_res = stats.fit(cus_dist, custom_diffs, bounds=[(0,5), (0,40), (1,10)])
    # cus_mean=cus_dist.mean(*cus_res.params)
    # cus_std=cus_dist.std(*cus_res.params)

    # Now want to overlay fit curves with raw disributions.
    plt.figure()

    n0, b0, _ = plt.hist(ms, bins='doane', alpha=0.75, edgecolor='black', density=True, color='y', 
    label='Buddings (n='+f'{len(ms)}) $\mu$='+str(np.round(bud_mean,2)) + ', $\sigma$='+str(np.round(bud_std,2)))
    budx=np.linspace(min(ms), max(ms), 100)
    plt.plot(budx, bud_dist.pdf(budx,*bud_res.params), 'k--', label='PDF fit')

    n1, b1, _ = plt.hist(md, bins='doane', alpha=0.75, edgecolor='black', density=True, color='b', 
    label='Divisions (n='+f'{len(md)}) $\mu$='+str(np.round(div_mean,2)) + ', $\sigma$='+str(np.round(div_std,2)))
    divx=np.linspace(min(md), max(md), 100)
    plt.plot(divx, div_dist.pdf(divx,*div_res.params), 'k--', label='PDF fit')

    # n2, b2, _ = plt.hist(custom_diffs, bins='doane', alpha=0.75, edgecolor='black', density=True, color='r', 
    # label='Custom Starts (n='+f'{len(custom_diffs)}) $\mu$='+str(np.round(np.mean(custom_diffs),2)))
    # cusx=np.linspace(min(custom_diffs), max(custom_diffs), 100)
    # plt.plot(cusx, cus_dist.pdf(cusx,*cus_res.params), 'k--', label='PDF fit')
    
    plt.title('PDFs and Fits of Waiting Time Distributions')
    plt.xlabel('Time (min)')
    plt.ylabel('Probability Density')
    plt.legend()
    plt.tight_layout()

    print('Budding mean = '+str(bud_dist.mean(*bud_res.params)))
    print('Budding std = '+str(bud_dist.std(*bud_res.params)))
    print('Division mean = '+str(div_dist.mean(*div_res.params)))
    print('Division std = '+str(div_dist.std(*div_res.params)))
    # print('Custom mean = '+str(cus_dist.mean(*cus_res.params)))
    # print('Custom std = '+str(cus_dist.std(*cus_res.params)))

    return bud_res.params, div_res.params


# %%
# junms, junmd = find_ccdiffs(juncpath, junaapath, camsamprate=2, consamprate=5)
# julms, julmd = find_ccdiffs(julcpath, julaapath, camsamprate=1, consamprate=3)
# ms, md = np.concatenate([junms, julms]), np.concatenate([julms, julmd])
# plot_distr(ms, md)