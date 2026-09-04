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
orgs=['ld','gl', 'vo']
cam_pxl_size=.108 #in microns
xbound=30

# %% 3c-ey2795 pathing
# acdc_paths=os.listdir(expmt_path+'/cell_measure')

#pathing fnc
#cpath is only path in cell_measure, 

# cpath=r'C:\Users\jglic\Downloads\07172026 myo1 high res\BF-timelapse_acdc_output_07172026_cpsam-d25-ms20.csv'
# maskpath=r'C:\Users\jglic\Downloads\07172026 myo1 high res\BF-timelapse_07172026_myo1-mlemon_glucose-2.0_fov1_segm-afftransf-expand.tif'
# myo1path=r'C:\Users\jglic\Downloads\07172026 myo1 high res\confyellow-jg_07172026_myo1-mlemon_glucose-2.0_fov1_despeckle_r=4.csv'
# aapath=r'C:\Users\jglic\Downloads\07172026 myo1 high res\autoassigntable_r=4.csv'

# cpath=r'C:\Users\jglic\Downloads\06102026 myo1 high res\BF-timelapse_acdc_output_cpsam-d25.csv'
# maskpath=r'C:\Users\jglic\Downloads\06102026 myo1 high res\BF-timelapse_cpsam-d25_afftransf-TEST2.tif'
# myo1path=r'C:\Users\jglic\Downloads\06102026 myo1 high res\smooth_SUM-confyellow-jg_06102026_ey2795-myo1_glucose-2.0_fov1_r=4.csv'
# aapath=r'C:\Users\jglic\Downloads\06102026 myo1 high res\autoassigntable_r=4.csv'

# cpath=r'C:\Users\jglic\Downloads\11052025 myo1\BF-timelapse_acdc_output_cpsam-d10-ms8.csv'
# maskpath=r'C:\Users\jglic\Downloads\11052025 myo1\BF-timelapse_cpsam-d10-ms8_afftransf_no-border.tif'
# myo1path=r'C:\Users\jglic\Downloads\11052025 myo1\eyrbow-yellowconfocal_zstack-avg-proj_myo1_fov2.csv'
# aapath=r'C:\Users\jglic\Downloads\11052025 myo1\autoassign-table_r=3.csv'

# %% Parse file name metadata
def parse_meta_organelle(name: str):
    """
    Args: Name is the stem of the ORGANELLE image file.
    Outptuts: Dictionary containing experiment metadata
    """
    #Unpack experiment metadata according to file naming convention.
    deconv, stk, time, organelle, date, strain, condition, field=name.split('_')
    # field,time=field_time.split('-')
    return {
        "organelle":  organelle,
        "date":       date,
        "strain":     strain,
        "condition":  condition,
        "field":      field,
        "time":       time[-1]
    }

def parse_meta_orgmeasure(name: str):
    """
    Args: Name is the stem of the ORGANELLE measure csv file.
    Outptuts: Dictionary containing experiment metadata
    """
    #Unpack experiment metadata according to file naming convention.
    organelle, date, strain, condition, field=name.split('_')
    # field,time=field_time.split('-')
    return {
        "organelle":  organelle,
        "date":       date,
        "strain":     strain,
        "condition":  condition,
        "field":      field,
    }
# %%
list_cell = [] #Initialize lists for cell and org csvs
list_in   = [] 

if not os.path.exists(newpath:=Path(expmt_path+'/cc_measure')):
    print('Creating folder ',str(expmt_path+'/cc_measure'))
    os.makedirs(newpath)

for path_in in Path(expmt_path+'/org_measure').glob('*.csv'):
    path_parts=path_in.stem.split("_")
    fov=path_parts[4][:4] 
    cell_parts=path_in.stem.split('-')
    cell_end="-".join(cell_parts[:3])[3:]
    path_acdc=Path(expmt_path+'/cell_measure')/f"BF-timelapse_{cell_end}_cpsam-d25-ms75.csv"
    # path_acdc=Path(expmt_path+'/cell_measure')/f"BF-timelapse_acdc_output_cpsam-d25-ms75.csv"

    list_in.append(path_in)
    list_cell.append(path_acdc)

# args = pd.DataFrame({
#     "path_in":   list_in,
#     "path_cell": list_cell
# })
# %% Normalization function
def norm(arr: np.ndarray):
    return arr/np.max(arr)
# %% Linear regression
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
def linreg(dataset: np.ndarray, n: int=1000, train_frac=0.5, plot_graph=False):
    """
    Args: 2D Array of data to perform regression upon, number of runs, fraction of data to use for training, optional graph.
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

# %% Calculate number of cc phases and details for a given cell
def extract_streaks(plist: list):
    """
    Args: List of cc phase annotations for a given cell.
    Outputs: Three lists together containing "complete" cc phase lengths, corresponding cc phase, and # of complete cc phases.
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

    return lens, phases, len(lens)
# %% find number of annotated cc transitions
def find_num_transitions(cell_rows: pd.DataFrame):
    """
    Args: Cell rows pertaining to a single Cell ID
    Outputs: The number of annotated cc transitions the cell undergoes
    """
    phase_list=cell_rows['cell_cycle_stage'].values
    lens,phases,num_phases=extract_streaks(phase_list)

    ind=0
    if cell_rows['frame_i'].values[0]!=0: #for G1 buds
        ind+=1
    
    return num_phases-1+ind
# %% Preprocess & clean up acdc_raw_output csv
def pp_acdc_output(acdc_path: str, remove_dead: bool=True, remove_excl: bool=True, bud_thresh: bool=True, norm_area_col: bool=True, approx_vol_col: bool=True) -> pd.DataFrame:
    """
    Args: Dataframe containing Cell-ACDC output metrics and analysis information.
    Outputs: Dataframe containing the same information sans dead or manually-excluded cells, maximum initial area bud threshold, addition of 
            normalized cell area column, and addition of vol approximation column. Multiple keyword boolean args to control individual operations. 
            All enabled by default.
    """
    df=pd.read_csv(acdc_path)
    df.columns=df.columns.str.strip() #remove leading and trailing spaces

    if remove_dead==True:
        dead_indices=df[df['is_cell_dead']==1].index
        if dead_indices.empty==True:
            df=df.drop(dead_indices)

    if remove_excl==True:
        excl_indices=df[df['is_cell_excluded']==1].index
        if excl_indices.empty==True:
            df=df.drop(excl_indices)

    initial_size_thresh=150 #pixels
    if bud_thresh==True:
        for cell_id in np.unique(df['Cell_ID'].values):
            if df.loc[df.Cell_ID==cell_id, 'frame_i'].values[0]>0 and df.loc[df.Cell_ID==cell_id, 'relationship'].values[0]=='bud' and df.loc[df.Cell_ID==cell_id, 'cell_area_pxl'].values[0]>initial_size_thresh:
                large_bud_indices=df.loc[df.Cell_ID==cell_id].index
                df=df.drop(large_bud_indices)

    for cell_id in np.unique(df['Cell_ID'].values):
        cell_rows=df.loc[df.Cell_ID==cell_id]
        if norm_area_col==True:
            norm_area=[(cell_rows['cell_area_pxl'].values[i]/np.max(cell_rows['cell_area_pxl'].values)) for i in range(len(cell_rows))]
            df.loc[df.Cell_ID==cell_id, 'cell_area_norm']=norm_area 
        if approx_vol_col==True:
            vol_estimate=[(df.loc[df.Cell_ID==cell_id, 'cell_area_pxl'].values[i]*df.loc[df.Cell_ID==cell_id, 'minor_axis_length'].values[i])*(cam_pxl_size**3) for i in range(len(cell_rows))]
            df.loc[df.Cell_ID==cell_id, 'approx_vol']=vol_estimate 
            vol_estimate_norm=vol_estimate/np.max(vol_estimate)
            df.loc[df.Cell_ID==cell_id, 'approx_vol_norm']=vol_estimate_norm
        # if add_transitions==True:
        #     num_transitions=find_num_transitions(cell_rows)
        #     df.loc[df.Cell_ID==cell_id, 'Transitions']= [num_transitions for i in range(len(cell_rows))]

    # output_name=Path(acdc_path).stem + '_pp.csv'
    # output_path=Path(acdc_path).parent/f"{output_name}"
    # df.to_csv(output_path, index=False)
    return df

# %% cc sort 
def cc_sort(acdc_path: str) -> (pd.DataFrame, pd.DataFrame):
    """
    Args: Path to raw acdc output 
    Outputs: Two cell metric dataframes (mothers AND daughters) distinguished by initial cell cycle position.
    """
    cell_df=pp_acdc_output(acdc_path)
    g1_cells=[]
    s_cells=[]
    nb_mothers=0
    nd_mothers=0
    for cell_id in np.unique(cell_df['Cell_ID'].values):
        cell_rows=cell_df[cell_df['Cell_ID']==cell_id].copy()
        if cell_rows['cell_cycle_stage'].values[0]=='G1':
            if len(cell_rows.loc[cell_rows.cell_cycle_stage=='S'])>0: #if G1 mother buds, keep
                cell_rows['cc_group']='gmom'
                g1_cells.append(cell_rows)

        elif cell_rows['relationship'].values[0]=='mother': 
            if len(cell_rows.loc[cell_rows.cell_cycle_stage=='G1'])>0: #if S mother divides, keep
                cell_rows['cc_group']='smom'
                s_cells.append(cell_rows)

        elif cell_rows['relationship'].values[0]=='bud': #same for buds
            if cell_rows['frame_i'].values[0]==0 and len(cell_rows.loc[cell_rows.cell_cycle_stage=='G1', 'frame_i'])>0:
                cell_rows['cc_group']='sbud'
                s_cells.append(cell_rows)
            elif cell_rows['frame_i'].values[0]>0:
                cell_rows['cc_group']='gbud'
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

# %% cc phase lengths for entire fov
def cc_lengths(cell_path: str):
    """
    Args: Path to cell segm output csv
    Outputs: Arrays of frame lengths of G1, S/G2/M phases and lists of complete cc-phase # distribution and composition.
    """
    g1_df,s_df=cc_sort(cell_path)
    g1_lens, s_lens, cc_lens =[],[],[]
    ind=0
    for cell_id in np.unique(g1_df.loc[g1_df.relationship=='mother','Cell_ID'].values): #for each G1 mother
        cell_rows=g1_df.loc[g1_df.Cell_ID==cell_id].copy()
        plist=cell_rows['cell_cycle_stage'].values
        lens,phases,num_phases = extract_streaks(plist)

        if len(lens)>2:
            del lens[0] #drop any streaks contacting the edges 
            del lens[-1]
            del phases[0]
            del phases[-1]
            for i in range(len(lens)):
                if phases[i]=='G1':
                    g1_lens.append(lens[i])
                else:
                    s_lens.append(lens[i])
        if len(lens)>1:
            for j in range(len(lens)):
                if len(lens)>(j+1): #now record lengths of full cell cycles
                    cc_lens.append(lens[j]+lens[j+1])
    
    for cell_id in np.unique(s_df.loc[s_df.relationship=='mother','Cell_ID'].values): #for each S mother
        cell_rows=s_df.loc[s_df.Cell_ID==cell_id].copy()
        plist=cell_rows['cell_cycle_stage'].values
        lens,phases,num_phases = extract_streaks(plist)
        if len(lens)>2:
            del lens[0] #drop any streaks contacting the edges 
            del lens[-1]
            del phases[0]
            del phases[-1]
            for i in range(len(lens)):
                if phases[i]=='G1':
                    g1_lens.append(lens[i])
                else:
                    s_lens.append(lens[i])
        if len(lens)>1:
            for j in range(len(lens)):
                if len(lens)>(j+1): #now record lengths of full cell cycles
                    cc_lens.append(lens[j]+lens[j+1])

    return g1_lens, s_lens, cc_lens

# g_lens,s_lens,lens_l=[],[],[]
# for path in acdc_paths:
#     gl,sl,lens_dist=cc_lengths(expmt_path+'/cell_measure'+'/'+path)
#     g_lens.append(gl)
#     s_lens.append(sl)
#     lens_l.append(lens_dist)
# g_lens=[element for innerList in g_lens for element in innerList]
# s_lens=[element for innerList in s_lens for element in innerList]
# lens_l=[element for innerList in lens_l for element in innerList]

# camsamprate=5
# plt.figure()
# plt.xlabel('# frames ('+str(camsamprate)+' min interval)')
# plt.title('Mother Cell CC Phase Annotation Lengths')
# bins=np.arange(0,max(cclen)+20,5)
# n0,b0,p0=plt.hist(slen,bins=bins, alpha=0.75, edgecolor='black', color='b', label='S/G2/M Phase n=('+f'{len(slen)})')
# n1,b1,p1=plt.hist(glen, bins=bins, alpha=0.75, edgecolor='black', color='y', label='G1 Phase n=('+f'{len(glen)})')
# n2,b2,p2=plt.hist(cclen, bins=bins, alpha=0.75, edgecolor='black', color='k', label='Total CC n=('+f'{len(cclen)})')
# plt.legend()

# %% full cc transition indices
def find_full_transitions(g1_df: pd.DataFrame, s_df: pd.DataFrame):
    """
    Args: cc_sort output -> dataframe of preprocessed cell info split in two according to starting cc phase
    Outputs: input dataframes modified to include cc-transitions indices and excludes cells that don't undergo 
                a fully-annotated cell cycle.
    """
    fullcc=0
    for cell_id in np.unique(g1_df['Cell_ID'].values): #for each G1 cell and their bud(s)   
        cell_rows=g1_df.loc[g1_df.Cell_ID==cell_id].copy()
        phase_list=cell_rows['cell_cycle_stage'].values
        lens,phases,num_phases = extract_streaks(phase_list)
        ind=0
        if cell_rows['relationship'].values[0]=='bud': #for G1 buds
            ind+=1
        num_transitions = num_phases-1+ind

        if num_transitions<3:
            g1_df=g1_df.drop(cell_rows.index)
        else:
            fullcc+=1
            checkpoint_frame=0
            for i in range(num_transitions):
                if i==0:
                    if cell_rows['relationship'].values[0]=='bud':
                        checkpoint_frame = cell_rows['frame_i'].values[0]
                    else: #for mothers
                        checkpoint_frame = lens[0]
                    column_name='Start_Index_1'

                else:
                    checkpoint_frame=checkpoint_frame+lens[i]
                
                    if i%2==0:
                        column_name='Start_Index_'+str(int(i/2 + 1))
                    else: #for S phase
                        column_name='Div_Index_'+str(int((i-1)/2 + 1))
                
                column_vals=[(cell_rows['frame_i'].values[j]-checkpoint_frame) for j in range(len(cell_rows))]
                g1_df.loc[g1_df.Cell_ID==cell_id, column_name]=column_vals

    for cell_id in np.unique(s_df['Cell_ID'].values): #for each S cell and their bud(s)   
        cell_rows=s_df.loc[s_df.Cell_ID==cell_id].copy()
        phase_list=cell_rows['cell_cycle_stage'].values
        lens,phases,num_phases = extract_streaks(phase_list)
        ind=0
        # if cell_rows['frame_i'].values[0]!=0: #for G1 buds
        #     ind+=1
        num_transitions = num_phases-1+ind

        if num_transitions<3:
            s_df=s_df.drop(cell_rows.index)
        else:
            fullcc+=1
            checkpoint_frame=0
            for i in range(num_transitions):
                checkpoint_frame=checkpoint_frame+lens[i]
                
                if i%2==0:
                    column_name='Div_Index_'+str(int(i/2 + 1))
                else: #for G1 phase
                    column_name='Start_Index_'+str(int((i-1)/2 + 1))
                
                column_vals=[(cell_rows['frame_i'].values[j]-checkpoint_frame) for j in range(len(cell_rows))]
                s_df.loc[s_df.Cell_ID==cell_id, column_name]=column_vals

    print('# full cell cycles = '+str(fullcc))
    return g1_df, s_df

# %% full cc organelle analysis
def analyze_fullcc(cell_dfpaths: list, org_dfpaths: list, save_csv: bool=False, cam_samprate: int=5, con_samprate: int=60):
    """
    Args: List of paths to raw acdc csvs and organelle measure csvs.
    Outputs: G1 and S/G2/M csvs containing geometric cell, cell cycle, and organelle information.
    """

    g1_dfs=[]
    s_dfs=[]

    for cell_path in cell_dfpaths:
        # celldf=pd.read_csv(cell_path)
        g1df_sort, sdf_sort = cc_sort(cell_path)
        g1_df, s_df = find_full_transitions(g1df_sort, sdf_sort)
        path_parts=cell_path.stem.split("_")
        date=path_parts[1]
        fov=path_parts[4]

        for cell_id in np.unique(g1_df['Cell_ID'].values):
            cell_rows=g1_df.loc[g1_df.Cell_ID==cell_id]
            start_1=cell_rows.loc[cell_rows.Start_Index_1==0, 'frame_i'].values[0]
            div_1=cell_rows.loc[cell_rows.Div_Index_1==0, 'frame_i'].values[0]
            start_2=cell_rows.loc[cell_rows.Start_Index_2==0, 'frame_i'].values[0]
            cc_length = (start_2 - start_1) * cam_samprate
            
            cell_metrics={
                "Cell_ID"           : [cell_id],
                "date"              : [date], 
                "fov"               : [fov],
                "relationship"      : [cell_rows['relationship'].values[0]],
                # "transitions"       : [cell_rows['transitions'].values[0]],
                "cc_length"         : [cc_length]
            }
            result=cell_metrics

            for org_path in org_dfpaths: 
                meta=parse_meta_orgmeasure(org_path.stem)
                org_label=meta['organelle']
                orgdf=pd.read_csv(org_path)
                org_rows=orgdf.loc[orgdf.idx_cell==cell_id]

                for frame in np.unique(org_rows['time'].values): #only do frames that have an org measurement
                    org_frame=org_rows.loc[org_rows.time==frame]
                    abs_time=org_frame['abs_time'].values[0]
                    frac_time = (abs_time - (start_1 * cam_samprate)) / (cc_length)
                    approx_vol=cell_rows.loc[cell_rows.time_minutes==int(abs_time), 'approx_vol'].values[0]
                    org_count=len(org_frame['idx-orga'].values)
                    org_vox=org_rows['volume-pixel'].values
                    if len(org_vox)>0:
                        org_vox_tot=np.sum(org_vox)
                    else:
                        org_vox_tot=org_vox[0]
                    org_vol_frac=(org_vox_tot*(cam_pxl_size**3)) / approx_vol

                    org_metrics={
                        "abs_time"         : [abs_time],
                        "frac_time"        : [frac_time],
                        "approx_vol"       : [approx_vol],
                        "cell_area_norm"   : [cell_rows.loc[cell_rows.time_minutes==int(abs_time), 'cell_area_norm'].values[0]],
                        "Relative_ID"      : [cell_rows.loc[cell_rows.time_minutes==int(abs_time), 'relative_ID'].values[0]],
                        org_label+"_vox"   : [org_vox_tot],
                        org_label+"_frac"  : [org_vol_frac],
                        org_label+"_count" : [org_count]
                            }
            
                    results = result | org_metrics
                    g1_dfs.append(pd.DataFrame(results))

        for cell_id in np.unique(s_df['Cell_ID'].values):
            cell_rows=s_df.loc[s_df.Cell_ID==cell_id]
            div_1=cell_rows.loc[cell_rows.Div_Index_1==0, 'frame_i'].values[0]
            start_1=cell_rows.loc[cell_rows.Start_Index_1==0, 'frame_i'].values[0]
            div_2=cell_rows.loc[cell_rows.Div_Index_2==0, 'frame_i'].values[0]
            cc_length = (div_2 - div_1) * cam_samprate
            
            cell_metrics={
                "Cell_ID"           : [cell_id],
                "date"              : [date], 
                "fov"               : [fov],
                "relationship"      : [cell_rows['relationship'].values[0]],
                # "transitions"       : [cell_rows['transitions'].values[0]],
                "cc_length"         : [cc_length]
            }
            result=cell_metrics

            for org_path in org_dfpaths: 
                meta=parse_meta_orgmeasure(org_path.stem)
                org_label=meta['organelle']
                orgdf=pd.read_csv(org_path)
                org_rows=orgdf.loc[orgdf.idx_cell==cell_id]

                for frame in np.unique(org_rows['time'].values): #only do frames that have an org measurement
                    org_frame=org_rows.loc[org_rows.time==frame]
                    abs_time=org_frame['abs_time'].values[0]
                    frac_time = (abs_time - (start_1 * cam_samprate)) / (cc_length)
                    approx_vol=cell_rows.loc[cell_rows.time_minutes==int(abs_time), 'approx_vol'].values[0]
                    org_count=len(org_frame['idx-orga'].values)
                    org_vox=org_rows['volume-pixel'].values
                    if len(org_vox)>0:
                        org_vox_tot=np.sum(org_vox)
                    else:
                        org_vox_tot=org_vox[0]
                    org_vol_frac=(org_vox_tot*(cam_pxl_size**3)) / approx_vol

                    org_metrics={
                        "abs_time"         : [abs_time],
                        "frac_time"        : [frac_time],
                        "approx_vol"       : [approx_vol],
                        "cell_area_norm"   : [cell_rows.loc[cell_rows.time_minutes==int(abs_time), 'cell_area_norm'].values[0]],
                        "Relative_ID"      : [cell_rows.loc[cell_rows.time_minutes==int(abs_time), 'relative_ID'].values[0]],
                        org_label+"_vox"   : [org_vox_tot],
                        org_label+"_frac"  : [org_vol_frac],
                        org_label+"_count" : [org_count]
                            }
            
                    results = result | org_metrics
                    s_dfs.append(pd.DataFrame(results))

        g1_org_df=pd.concat(g1_dfs, ignore_index=True)
        s_org_df=pd.concat(s_dfs, ignore_index=True)

    g1_org_df.sort_values(by=['frac_time', 'Cell_ID'], inplace=True, ignore_index=True)
    s_org_df.sort_values(by=['frac_time', 'Cell_ID'], inplace=True, ignore_index=True)

    if save_csv==True:
        g1_org_df.to_csv(Path(expmt_path+'\cc_measure')/f"{meta['date']}.csv",index=False)
        s_org_df.to_csv(Path(expmt_path+'\cc_measure')/f"{meta['date']}.csv",index=False)


    return g1_org_df, s_org_df
# %% plot full cc org dynamics
def extract_frame_metrics(cellrows: pd.DataFrame, org: str, metric: str):
    """
    Args: (frac, count, vox)
    Outputs: Two 1D arrays containing fractional time x-coords and chosen metric y-coords.
    """
    frac_times = cellrows['frac_time'].values
    metric_arr = cellrows[org+'_'+metric].values

    return frac_times, metric_arr

def plot_fullcc(g1_org_df: pd.DataFrame, s_org_df: pd.DataFrame, org_label: str, xl: float=-.5, xr: float=1.5, ms: int=5):
    """
    Args: Two dataframes for either G1 or S/G2/M; x-axis bound param, set markersize.
    Outputs: Plots of org metric vs cc positions for either df. Option to save plots.
    """ 
    metrics=['frac', 'count']
    fig,axes=plt.subplots(nrows=2,ncols=len(metrics))
    for cell in np.unique(g1_org_df['Cell_ID'].values):
        cell_rows=g1_org_df.loc[g1_org_df.Cell_ID==cell]
        if cell_rows['relationship'].values[0]=='mother':
            markerstyle='o'
            color='k'
        else:
            markerstyle='s'
            color='r'

        # for i in range(len(orgs)): 
        #     org_label=orgs[i]

        for j in range(len(metrics)):
            times, m_arr = extract_frame_metrics(cell_rows, org_label, metrics[j])
            axes[0, j].errorbar(times, m_arr, yerr=0, ls='none', c=color,marker=markerstyle,ms=ms)
            axes[0, j].set_xlabel('Normalized CC Position (Bud to Bud)')
        axes[0, 0].set_title('G1 '+org_label+' Vol Frac vs CC position')
        axes[0, 1].set_title('G1 '+org_label+' Count vs CC position')
        

    for cell in np.unique(s_org_df['Cell_ID'].values):
        cell_rows=s_org_df.loc[s_org_df.Cell_ID==cell]
        if cell_rows['relationship'].values[0]=='mother':
            markerstyle='o'
            color='k'
        else:
            markerstyle='s'
            color='r'

        # for i in range(len(orgs)): 
        #     org_label=orgs[i]

        for j in range(len(metrics)):
            times, m_arr = extract_frame_metrics(cell_rows, org_label, metrics[j])
            axes[1, j].errorbar(times, m_arr, yerr=0, ls='none', c=color,marker=markerstyle,ms=ms)
            axes[1, j].set_xlabel('Normalized CC Position (Div to Div)')
        axes[1, 0].set_title('S '+org_label+' Vol Frac vs CC position')
        axes[1, 1].set_title('S '+org_label+' Count vs CC position')


    fig.tight_layout()
    return

# %% ALL cc transition indices
def find_all_transitions(g1_df: pd.DataFrame, s_df: pd.DataFrame):
    """
    Args: cc_sort output Dataframes containing cells in G1 and S/G2/M phases respectively at frame 0.
    Outputs: Same dataframes with new columns for cc transition indexes and # cc transitions.
    """
    fullcc=0

    for cell_id in np.unique(g1_df['Cell_ID'].values): #for each G1 cell and their buds   
        cell_rows=g1_df.loc[g1_df.Cell_ID==cell_id].copy()
        plist=cell_rows['cell_cycle_stage'].values
        lens,phases,num_phases = extract_streaks(plist)
        ind=0
        if cell_rows['relationship'].values[0]=='bud': #for G1 buds
            ind+=1
        num_transitions = num_phases-1+ind

        if num_transitions==1: #those cells with ONE annotated cc transition
            if cell_rows['relationship'].values[0]=='bud':
                start_frame_1 = cell_rows['frame_i'].values[0]
            else: #for mothers
                start_frame_1 = lens[0]
            start_index_1=[(cell_rows['frame_i'].values[i]-start_frame_1) for i in range(len(cell_rows))]
            g1_df.loc[g1_df.Cell_ID==cell_id, 'Start_Index_1']=start_index_1
        
        elif num_transitions==2:
            if cell_rows['relationship'].values[0]=='bud':
                start_frame_1 = cell_rows['frame_i'].values[0]
                div_frame_1 = start_frame_1 + lens[0]
            else: #for mothers
                start_frame_1 = lens[0]
                div_frame_1 = start_frame_1 + lens[1]
            start_index_1=[(cell_rows['frame_i'].values[i]-start_frame_1) for i in range(len(cell_rows))]
            g1_df.loc[g1_df.Cell_ID==cell_id, 'Start_Index_1']=start_index_1
            div_index_1=[(cell_rows['frame_i'].values[i]-div_frame_1) for i in range(len(cell_rows))]
            g1_df.loc[g1_df.Cell_ID==cell_id, 'Div_Index_1']=div_index_1

        elif num_transitions>=3:
            fullcc+=1
            # four cc index columns, two for each annotated transition
            if cell_rows['relationship'].values[0]=='bud':
                start_frame_1 = cell_rows['frame_i'].values[0]
                div_frame_1 = start_frame_1 + lens[0]
                start_frame_2 = div_frame_1 + lens[1] 
            else: #for mothers
                start_frame_1 = lens[0]
                div_frame_1 = start_frame_1 + lens[1]
                start_frame_2 = div_frame_1 + lens[2] 

            start_index_1=[(cell_rows['frame_i'].values[i]-start_frame_1) for i in range(len(cell_rows))]
            g1_df.loc[g1_df.Cell_ID==cell_id, 'Start_Index_1']=start_index_1
            div_index_1=[(cell_rows['frame_i'].values[i]-div_frame_1) for i in range(len(cell_rows))]
            g1_df.loc[g1_df.Cell_ID==cell_id, 'Div_Index_1']=div_index_1
            start_index_2=[(cell_rows['frame_i'].values[i]-start_frame_2) for i in range(len(cell_rows))]
            g1_df.loc[g1_df.Cell_ID==cell_id, 'Start_Index_2']=start_index_2

    for cell_id in np.unique(s_df['Cell_ID'].values): #for each S/G2/M cell
        cell_rows=s_df.loc[s_df.Cell_ID==cell_id].copy()
        div_frame=cell_rows.loc[cell_rows.cell_cycle_stage=='G1', 'frame_i'].values[0]
        div_index=[(cell_rows['frame_i'].values[i]-div_frame) for i in range(len(cell_rows))]

        cell_rows.loc[cell_rows.Cell_ID==cell_id, 'Div_Index']=div_index
        plist=cell_rows['cell_cycle_stage'].values
        lens,phases,num_phases = extract_streaks(plist)
        ind=0
        num_transitions = num_phases-1+ind

        if num_transitions==1:
            div_frame_1 = lens[0]
            div_index_1=[(cell_rows['frame_i'].values[i]-div_frame_1) for i in range(len(cell_rows))]
            s_df.loc[s_df.Cell_ID==cell_id, 'Div_Index_1']=div_index_1

        elif num_transitions==2:
            div_frame_1 = lens[0]
            start_frame_1 = div_frame_1 + lens[1]
            div_index_1=[(cell_rows['frame_i'].values[i]-div_frame_1) for i in range(len(cell_rows))]
            s_df.loc[s_df.Cell_ID==cell_id, 'Div_Index_1']=div_index_1
            start_index_1=[(cell_rows['frame_i'].values[i]-start_frame_1) for i in range(len(cell_rows))]
            s_df.loc[s_df.Cell_ID==cell_id, 'Start_Index_1']=start_index_1
        
        elif num_transitions>=3:
            fullcc+=1
            div_frame_1 = lens[0]
            start_frame_1 = div_frame_1 + lens[1]
            div_frame_2 = start_frame_1 + lens[2] 
            div_index_1=[(cell_rows['frame_i'].values[i]-div_frame_1) for i in range(len(cell_rows))]
            s_df.loc[s_df.Cell_ID==cell_id, 'Div_Index_1']=div_index_1
            start_index_1=[(cell_rows['frame_i'].values[i]-start_frame_1) for i in range(len(cell_rows))]
            s_df.loc[s_df.Cell_ID==cell_id, 'Start_Index_1']=start_index_1
            div_index_2=[(cell_rows['frame_i'].values[i]-div_frame_2) for i in range(len(cell_rows))]
            s_df.loc[s_df.Cell_ID==cell_id, 'Div_Index_2']=div_index_2

    print('# full cell cycles = '+str(fullcc))
    return g1_df, s_df

#%% ALL cc organelle analysis
def cc_org_analysis(cell_dfpaths: list, pre_org_dfpaths: list, vol_frac=False) -> (pd.DataFrame,pd.DataFrame):
    """
    Args: Lists of paths to csv files containing extracted cell and organelle information respectively.
    Outputs: Two dataframes containing selected cell and organelle metrics, and relative cc transition frames (G1/S 
            for G1 df, M/G1 for S/G2/M df). Optional plots
    """
    g1_dfs, s_dfs, s_b = [],[],[]
    for i in range(len(cell_dfpaths)):
        g1_df,s_df=cc_sort(cell_dfpaths[i]) #prelim sort and process
        
        g1_moms, s_moms, g1_buds, s_buds = find_all_transitions(g1_df, s_df) #find annotated cc transitions
        path_parts=cell_dfpaths[i].stem.split("_")
        date=path_parts[1]
        fov=path_parts[4]
        fov_org_csvs=[p for p in pre_org_dfpaths if p.stem.split("_")[1]==date and p.stem.split("_")[4][:4]==fov] #parse out org csvs for this fov

        for cell_id in np.unique(g1_moms['Cell_ID'].values): #for each G1 mom
            transitions=g1_moms.loc[g1_moms.Cell_ID==cell_id, 'Transitions'].values[0]
            pre_start_offset=g1_moms.loc[g1_moms.Cell_ID==cell_id, 'Start_Index'].values[0]
            if transitions==1:
                post_label="post_start_offset"
                post_offset=g1_moms.loc[g1_moms.Cell_ID==cell_id, 'Start_Index'].values[-1]
            elif transitions==2:
                post_label="post_div_offset"
                post_offset=g1_moms.loc[g1_moms.Cell_ID==cell_id, 'Div_Index'].values[-1]
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
                    pre_cell_vol=approx_vol(g1_moms, cell_id, 0)
                    # pre_cell_vol=g1_moms.loc[g1_moms.Cell_ID==cell_id, 'cell_vol_vox'].values[0]
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
                    post_cell_vol=approx_vol(g1_moms, cell_id, -1)
                    # post_cell_vol=g1_moms.loc[g1_moms.Cell_ID==cell_id, 'cell_vol_vox'].values[-1]
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
            g1_dfs.append(pd.DataFrame(result))

        for cell_id in np.unique(g1_buds['Cell_ID'].values): #for each G1 bud
            transitions=g1_buds.loc[g1_buds.Cell_ID==cell_id, 'Transitions'].values[0]
            pre_start_offset=g1_buds.loc[g1_buds.Cell_ID==cell_id, 'Start_Index'].values[0]
            if transitions==1:
                post_label="post_start_offset"
                post_offset=g1_buds.loc[g1_buds.Cell_ID==cell_id, 'Start_Index'].values[-1]
            elif transitions==2:
                post_label="post_div_offset"
                post_offset=g1_buds.loc[g1_buds.Cell_ID==cell_id, 'Div_Index'].values[-1]
            else:
                continue
            cell_metrics={
                "idx_cell"          : [cell_id],
                "date"              : [date], 
                "fov"               : [fov],
                "relationship"      : ["bud"],
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
                    pre_cell_vol=approx_vol(g1_buds, cell_id, 0)    
                    # pre_cell_vol=g1_buds.loc[g1_buds.Cell_ID==cell_id, 'cell_vol_vox'].values[0]
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
                    post_cell_vol=approx_vol(g1_buds, cell_id, -1) 
                    # post_cell_vol=g1_buds.loc[g1_buds.Cell_ID==cell_id, 'cell_vol_vox'].values[-1]
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
            g1_dfs.append(pd.DataFrame(result))

        for cell_id in np.unique(s_moms['Cell_ID'].values): #for each S mom
            pre_div_offset=s_moms.loc[s_moms.Cell_ID==cell_id, 'Div_Index'].values[0]
            transitions=s_moms.loc[s_moms.Cell_ID==cell_id, 'Transitions'].values[0]
            if transitions==1:
                post_label="post_div_offset"
                post_offset=s_moms.loc[s_moms.Cell_ID==cell_id, 'Div_Index'].values[-1]
            elif transitions==2:
                post_label="post_start_offset"
                post_offset=s_moms.loc[s_moms.Cell_ID==cell_id, 'Start_Index'].values[-1]
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
                    pre_cell_vol=approx_vol(s_moms, cell_id, 0) 
                    # pre_cell_vol=s_moms.loc[s_moms.Cell_ID==cell_id, 'cell_vol_vox'].values[0]
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
                    post_cell_vol=approx_vol(s_moms, cell_id, -1) 
                    # post_cell_vol=s_moms.loc[s_moms.Cell_ID==cell_id, 'cell_vol_vox'].values[0]
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

        for cell_id in np.unique(s_buds['Cell_ID'].values): #for each S bud
            pre_div_offset=s_buds.loc[s_buds.Cell_ID==cell_id, 'Div_Index'].values[0]
            transitions=s_buds.loc[s_buds.Cell_ID==cell_id, 'Transitions'].values[0]
            if transitions==1:
                post_label="post_div_offset"
                post_offset=s_buds.loc[s_buds.Cell_ID==cell_id, 'Div_Index'].values[-1]
            elif transitions==2:
                post_label="post_start_offset"
                post_offset=s_buds.loc[s_buds.Cell_ID==cell_id, 'Start_Index'].values[-1]
            else:
                continue
            cell_metrics={
                "idx_cell"          : [cell_id],
                "date"              : [date], 
                "fov"               : [fov],
                "relationship"      : ["bud"],
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
                    pre_cell_vol=approx_vol(s_buds, cell_id, 0) 
                    # pre_cell_vol=s_buds.loc[s_buds.Cell_ID==cell_id, 'cell_vol_vox'].values[0]
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
                    post_cell_vol=approx_vol(s_buds, cell_id, -1) 
                    # post_cell_vol=s_buds.loc[s_buds.Cell_ID==cell_id, 'cell_vol_vox'].values[-1]
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

    g1_org_df=pd.concat(g1_dfs, ignore_index=True)
    s_org_df=pd.concat(s_dfs, ignore_index=True)
    
    return g1_org_df, s_org_df
# %% Plot ALL cc org dynamics
def plot_ccorg(g1_org_df: pd.DataFrame, s_org_df: pd.DataFrame, save_fig=False, xbound: int=30,ms=3):
    """
    Args: Two dataframes for either G1 or S/G2/M. Bool to save figure, x-axis bound param, set markersize.
    Outputs: Plots of org metric vs cc positions for either df. Option to save plots.
    """ 
    g1_moms=g1_org_df.loc[g1_org_df.relationship=='mother']
    g1_buds=g1_org_df.loc[g1_org_df.relationship=='bud']
    s_moms=s_org_df.loc[s_org_df.relationship=='mother']
    s_buds=s_org_df.loc[s_org_df.relationship=='bud']
    #Plot desired organelle metric as function of cc position wrt reference cc checkpoint for both groups.
    fig,axes=plt.subplots(nrows=len(orgs),ncols=2)
    for i in range(len(orgs)): #should i keep the four offset columns or distinguish using transitions count in here?
        org_label=orgs[i]
        #all the pre's
        axes[i, 0].errorbar(g1_moms['pre_start_offset'].values, g1_moms[org_label+"_vol_pre"].values, yerr=org_err[org_label]*g1_moms[org_label+"_vol_pre"].values, ls='none', c='m',marker='o',ms=ms)
        axes[i, 1].errorbar(s_moms['pre_div_offset'].values, s_moms[org_label+"_vol_pre"].values, yerr=org_err[org_label]*s_moms[org_label+"_vol_pre"].values,ls='none', c='m',marker='o',ms=ms) 
        axes[i, 1].errorbar(s_buds['pre_div_offset'].values, s_buds[org_label+"_vol_pre"].values, yerr=org_err[org_label]*s_buds[org_label+"_vol_pre"].values,ls='none', c='r',marker='s',ms=ms) 
        #then then the 1-transitions
        axes[i, 0].errorbar(g1_moms.loc[g1_moms.transitions==1,'post_start_offset'].values,g1_moms.loc[g1_moms.transitions==1,org_label+"_vol_post"].values, yerr=org_err[org_label]*g1_moms.loc[g1_org_df.transitions==1,org_label+"_vol_post"].values,ls='none', c='m',marker='o',ms=ms)
        axes[i, 1].errorbar(s_moms.loc[s_moms.transitions==1,'post_div_offset'].values, s_moms.loc[s_moms.transitions==1,org_label+"_vol_post"].values, yerr=org_err[org_label]*s_moms.loc[s_moms.transitions==1,org_label+"_vol_post"].values,ls='none', c='m',marker='o',ms=ms)
        axes[i, 0].errorbar(g1_buds.loc[g1_buds.transitions==1,'post_start_offset'].values,g1_buds.loc[g1_buds.transitions==1,org_label+"_vol_post"].values, yerr=org_err[org_label]*g1_buds.loc[g1_org_df.transitions==1,org_label+"_vol_post"].values,ls='none', c='r',marker='s',ms=ms)
        axes[i, 1].errorbar(s_buds.loc[s_buds.transitions==1,'post_div_offset'].values, s_buds.loc[s_buds.transitions==1,org_label+"_vol_post"].values, yerr=org_err[org_label]*s_buds.loc[s_buds.transitions==1,org_label+"_vol_post"].values,ls='none', c='r',marker='s',ms=ms)
        
        #then the 2-transitions
        axes[i, 1].errorbar(g1_moms.loc[g1_moms.transitions==2,'post_div_offset'].values,g1_moms.loc[g1_moms.transitions==2,org_label+"_vol_post"].values, yerr=org_err[org_label]*g1_moms.loc[g1_moms.transitions==2,org_label+"_vol_post"].values,ls='none', c='m',marker='o',ms=ms)
        axes[i, 1].errorbar(g1_buds.loc[g1_buds.transitions==2,'post_div_offset'].values,g1_buds.loc[g1_buds.transitions==2,org_label+"_vol_post"].values, yerr=org_err[org_label]*g1_buds.loc[g1_buds.transitions==2,org_label+"_vol_post"].values,ls='none', c='r',marker='s',ms=ms)
        axes[i, 0].errorbar(s_moms.loc[s_moms.transitions==2,'post_start_offset'].values, s_moms.loc[s_moms.transitions==2,org_label+"_vol_post"].values, yerr=org_err[org_label]*s_moms.loc[s_moms.transitions==2,org_label+"_vol_post"].values,ls='none', c='m',marker='o',ms=ms)
        axes[i, 0].errorbar(s_buds.loc[s_buds.transitions==2,'post_start_offset'].values, s_buds.loc[s_buds.transitions==2,org_label+"_vol_post"].values, yerr=org_err[org_label]*s_buds.loc[s_buds.transitions==2,org_label+"_vol_post"].values,ls='none', c='r',marker='s',ms=ms)

        
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
def binned_org(g1_org_df: pd.DataFrame, s_org_df: pd.DataFrame,plot_graph: bool=True,ms=3):
    """
    Args: Dataframes of extracted cell and organelle metrics for each cc group. Option to plot.
    Outputs: Graphs displaying org vs cc trends, binned along the cc position axis.
    """
    all_cells=pd.concat([g1_org_df,s_org_df], ignore_index=True)
    moms=all_cells.loc[all_cells.relationship=="mother"]
    buds=all_cells.loc[all_cells.relationship=="bud"]
    if plot_graph==True:
        plt_ind=0
        fig,axes=plt.subplots(nrows=len(orgs),ncols=2,sharex=True)

    scores={}
    for org in orgs:
        mom_start_out=[]
        for time in np.unique(moms['pre_start_offset'].values):
            bin_vals=moms.loc[moms.pre_start_offset==time,org+"_vol_pre"].values
            bin_vals=[entry for entry in bin_vals if ~np.isnan(entry)]
            if len(bin_vals)>0:
                bin_avg=np.mean(bin_vals)
                bin_std=np.std(bin_vals)
                mom_start_out.append([time,bin_avg,bin_std])
        for time in np.unique(moms['post_start_offset'].values):
            bin_vals=moms.loc[moms.post_start_offset==time,org+"_vol_post"].values
            bin_vals=[entry for entry in bin_vals if ~np.isnan(entry)]
            if len(bin_vals)>0:
                bin_avg=np.mean(bin_vals)
                bin_std=np.std(bin_vals)
                mom_start_out.append([time,bin_avg,bin_std])

        mom_div_out=[]
        for time in np.unique(moms['pre_div_offset'].values):
            bin_vals=moms.loc[moms.pre_div_offset==time,org+"_vol_pre"].values
            bin_vals=[entry for entry in bin_vals if ~np.isnan(entry)]
            if len(bin_vals)>0:
                bin_avg=np.mean(bin_vals)
                bin_std=np.std(bin_vals)
                mom_div_out.append([time,bin_avg,bin_std])
        for time in np.unique(moms['post_div_offset'].values):
            bin_vals=moms.loc[moms.post_div_offset==time,org+"_vol_post"].values
            bin_vals=[entry for entry in bin_vals if ~np.isnan(entry)]
            if len(bin_vals)>0:
                bin_avg=np.mean(bin_vals)
                bin_std=np.std(bin_vals)
                mom_div_out.append([time,bin_avg,bin_std])
        
        bud_start_out=[]
        for time in np.unique(buds['pre_start_offset'].values):
            bin_vals=buds.loc[buds.pre_start_offset==time,org+"_vol_pre"].values
            bin_vals=[entry for entry in bin_vals if ~np.isnan(entry)]
            if len(bin_vals)>0:
                bin_avg=np.mean(bin_vals)
                bin_std=np.std(bin_vals)
                bud_start_out.append([time,bin_avg,bin_std])
        for time in np.unique(buds['post_start_offset'].values):
            bin_vals=buds.loc[buds.post_start_offset==time,org+"_vol_post"].values
            bin_vals=[entry for entry in bin_vals if ~np.isnan(entry)]
            if len(bin_vals)>0:
                bin_avg=np.mean(bin_vals)
                bin_std=np.std(bin_vals)
                bud_start_out.append([time,bin_avg,bin_std])

        bud_div_out=[]
        for time in np.unique(buds['pre_div_offset'].values):
            bin_vals=buds.loc[buds.pre_div_offset==time,org+"_vol_pre"].values
            bin_vals=[entry for entry in bin_vals if ~np.isnan(entry)]
            if len(bin_vals)>0:
                bin_avg=np.mean(bin_vals)
                bin_std=np.std(bin_vals)
                bud_div_out.append([time,bin_avg,bin_std])
        for time in np.unique(buds['post_div_offset'].values):
            bin_vals=buds.loc[buds.post_div_offset==time,org+"_vol_post"].values
            bin_vals=[entry for entry in bin_vals if ~np.isnan(entry)]
            if len(bin_vals)>0:
                bin_avg=np.mean(bin_vals)
                bin_std=np.std(bin_vals)
                bud_div_out.append([time,bin_avg,bin_std])

        m_start_out=np.array(mom_start_out)
        m_div_out=np.array(mom_div_out)
        b_start_out=np.array(bud_start_out)
        b_div_out=np.array(bud_div_out)
        # g1_out,s_out=np.array(g1_out),np.array(s_out)

        #regression analysis here ###############
        # g1_scores=linreg(g1_out[:,0:2])
        # s_scores=linreg(s_out[:,0:2])
        # score={
        #     org:[g1_scores,s_scores]
        # }
        # scores=scores|score

        if plot_graph==True:
            axes[plt_ind,0].errorbar(m_start_out[:,0],m_start_out[:,1],yerr=m_start_out[:,2],ls='none', c='m',marker='o', ms=ms)
            axes[plt_ind,1].errorbar(m_div_out[:,0],m_div_out[:,1],yerr=m_div_out[:,2],ls='none', c='m',marker='o', ms=ms)
            axes[plt_ind,0].errorbar(b_start_out[:,0],b_start_out[:,1],yerr=b_start_out[:,2],ls='none', c='r',marker='s', ms=ms)
            axes[plt_ind,1].errorbar(b_div_out[:,0],b_div_out[:,1],yerr=b_div_out[:,2],ls='none', c='r',marker='s', ms=ms)
            
            axes[plt_ind,0].set_title('G1 '+org+' vs CC position')
            # axes[plt_ind,0].set_xlabel('Frames relative to Start (5 min interval)')
            # axes[plt_ind,0].set_ylabel('Binned'+org+' volume')
            axes[plt_ind,0].set_xlim(-xbound,xbound)
            # axes[plt_ind,0].set_ylim(-.05,1)
            axes[plt_ind,1].set_title('S/G2/M '+org+' vs CC position')
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
    return

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
    g1_moms,s_moms,s_buds=find_all_transitions(g1_df,s_df)
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

# %% Myo1 dynamics analysis
def ticks(myo1_df: pd.DataFrame, thresh: float=0.5):
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

        on_indices=np.where(vals>thresh)[0]
        off_indices=np.where(vals<thresh)[0]
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
                    print(cell_id, real_time, frame_i, mask_frame)
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
def analyze_myo1size(cell_path: str, aa_path: str, thresh: float=0.5, cell_metric: str='cell_area_norm', myo1_metric: str='norm_int'):
    """
    Args: Paths to dataframes containing cell segm & annot outputs and myo1 intensity df with auto-matched cell ids.
    Outputs: Df containing myo1 and cell metrics. By default, deals with normalized values. For raw use: ['cell_area_pxl', 'sum']
    """

    #so like, find ticks from myo1df, use detected start/div frames to choose alignment point, extract corresponding cell area frames from celldf, plot aligned
    # then, also somehow average each profile by frame -- on that note, need to plot both sampling rates on same graph -> need to take along x-values.

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

# %% Plot myo1 intensity profiles
def plot_myo1(mdf: pd.DataFrame, con_samprate: int=5, ms: int=5):
    """
    Args: df containing processed cell and myo1 data; misc params for plotting
    Outputs: Plot displaying myo1 intensity as fnc of time.
    """
    fig,axes=plt.subplots(nrows=1,ncols=2,sharex=True)

    myo1sta=mdf.loc[mdf.cat=='myo1-sta']
    for time in myo1sta['rel_frame'].values: 
        y=myo1sta.loc[myo1sta.rel_frame==time, 'vals'].values[0]
        x=np.array([int(time)*con_samprate for i in range(len(y))])
        axes[0].scatter(x, y,s=5,c='k',marker='o')
        # axes[0].set_ylim(-5,5)
    
    myo1div=mdf.loc[mdf.cat=='myo1-div']
    for time in myo1div['rel_frame'].values:      
        y=myo1div.loc[myo1div.rel_frame==time, 'vals'].values[0]
        x=np.array([int(time)*con_samprate for i in range(len(y))])
        axes[1].scatter(x, y,s=5,c='k',marker='o')

    axes[0].set_title('Myo1 Intensity vs CC Position \n about Bud Neck Formation')
    axes[0].set_xlabel('Time Relative to Bud Neck Formation (min)')
    axes[0].set_ylabel('Estimated Cell Volume (um^3)')
    axes[1].set_title('Myo1 Intensity vs CC Position \n about Cell Division')
    axes[1].set_xlabel('Time Relative to Division (min)')
    plt.tight_layout()

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
def plot_myo1_overlay(mdf: pd.DataFrame, cam_samprate: int=2, con_samprate: int=5, bin_ss: bool=True, ms: int=5):
    """
    Args: df containing processed cell and myo1 data; misc params for plotting
    Outputs: Plot displaying binned myo1 intensity and cell size as function of myo1 dynamics.
    """
    fig,axes=plt.subplots(nrows=1,ncols=2,sharex=True)
    axes[0].errorbar(mdf.loc[mdf.cat=='myo1-sta', 'rel_frame'].values*con_samprate, mdf.loc[mdf.cat=='myo1-sta', 'avg'].values, yerr=mdf.loc[mdf.cat=='myo1-sta', 'std'].values,ls='none', c='m',marker='o', ms=ms)
    axes[1].errorbar(mdf.loc[mdf.cat=='myo1-div', 'rel_frame'].values*con_samprate, mdf.loc[mdf.cat=='myo1-div', 'avg'].values, yerr=mdf.loc[mdf.cat=='myo1-div', 'std'].values,ls='none', c='m',marker='o', ms=ms)
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

    axes[0].set_title('Cell Volume vs CC Position \n about Bud Neck Formation')
    axes[0].set_xlabel('Time Relative to Bud Neck Formation (min)')
    axes[0].set_ylabel('Average Estimated Cell Volume Normalized (vox)')
    axes[1].set_title('Cell Volume vs CC Position \n about Division')
    axes[1].set_xlabel('Time Relative to Division (min)')
    ax2y=axes[1].twinx()
    ax2y.set_ylabel('Average Normalized Myo1 Intensity', color='m')
    plt.tight_layout()

# %% finding cc differences
import math
from scipy import stats
def find_ccdiffs(cell_path: str, aamyo1_path: str, thresh: float=0.5, camsamprate: int=2, consamprate: int=5):
    """
    Args: Paths to cell segm output csv, myo1 img csv. Param for tick function and frame to min unit conversions.
    Outputs: Tuple of arrays containing the frame differences in cc progression between myo1 signal and manual annotations
    """
    # nov0525_ignore=[490,409,472,192,323,304,434,468,468,482,524,453,465]
    # ignore_list=[]
    # acdcdf=pd.read_csv(cell_path)
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
                            # if start_diff_1>32:
                            #     print('Start'+str((cell_id, roi_id, start_diff_1)))
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
                                # if start_diff_2>32:
                                #     print('Start'+str((cell_id, roi_id, start_diff_2)))
                        

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
                            # if 45>div_diff_1>22:
                            #     print('Div'+str((cell_id, roi_id, div_diff_1)))
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

                            #    if 45>div_diff_2>22:
                            #     print('Div'+str((cell_id, roi_id, div_diff_2)))

    print(len(myo1start), len(myo1div))
    diff_df=pd.concat(dfs, ignore_index=True)
    # return myo1start, myo1div
    return diff_df
    # return all_cells, ticksdf
    # return g1_moms
# %% Plot cc uncertainty histograms
def plot_distr(diff_df: pd.DataFrame):
    """
    Args: Dataframe containing difference measurements between annotations and myo1 dynamics.
    Outputs: Two plots illustrating the distribution of said measurements.
    """
    ms=diff_df.loc[diff_df.cc_group=='mom', 'start_diff'].dropna().values
    bs=diff_df.loc[diff_df.cc_group=='bud', 'start_diff'].dropna().values
    md=diff_df.loc[diff_df.cc_group=='mom', 'div_diff'].dropna().values
    bd=diff_df.loc[diff_df.cc_group=='bud', 'div_diff'].dropna().values
    
    # ms = [i for i in ms if 45>i]
    # md = [i for i in md if 30>i]
    # ms, md=start_arr, div_arr
    # fig,axes=plt.subplots(nrows=1,ncols=2,sharex=True)
    plt.figure()
    plt.xlabel('Time min')
    plt.title('Relative difference between myo1 on/off and division annotations')
    n0,b0,p0=plt.hist(md,bins='auto', alpha=0.75, edgecolor='black', color='b', label='Divisions (n='+f'{len(md)+len(bd)})')
    # n1,b1,p1=plt.hist(bd,bins='auto', alpha=0.75, edgecolor='black', color='orange', label='buds')
    plt.legend()
    plt.xlim(0,30)
    plt.vlines(np.mean(md), ymin=0, ymax=max(n0),colors='k', linestyle='dashed')
    # plt.vlines(np.mean(bd), ymin=0, ymax=max(n0),colors='orange', linestyle='dashed')

    plt.figure()
    plt.xlabel('Time (min)')
    plt.title('Relative difference between myo1 on and bud emergence annotations})')
    n2,b2,p2=plt.hist(np.concatenate([ms,bs]),bins='auto', alpha=0.75, edgecolor='black', color='b', label='Bud Emergence (n='+f'{len(ms)+len(bs)})')
    # n2,b2,p2=plt.hist(ms,bins='auto', alpha=0.75, edgecolor='black', color='b', label='moms')
    # n3,b3,p3=plt.hist(bs,bins='auto', alpha=0.75, edgecolor='black', color='orange', label='buds')
    plt.legend()
    plt.xlim(0,45)
    plt.vlines(np.mean(np.concatenate([ms,bs])), ymin=0, ymax=max(n2),colors='k', linestyle='dashed')
    # plt.vlines(np.mean(bs), ymin=0, ymax=max(n2),colors='orange', linestyle='dashed')
    

# %% Run the codes
# ggg,sss=analyze_fullcc(np.unique(list_cell), list_in)
# for org in orgs:
#     plot_fullcc(ggg,sss,org)