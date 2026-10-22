# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
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
zstepsize=.200 #in microns
xbound=30

divtime_shift = 15.07 #minutes #0.4 thresh
divtime_uncertainty = 4.55 #minutes

# %% pathing
# acdc_paths=os.listdir(expmt_path+'/cell_measure')

# julcpath=r'C:\Users\jglic\Downloads\07172026 myo1 high res\cell_measure\BF-timelapse_acdc_output_07172026_cpsam-d25-ms20.csv'
# julmaskpath=r'C:\Users\jglic\Downloads\07172026 myo1 high res\BF-timelapse_07172026_myo1-mlemon_glucose-2.0_fov1_segm-afftransf-expand.tif'
# julmyo1path=r'C:\Users\jglic\Downloads\07172026 myo1 high res\confyellow-jg_07172026_myo1-mlemon_glucose-2.0_fov1_despeckle_r=4.csv'
# julaapath=r'C:\Users\jglic\Downloads\07172026 myo1 high res\autoassigntable_r=4_norm.csv'

# juncpath=r'C:\Users\jglic\Downloads\06102026 myo1 high res\BF-timelapse_acdc_output_cpsam-d25.csv'
# junmaskpath=r'C:\Users\jglic\Downloads\06102026 myo1 high res\BF-timelapse_cpsam-d25_afftransf-TEST2.tif'
# junmyo1path=r'C:\Users\jglic\Downloads\06102026 myo1 high res\SUM-confyellow-jg_06102026_ey2795-myo1_glucose-2.0_fov1_despeckle_r=4.csv'
# junaapath=r'C:\Users\jglic\Downloads\06102026 myo1 high res\autoassigntable_r=4_norm.csv'

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
    tags=name.split('_')
    if tags[0]=='deconv':
        deconv, stk, time, organelle, date, strain, condition, field, probtag=name.split('_')
    else:
        stk, time, organelle, date, strain, condition, field, probtag=name.split('_')
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
# %% misc functions
#ceiling division
def ceildiv(a, b):
    return -(a // -b)
# Normalization function
def norm(arr: np.ndarray):
    return arr/np.max(arr)
# %% Generate input file path lists
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
    # print(path_acdc)

    list_in.append(path_in)
    list_cell.append(path_acdc)

# args = pd.DataFrame({
#     "path_in":   list_in,
#     "path_cell": list_cell
# })
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
        if dead_indices.empty==False:
            df=df.drop(dead_indices)

    if remove_excl==True:
        excl_indices=df[df['is_cell_excluded']==1].index
        if excl_indices.empty==False:
            df=df.drop(excl_indices)

    initial_size_thresh=250 #px^2   #### I guess better way to do this is to ignore that Start but still utilize future information. 
    large_bud_tally=0
    if bud_thresh==True:
        for cell_id in np.unique(df['Cell_ID'].values):
            if df.loc[df.Cell_ID==cell_id, 'frame_i'].values[0]>0 and df.loc[df.Cell_ID==cell_id, 'relationship'].values[0]=='bud' and df.loc[df.Cell_ID==cell_id, 'cell_area_pxl'].values[0]>initial_size_thresh:
                large_bud_indices=df.loc[df.Cell_ID==cell_id].index
                df=df.drop(large_bud_indices)
                large_bud_tally+=1
    # print(f"# bud rows dropped = {large_bud_tally}")

    for cell_id in np.unique(df['Cell_ID'].values):
        cell_rows=df.loc[df.Cell_ID==cell_id]
        if norm_area_col==True:
            norm_area=[(cell_rows['cell_area_pxl'].values[i]/np.max(cell_rows['cell_area_pxl'].values)) for i in range(len(cell_rows))]
            df.loc[df.Cell_ID==cell_id, 'cell_area_norm']=norm_area 
        if approx_vol_col==True:
            vol_estimate=[(df.loc[df.Cell_ID==cell_id, 'cell_area_pxl'].values[i]*df.loc[df.Cell_ID==cell_id, 'minor_axis_length'].values[i])*(cam_pxl_size**2)*(zstepsize) for i in range(len(cell_rows))]
            df.loc[df.Cell_ID==cell_id, 'approx_vol']=vol_estimate 
            vol_estimate_norm=vol_estimate/np.max(vol_estimate)
            df.loc[df.Cell_ID==cell_id, 'approx_vol_norm']=vol_estimate_norm

    # output_name=Path(acdc_path).stem + '_pp.csv'
    # output_path=Path(acdc_path).parent/f"{output_name}"
    # df.to_csv(output_path, index=False)
    return df
# %% Linear regression
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
def linreg(dataset: np.ndarray, n: int=1000, train_frac=0.5):
    """
    Args: 2D Array of data to perform regression upon, number of runs, fraction of data to use for training.
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
                print(f"Error: Cell state not recognized for bud {cell_id}")
        else:
            print(f"Error: Cell state not recognized for cell {cell_id}")

    g1_output=pd.concat(g1_cells, ignore_index=True)
    s_output=pd.concat(s_cells, ignore_index=True)

    return g1_output, s_output
    
# for path in acdc_paths:
#     gg,ss=cc_sort(expmt_path+'/cell_measure'+'/'+path)
#     print(len(np.unique(gg.loc[gg.relationship=='mother','Cell_ID'].values)))
#     print(len(np.unique(ss.loc[ss.relationship=='mother','Cell_ID'].values)))

# %% cc phase lengths for entire fov
def cc_lengths(cell_path: str, camsamprate: int=5):
    """
    Args: Path to cell segm output csv
    Outputs: Arrays of frame lengths of G1, S/G2/M phases and lists of complete cc-phase # distribution and composition.
    """
    g1_df,s_df=cc_sort(cell_path)
    cdf=pd.concat([g1_df, s_df])
    g1_lens, s_lens, cc_lens =[],[],[]
    ratios=[]
    ind=0
    for cell_id in np.unique(cdf.loc[cdf.relationship=='mother','Cell_ID'].values): #for each G1 mother
        cell_rows=cdf.loc[cdf.Cell_ID==cell_id].copy()
        plist=cell_rows['cell_cycle_stage'].values
        lens,phases,num_phases = extract_streaks(plist)

        if len(lens)>2: #lengths of individual phases
            del lens[0] #drop any streaks contacting the edges 
            del lens[-1]
            del phases[0]
            del phases[-1]
            for i in range(len(lens)):
                if phases[i]=='G1':
                    g1_lens.append(lens[i])
                else:
                    s_lens.append(lens[i])

        if 2<=len(lens)<=3: #now record lengths of full cell cycles
            cc_lens.append(lens[0]+lens[1])
        if 4<=len(lens)<=5:
            cc_lens.append(lens[0]+lens[1])
            cc_lens.append(lens[2]+lens[3])
        if 6<=len(lens)<=7:
            cc_lens.append(lens[0]+lens[1])
            cc_lens.append(lens[2]+lens[3])
            cc_lens.append(lens[4]+lens[5])

    return np.array(g1_lens)*camsamprate, np.array(s_lens)*camsamprate, np.array(cc_lens)*camsamprate, ratios

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
            fullcc+=(num_transitions//4)
            checkpoint_frame=0
            for i in range(num_transitions):
                if i==0:
                    if cell_rows['relationship'].values[0]=='bud':
                        checkpoint_frame = cell_rows['frame_i'].values[0]
                    else: #for mothers
                        checkpoint_frame = lens[0]
                    column_name='Start_Index_1'
                else:
                    if cell_rows['relationship'].values[0]=='bud':
                        checkpoint_frame=checkpoint_frame+lens[i-1]
                    else: #for mothers
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
        num_transitions = num_phases-1+ind

        if num_transitions<3:
            s_df=s_df.drop(cell_rows.index)
        else:
            fullcc+=(num_transitions//3)
            checkpoint_frame=0
            for i in range(num_transitions):
                checkpoint_frame=checkpoint_frame+lens[i]
                if i%2==0:
                    column_name='Div_Index_'+str(int(i/2 + 1))
                else: #for G1 phase
                    column_name='Start_Index_'+str(int((i-1)/2 + 1))
                
                column_vals=[(cell_rows['frame_i'].values[j]-checkpoint_frame) for j in range(len(cell_rows))]
                s_df.loc[s_df.Cell_ID==cell_id, column_name]=column_vals

    print(f"# full cell cycles = {fullcc}")

    # return g1_df, s_df
    return pd.concat([g1_df, s_df], ignore_index=True)

# %% calc cc lengths, position, and phase metrics
def calc_ccmetrics(cell_path: str, save_csv: bool=False, cam_samprate: int=5, con_samprate: int=60):
    """
    Args: Cell-ACDC Output csv path, other args.
    Outputs: Dictionary containing selected cc metric info.
    """
    cdfs=[]
    dds,bbs=[],[]
    # for cell_path in cell_dfpaths:
        # celldf=pd.read_csv(cell_path)
    # cell_path=cell_dfpaths[0]
    g1df_sort, sdf_sort = cc_sort(cell_path)
    celldf = find_full_transitions(g1df_sort, sdf_sort)
    path_parts=Path(cell_path).stem.split("_")
    date=path_parts[1]
    fov=path_parts[4]

    for cell_id in np.unique(celldf['Cell_ID'].values):
        cell_rows=celldf.loc[celldf.Cell_ID==cell_id]
        cell_metrics={
            "Cell_ID"           : [cell_id],
            "date"              : [date], 
            "fov"               : [fov]
        }

        start_1=cell_rows.loc[cell_rows.Start_Index_1==0, 'frame_i']
        if start_1.empty==False:
            start_1_frame=start_1.values[0]
        div_1_frame=cell_rows.loc[cell_rows.Div_Index_1==0, 'frame_i'].values[0]
        div_2=cell_rows.loc[cell_rows.Div_Index_2==0, 'frame_i']
        if div_2.empty==False:
            div_2_frame=div_2.values[0]
        start_2=cell_rows.loc[cell_rows.Start_Index_2==0, 'frame_i']
        if start_2.empty==False: #for use with cells that fully pass Start_1 to Start_2 but not Div_2
            start_2_frame=start_2.values[0]
        div_3=cell_rows.loc[cell_rows.Div_Index_3==0, 'frame_i']
        if div_3.empty==False:
            div_3_frame=div_3.values[0]
        start_3=cell_rows.loc[cell_rows.Start_Index_3==0, 'frame_i']
        if start_3.empty==False: 
            start_3_frame=start_3.values[0]
        # div_4=cell_rows.loc[cell_rows.Div_Index_4==0, 'frame_i']
        # if div_4.empty==False:
        #     div_4_frame=div_4.values[0]

        min_frame=ceildiv(np.min(cell_rows['time_minutes'].values), con_samprate)
        max_frame=ceildiv(np.max(cell_rows['time_minutes'].values), con_samprate)
        for frame in np.arange(min_frame, max_frame, 1):
            abs_time=int(frame*con_samprate)
            cell_frames=cell_rows.loc[cell_rows.time_minutes==int(abs_time), 'frame_i']
            if cell_frames.empty==False: #if cell existed at conf timepoint
                relationship=cell_rows.loc[cell_rows.time_minutes==int(abs_time), 'relationship'].values[0]
                approx_vol=cell_rows.loc[cell_rows.time_minutes==int(abs_time), 'approx_vol'].values[0]
                cell_frame=cell_frames.values[0]
                g1_ratio= -100000 #dummy assign
                s_ratio=  -100000

                if div_2.empty==False and cell_frame<=div_1_frame:
                    cc_length = (div_2_frame - div_1_frame) * cam_samprate
                    if (div_1_frame-cell_frame)<cc_length: #if w/in one cc before div_1 
                        frac_time = (abs_time - (div_1_frame * cam_samprate)) / (cc_length)  
                        frac_time = 1 - (np.abs(frac_time) % 1)
                        if start_1_frame<div_1_frame:
                            g1_ratio=(start_2_frame - div_1_frame)*cam_samprate/cc_length
                            s_ratio= (div_2_frame - start_2_frame)*cam_samprate/cc_length
                        else:
                            g1_ratio=(start_1_frame - div_1_frame)*cam_samprate/cc_length
                            s_ratio= (div_2_frame - start_1_frame)*cam_samprate/cc_length
                        if relationship=='bud':
                            cc_num=0
                        else:
                            cc_num=0   
                    else:
                        continue
                # elif div_2.empty==True:
                #     cc_length = (start_2_frame - start_1_frame) * cam_samprate
                #     frac_time = (abs_time - (div_1_frame * cam_samprate)) / (cc_length) 
                #     if -1<frac_time<0:
                #         frac_time = 1 - (np.abs(frac_time) % 1)
                #     elif 1<frac_time<(start_2_frame+(start_2_frame - start_1_frame)):
                #         frac_time = frac_time % 1
                #     elif frac_time<-1 or frac_time>(start_2_frame+(start_2_frame - start_1_frame)):
                #         continue
                #     cc_num=1  
                    # bbs.append(cc_length)
                elif div_2.empty==False and div_1_frame<cell_frame<=div_2_frame:
                    cc_length = (div_2_frame - div_1_frame) * cam_samprate
                    frac_time = (abs_time - (div_1_frame * cam_samprate)) / (cc_length)  
                    if start_1_frame>div_1_frame:
                        g1_ratio=(start_1_frame - div_1_frame)*cam_samprate/cc_length
                        s_ratio= (div_2_frame - start_1_frame)*cam_samprate/cc_length
                    else:
                        g1_ratio=(start_2_frame - div_1_frame)*cam_samprate/cc_length
                        s_ratio= (div_2_frame - start_2_frame)*cam_samprate/cc_length
                    cc_num=1
                elif div_3.empty==False and div_2_frame<cell_frame<=div_3_frame:
                    cc_length = (div_3_frame - div_2_frame) * cam_samprate
                    frac_time = (abs_time - (div_2_frame * cam_samprate)) / (cc_length)
                    if start_2.empty==False and start_2_frame>div_2_frame:
                        g1_ratio=(start_2_frame - div_2_frame)*cam_samprate/cc_length
                        s_ratio= (div_3_frame - start_2_frame)*cam_samprate/cc_length
                    else:
                        g1_ratio=(start_3_frame - div_2_frame)*cam_samprate/cc_length
                        s_ratio= (div_3_frame - start_3_frame)*cam_samprate/cc_length
                    cc_num=2
                # elif div_4.empty==False and div_3_frame<cell_frame<=div_4_frame:
                #     cc_length = (div_4_frame - div_3_frame) * cam_samprate
                #     frac_time = (abs_time - (div_3_frame * cam_samprate)) / (cc_length)
                #     cc_num=3
                # elif div_4.empty==False and cell_frame>div_4_frame:
                #     cc_length = (div_4_frame - div_3_frame) * cam_samprate
                #     frac_time = (abs_time - (div_3_frame * cam_samprate)) / (cc_length)
                #     frac_time = frac_time % 1
                #     cc_num=4
                else:
                    continue
                real_time = frac_time * cc_length

                more_cell_metrics={                        
                    "relationship"     : [relationship],
                    "abs_time"         : [abs_time],
                    "frac_time"        : [frac_time],
                    "cc_length"        : [cc_length],
                    "real_time"        : [real_time],
                    "cc_num"           : [cc_num],
                    "approx_vol"       : [approx_vol],
                    "cell_area_norm"   : [cell_rows.loc[cell_rows.time_minutes==int(abs_time), 'cell_area_norm'].values[0]],
                    "Relative_ID"      : [cell_rows.loc[cell_rows.time_minutes==int(abs_time), 'relative_ID'].values[0]],
                    "g1_ratio"         : [g1_ratio],
                    "s_ratio"          : [s_ratio]
                    }
                results=cell_metrics | more_cell_metrics
                cdfs.append(pd.DataFrame(results))
                
        if len(cdfs)>0:
            out_df=pd.concat(cdfs, ignore_index=True)

    #somewhat arbitrary cutoffs, but based upon normal grow periods and cc length distributions.
    # excl_indices=out_df.loc[(out_df.cc_length<75) | (out_df.cc_length>180)].index 
    # if excl_indices.empty==False:
    #     out_df=out_df.drop(excl_indices)
    out_df.sort_values(by=["frac_time", "Cell_ID"], inplace=True, ignore_index=True)

    if save_csv==True:
        out_df.to_csv(Path(expmt_path+'\cc_measure')/f"orgdf_{meta['date']}.csv",index=False)
    return out_df
# %% full cc organelle analysis
def analyze_fullcc(cell_dfpaths: list, org_dfpaths: list, save_csv: bool=False, cam_samprate: int=5, con_samprate: int=60):
    """
    Args: List of paths to raw acdc csvs and organelle measure csvs.
    Outputs: One df containing geometric cell, cell cycle, and organelle information. Option to save to csv file.
    """
    cdfs=[]
    dds,bbs=[],[]
    # for cell_path in cell_dfpaths:
        # celldf=pd.read_csv(cell_path)
    cell_path=cell_dfpaths[0]
    g1df_sort, sdf_sort = cc_sort(cell_path)
    celldf = find_full_transitions(g1df_sort, sdf_sort)
    path_parts=cell_path.stem.split("_")
    date=path_parts[1]
    fov=path_parts[4]

    for cell_id in np.unique(celldf['Cell_ID'].values):
        cell_rows=celldf.loc[celldf.Cell_ID==cell_id]
        cell_metrics={
            "Cell_ID"           : [cell_id],
            "date"              : [date], 
            "fov"               : [fov]
        }

        start_1_frame=cell_rows.loc[cell_rows.Start_Index_1==0, 'frame_i'].values[0] #for use in filtering negative time frac buds
        div_1_frame=cell_rows.loc[cell_rows.Div_Index_1==0, 'frame_i'].values[0]
        div_2=cell_rows.loc[cell_rows.Div_Index_2==0, 'frame_i']
        if div_2.empty==False:
            div_2_frame=div_2.values[0]
        start_2=cell_rows.loc[cell_rows.Start_Index_2==0, 'frame_i']
        if start_2.empty==False: #for use with cells that fully pass Start_1 to Start_2 but not Div_2
            start_2_frame=start_2.values[0]
        div_3=cell_rows.loc[cell_rows.Div_Index_3==0, 'frame_i']
        if div_3.empty==False:
            div_3_frame=div_3.values[0]
        start_3=cell_rows.loc[cell_rows.Start_Index_3==0, 'frame_i']
        if start_3.empty==False: 
            start_3_frame=start_3.values[0]
        div_4=cell_rows.loc[cell_rows.Div_Index_4==0, 'frame_i']
        if div_4.empty==False:
            div_4_frame=div_4.values[0]

        min_frame=ceildiv(np.min(cell_rows['time_minutes'].values), con_samprate)
        max_frame=ceildiv(np.max(cell_rows['time_minutes'].values), con_samprate)
        for frame in np.arange(min_frame, max_frame, 1):
            abs_time=int(frame*con_samprate)
            cell_frames=cell_rows.loc[cell_rows.time_minutes==int(abs_time), 'frame_i']
            if cell_frames.empty==False: #if cell existed at conf timepoint
                relationship=cell_rows.loc[cell_rows.time_minutes==int(abs_time), 'relationship'].values[0]
                approx_vol=cell_rows.loc[cell_rows.time_minutes==int(abs_time), 'approx_vol'].values[0]
                cell_frame=cell_frames.values[0]
                g1_ratio= -1 #dummy assign
                s_ratio=  -1

                if div_2.empty==False and cell_frame<=div_1_frame:
                    cc_length = (div_2_frame - div_1_frame) * cam_samprate
                    if (div_1_frame-cell_frame)<cc_length: #if w/in one cc before div_1 
                        frac_time = (abs_time - (div_1_frame * cam_samprate)) / (cc_length)  
                        frac_time = 1 - (np.abs(frac_time) % 1)
                        if start_1_frame<div_1_frame:
                            g1_ratio=(start_2_frame - div_1_frame)*cam_samprate/cc_length
                            s_ratio= (div_2_frame - start_2_frame)*cam_samprate/cc_length
                        else:
                            g1_ratio=(start_1_frame - div_1_frame)*cam_samprate/cc_length
                            s_ratio= (div_2_frame - start_1_frame)*cam_samprate/cc_length
                        if relationship=='bud':
                            cc_num=0
                        else:
                            cc_num=0
                        
                    else:
                        continue
                # elif div_2.empty==True:
                #     cc_length = (start_2_frame - start_1_frame) * cam_samprate
                #     frac_time = (abs_time - (div_1_frame * cam_samprate)) / (cc_length) 
                #     if -1<frac_time<0:
                #         frac_time = 1 - (np.abs(frac_time) % 1)
                #     elif 1<frac_time<(start_2_frame+(start_2_frame - start_1_frame)):
                #         frac_time = frac_time % 1
                #     elif frac_time<-1 or frac_time>(start_2_frame+(start_2_frame - start_1_frame)):
                #         continue
                #     cc_num=1  
                    # bbs.append(cc_length)
                elif div_2.empty==False and div_1_frame<cell_frame<=div_2_frame:
                    cc_length = (div_2_frame - div_1_frame) * cam_samprate
                    frac_time = (abs_time - (div_1_frame * cam_samprate)) / (cc_length)  
                    if start_1_frame>div_1_frame:
                        g1_ratio=(start_1_frame - div_1_frame)*cam_samprate/cc_length
                        s_ratio= (div_2_frame - start_1_frame)*cam_samprate/cc_length
                    else:
                        g1_ratio=(start_2_frame - div_1_frame)*cam_samprate/cc_length
                        s_ratio= (div_2_frame - start_2_frame)*cam_samprate/cc_length
                    cc_num=1
                elif div_3.empty==False and div_2_frame<cell_frame<=div_3_frame:
                    cc_length = (div_3_frame - div_2_frame) * cam_samprate
                    frac_time = (abs_time - (div_2_frame * cam_samprate)) / (cc_length)
                    if start_2.empty==False and start_2_frame>div_2_frame:
                        g1_ratio=(start_2_frame - div_2_frame)*cam_samprate/cc_length
                        s_ratio= (div_3_frame - start_2_frame)*cam_samprate/cc_length
                    else:
                        g1_ratio=(start_3_frame - div_2_frame)*cam_samprate/cc_length
                        s_ratio= (div_3_frame - start_3_frame)*cam_samprate/cc_length
                    cc_num=2
                elif div_4.empty==False and div_3_frame<cell_frame<=div_4_frame:
                    cc_length = (div_4_frame - div_3_frame) * cam_samprate
                    frac_time = (abs_time - (div_3_frame * cam_samprate)) / (cc_length)
                    cc_num=3
                elif div_4.empty==False and cell_frame>div_4_frame:
                    cc_length = (div_4_frame - div_3_frame) * cam_samprate
                    frac_time = (abs_time - (div_3_frame * cam_samprate)) / (cc_length)
                    frac_time = frac_time % 1
                    cc_num=4
                else:
                    continue
                real_time = frac_time * cc_length

                more_cell_metrics={                        
                    "relationship"     : [relationship],
                    "abs_time"         : [abs_time],
                    "frac_time"        : [frac_time],
                    "cc_length"        : [cc_length],
                    "real_time"        : [real_time],
                    "cc_num"           : [cc_num],
                    "approx_vol"       : [approx_vol],
                    "cell_area_norm"   : [cell_rows.loc[cell_rows.time_minutes==int(abs_time), 'cell_area_norm'].values[0]],
                    "Relative_ID"      : [cell_rows.loc[cell_rows.time_minutes==int(abs_time), 'relative_ID'].values[0]],
                    "g1_ratio"         : [g1_ratio],
                    "s_ratio"          : [s_ratio]
                    }
                results=cell_metrics | more_cell_metrics

                for org_path in org_dfpaths: 
                    orgdf=pd.read_csv(org_path)
                    org_rows=orgdf.loc[orgdf.idx_cell==cell_id].loc[orgdf.abs_time==abs_time]
                    org_vox=org_rows['volume-pixel'].dropna().values #drop any NaNs
                    # if len(org_vox)>0: #if org measurement exists at this time
                    meta=parse_meta_orgmeasure(org_path.stem)
                    org_label=meta['organelle']
                    org_frame=org_rows.loc[org_rows.abs_time==abs_time]
                    org_count=len(org_frame['idx-orga'].values)
                    if len(org_vox)==0:
                        org_vox_tot=0
                    elif len(org_vox)==1:
                        org_vox_tot=org_vox[0]
                    else:
                        org_vox_tot=np.sum(org_vox)
                    org_vol_frac=(org_vox_tot*(cam_pxl_size**2)*(zstepsize)) / approx_vol

                    if org_label=='vo' and len(org_vox)>0:
                        eccentricity=org_rows['eccentricity'].values[0]
                    else:
                        eccentricity=0
                    org_metrics={
                        org_label+"_vol"   : [org_vox_tot],
                        org_label+"_frac"  : [org_vol_frac],
                        org_label+"_count" : [org_count],
                        org_label+"_ecc"   : [eccentricity]
                        }
                
                    results = results | org_metrics
                    # else:
                    #     continue
                cdfs.append(pd.DataFrame(results))
                
        if len(cdfs)>0:
            out_df=pd.concat(cdfs, ignore_index=True)

    #somewhat arbitrary cutoffs, but based upon normal grow periods and cc length distributions.
    excl_indices=out_df.loc[(out_df.cc_length<75) | (out_df.cc_length>180)].index 
    if excl_indices.empty==False:
        out_df=out_df.drop(excl_indices)
    out_df.sort_values(by=["frac_time", "Cell_ID"], inplace=True, ignore_index=True)

    if save_csv==True:
        out_df.to_csv(Path(expmt_path+'\cc_measure')/f"orgdf_{meta['date']}.csv",index=False)
    return out_df
    # return bbs, dds

# %% plot full cc org dynamics
def extract_frame_metrics(cellrows: pd.DataFrame, org: str, metric: str):
    """
    Args: processed org df, org label, metric (frac, count, vol, ecc)
    Outputs: Two 1D arrays containing fractional time x-coords and chosen metric y-coords.
    """
    frac_times = cellrows['frac_time'].values
    real_times = cellrows['real_time'].values
    metric_arr = cellrows[org+'_'+metric].values

    return frac_times, real_times, metric_arr

def plot_fullcc(org_df: pd.DataFrame, org_label: str, xl: float=-.5, xr: float=1.5, ms: int=2):
    """
    Args: One df containing both G1 or S/G2/M info, organelle label, x-axis bound param, set markersize.
    Outputs: Plots of org metric vs cc position.
    """ 
    # metrics=['frac', 'count']
    if org_label=='vo':
        metrics=['frac', 'vol', 'ecc']
    else:
        metrics=['frac', 'vol', 'count']

    fig,axes=plt.subplots(nrows=2,ncols=len(metrics))
    tally=0
    for cell in np.unique(org_df['Cell_ID'].values):
        cell_rows=org_df.loc[org_df.Cell_ID==cell]
        date=cell_rows["date"].values[0]
        mom_rows=cell_rows.loc[cell_rows.relationship=="mother"]
        bud_rows=cell_rows.loc[cell_rows.relationship=="bud"]

        # if cell_rows['relationship'].values[0]=='mother':
        #     markerstyle='o'
        #     color='k'
        # else:
        #     markerstyle='s'
        #     color='r'

        for j in range(len(metrics)):
            frac_times, real_times, m_arr = extract_frame_metrics(mom_rows, org_label, metrics[j])
            if len(m_arr)>0:
                axes[0,j].errorbar(frac_times, m_arr, yerr=0, ls='none', c='b',marker='o',ms=ms)
                axes[1, j].errorbar(real_times, m_arr, yerr=0, ls='none', c='b',marker='o',ms=ms)
                if j==1:
                    tally+=len(m_arr)
            frac_times, real_times, m_arr = extract_frame_metrics(bud_rows, org_label, metrics[j])
            if len(m_arr)>0:
                axes[0,j].errorbar(frac_times, m_arr, yerr=0, ls='none', c='r',marker='s',ms=ms)
                axes[1, j].errorbar(real_times, m_arr, yerr=0, ls='none', c='r',marker='s',ms=ms)
                if j==1:
                    tally+=len(m_arr)

    for j in range(len(metrics)):
        axes[0,j].set_xlabel('Normalized CC \n Position (Div to Div)')
        axes[1,j].set_xlabel('Time Relative \n to Previous Division (min)')
        axes[0,j].set_title(f"Unbinned {org_label} {metrics[j]} \n vs CC position")

    fig.tight_layout()
    print(f"{org} n: {tally} data points")
    return
# %% Binned full cc org dynamics
def fullccbin(org_df: pd.DataFrame, org: str, binwidth: float):
    """
    Args: One df containing both G1 or S/G2/M info, binwidth in fractional units
    Outputs: Two output arrays containing [frac_time, org_vox, bin_std] for mother and daughter respectively.
    """
    bin_edges=np.linspace(np.min([min(org_df['frac_time'].values),0]), np.max([max(org_df['frac_time'].values),1]), int(1//binwidth))
    mom_arr, bud_arr=[],[]
    for i in range(len(bin_edges)-1):
        l_bound, r_bound = bin_edges[i], bin_edges[i+1]
        fracbin_rows=org_df.loc[(org_df.frac_time>=l_bound) & (org_df.frac_time<=r_bound)]
        mom_rows=fracbin_rows.loc[fracbin_rows.relationship=="mother"]
        bud_rows=fracbin_rows.loc[fracbin_rows.relationship=="bud"]
        frac_time=(l_bound+r_bound)/2

        m_arr=mom_rows[org+'_vol'].dropna().values #Drop NaNs from stats
        if len(m_arr)>0:
            mom_arr.append([frac_time, np.mean(m_arr), np.std(m_arr)])
        m_arr=bud_rows[org+'_vol'].dropna().values
        if len(m_arr)>0:
            bud_arr.append([frac_time, np.mean(m_arr), np.std(m_arr)])

    return np.array(mom_arr), np.array(bud_arr)
# %% plot binned full cc org dynamics
def plot_fullccbin(org_df: pd.DataFrame, org: str, bin_text: bool=False, binwidth: int=.05, ms: int=5):
    """
    Args: One df containing both G1 or S/G2/M info, organelle label, ARG ABOUT BIN WIDTH OR RAW VS FRAC TIME
    Outputs: Plots of binned org metric vs cc position.
    """
    #MAYBE ADD cc_num AND/OR DATE COLOR CODE INTO HERE
    #bin frac or real times by 5-10%, avg, plot error scatter

    if org=='vo':
        metrics=['frac', 'vol', 'ecc']
    else:
        metrics=['frac', 'vol', 'count']

    fig,axes=plt.subplots(nrows=2,ncols=len(metrics))
    bin_edges=np.linspace(np.min([min(org_df['frac_time'].values),0]), np.max([max(org_df['frac_time'].values),1]), int(1//binwidth))
    tally=0
    out_arr=[]
    for i in range(len(bin_edges)-1):
        l_bound, r_bound = bin_edges[i], bin_edges[i+1]
        fracbin_rows=org_df.loc[(org_df.frac_time>=l_bound) & (org_df.frac_time<=r_bound)]
        mom_rows=fracbin_rows.loc[fracbin_rows.relationship=="mother"]
        bud_rows=fracbin_rows.loc[fracbin_rows.relationship=="bud"]
        frac_time=(l_bound+r_bound)/2

        for j in range(len(metrics)):
            m_arr=mom_rows[org+'_'+metrics[j]].dropna().values #Drop NaNs from stats
            if len(m_arr)>0:
                axes[0,j].errorbar(frac_time, np.mean(m_arr), yerr=np.std(m_arr), xerr=0, ls='none', c='b',marker='o',ms=ms)
                out_arr.append([frac_time, np.mean(m_arr), np.std(m_arr)])
                if bin_text==True:
                    axes[0,j].text(frac_time, np.mean(m_arr)*1.25, int(len(m_arr)), ha="center", fontsize="small")            
                if j==1:
                    tally+=int(len(m_arr))
            m_arr=bud_rows[org+'_'+metrics[j]].dropna().values
            if len(m_arr)>0:
                axes[0,j].errorbar(frac_time, np.mean(m_arr), yerr=np.std(m_arr), xerr=0, ls='none', c='r',marker='s',ms=ms)
                if bin_text==True:
                    axes[0,j].text(frac_time, np.mean(m_arr)*1.25, int(len(m_arr)), ha="center", fontsize="small")
                if j==1:
                    tally+=int(len(m_arr))

    bin_edges=np.linspace(np.min([min(org_df['frac_time'].values), 0]), max(org_df['cc_length'].values), int(1//(binwidth)))
    for i in range(len(bin_edges)-1):
        l_bound, r_bound = bin_edges[i], bin_edges[i+1]
        realbin_rows=org_df.loc[(org_df.real_time>=l_bound) & (org_df.real_time<=r_bound)]
        mom_rows=realbin_rows.loc[realbin_rows.relationship=="mother"]
        bud_rows=realbin_rows.loc[realbin_rows.relationship=="bud"]
        real_time=(l_bound+r_bound)/2

        for j in range(len(metrics)):
            m_arr=mom_rows[org+'_'+metrics[j]].dropna().values #Drop NaNs from stats
            if len(m_arr)>0:
                axes[1,j].errorbar(real_time, np.mean(m_arr), yerr=np.std(m_arr), xerr=0, ls='none', c='b',marker='o',ms=ms)
                if bin_text==True:
                    axes[1,j].text(real_time, np.mean(m_arr)*1.25, int(len(m_arr)), ha="center", fontsize="small")
            m_arr=bud_rows[org+'_'+metrics[j]].dropna().values
            if len(m_arr)>0:
                axes[1,j].errorbar(real_time, np.mean(m_arr), yerr=np.std(m_arr), xerr=0, ls='none', c='r',marker='s',ms=ms)
                if bin_text==True:
                    axes[1,j].text(real_time, np.mean(m_arr)*1.25, int(len(m_arr)), ha="center", fontsize="small")

    for j in range(len(metrics)):
        axes[0,j].set_xlabel('Normalized CC \n Position (Div to Div)')
        axes[1,j].set_xlabel('Time Relative \n to Previous Division (min)')
        axes[0,j].set_title(f"Binned {org} {metrics[j]} \n vs CC position")
    fig.tight_layout()

    print(f"Number of {org} data points in binned graphs: {tally}")
    return np.array(out_arr)
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
    return pd.concat([g1_df, s_df])

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
        
        ###############g1_moms, s_moms, g1_buds, s_buds = find_all_transitions(g1_df, s_df) #find annotated cc transitions
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

# %% Plot binned all cc org dynamics
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

# %% plot cc sizes about transition ref points
def plot_ccsize(cell_dfpath: str, norm_size: bool=False, xbound: int=30):
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

# %% Run full cc org analysis
# odf=analyze_fullcc(np.unique(list_cell), list_in)
# for org in orgs:
#     #plot_fullcc(odf,org)
#     plot_fullccbin(odf, org)
# %% Plot cc lengths from analyze_fullcc() output
# plt.figure()
# date=odf['date'].values[0]
# # n,b,p=plt.hist(odf['cc_length'].values, bins='auto',edgecolor='black')
# n,b,p=plt.hist(odf.loc[odf.cc_num==1, 'cc_length'].values, bins='auto', color='b',alpha=0.75,edgecolor='black', label=f"CC #1")
# n,b,p=plt.hist(odf.loc[odf.cc_num==2, 'cc_length'].values, bins='auto', color='y',alpha=0.75,edgecolor='black', label=f"CC #2")
# plt.xlabel('Time (min)')
# plt.ylabel('Number of Cell Cycles')
# plt.title(f"Raw CC Lengths by CC No. for {date}")
# plt.legend(loc='best')

# %% Plot org metric vs cc position with cc length as color gradient
# fig,axes=plt.subplots(nrows=1,ncols=1)
# metrics=['count']
# org_label='ld'
# ms=5
# for cell in np.unique(odf['Cell_ID'].values):
#     cell_rows=odf.loc[odf.Cell_ID==cell]
#     date=cell_rows["date"].values[0]
#     mom_rows=cell_rows.loc[cell_rows.relationship=="mother"]
#     bud_rows=cell_rows.loc[cell_rows.relationship=="bud"]
#     frac_times, real_times, m_arr = extract_frame_metrics(mom_rows, org_label, metrics[0])
#     if len(m_arr)>0:
#         grad=m_arr
#         axes.scatter(frac_times, m_arr,s=5, c=grad,marker='o', cmap='Greys', edgecolors='k')
#     frac_times, real_times, m_arr = extract_frame_metrics(bud_rows, org_label, metrics[0])
#     if len(m_arr)>0:
#         grad=m_arr
#         mapp=axes.scatter(frac_times, m_arr,s=5, c=grad,marker='v', cmap='Greys', edgecolors='k')
# axes.set_xlabel('Normalized CC \n Position (Div to Div)')
# axes.set_ylabel(f"{org_label} {metrics[0]}")
# fig.colorbar(mapp,label='CC Length')
# https://matplotlib.org/stable/users/explain/colors/colorbar_only.html

# %% Calc and plot diffs between concurrent cell cycles for each cell.
# diffs, cc1, cc2=[],[],[]
# date=odf['date'].values[0]
# for cell_id in np.unique(odf['Cell_ID'].values):
#     cell_rows=odf.loc[odf.Cell_ID==cell_id]
#     cc_1=cell_rows.loc[cell_rows.cc_num==1, 'cc_length']
#     cc_2=cell_rows.loc[cell_rows.cc_num==2, 'cc_length']
#     if cc_2.empty==False and cc_1.empty==False:
        # cc1.append(cc_1.values[0])
        # cc2.append(cc_2.values[0])
#         diffs.append(cc_2.values[0]-cc_1.values[0])
# fig,axes=plt.subplots(nrows=1,ncols=2)
# axes[0].hist(diffs, bins='auto',edgecolor='black')
# axes[0].set_xlabel('Time (min)')
# axes[0].set_ylabel('CC_2 - CC_1 (min)')
# axes[1].scatter(cc1, cc2)
# axes[1].set_xlabel('CC 1 (min)')
# axes[1].set_ylabel('CC 2 (min)')
# fig.suptitle(f"CC Differences by CC Num for {date}")
# fig.tight_layout()

# %% Extract and plot phase comparisons.
# undf=odf.drop_duplicates(subset=['Cell_ID', 'cc_length']) #drops if both column entries same
# undf=undf.drop(undf[undf.g1_ratio<-100].index | undf[undf.s_ratio<-100].index)
## Plot G1 vs Total CC length scatter
# plt.figure()
# shapes=['o','s','v']
# colors=['m','r','k']
# for cc in range(max(undf['cc_num'].values)+1):
#     g1s=undf.loc[undf.cc_num==cc,'g1_ratio'].values
#     ccs=undf.loc[undf.cc_num==cc,'cc_length'].values
#     plt.scatter(ccs,g1s, c=colors[cc], s=20, marker=shapes[cc], label=f"cc {cc}, n ={len(g1s)}")
# plt.xlabel('Total CC Length')
# plt.ylabel('G1 Phase Ratio')
# plt.title(f"{undf['date'].values[0]}")
# plt.legend(loc='best')
# plt.ylim(0,1)
## Plot G1 vs Total CC length scatter BINNED
# plt.figure()
# shapes=['o','s','v']
# colors=['m','r','k']
# for cc in range(max(undf['cc_num'].values)+1):
#     cc_rows=undf.loc[undf.cc_num==cc]
#     cc_lengths=cc_rows['cc_length'].values
#     for cc_length in np.unique(np.sort(cc_lengths)):
#         g1s=cc_rows.loc[cc_rows.cc_length==cc_length,'g1_ratio'].values
#         plt.errorbar(cc_length,np.mean(g1s),np.std(g1s),c=colors[cc], ms=10, marker=shapes[cc])
#         # plt.text(cc_length, np.mean(g1s)*1.15, int(len(g1s)), ha="center", fontsize="large")
# plt.xlabel('Total CC Length')
# plt.ylabel('G1 Phase Ratio')
# plt.title(f"{undf['date'].values[0]}")
# plt.ylim(0,1)
# cc0 = mpatches.Patch(color=colors[0], label="cc 0")
# cc1 = mpatches.Patch(color=colors[1], label="cc 1")
# cc2 = mpatches.Patch(color=colors[2], label="cc 2")
# plt.legend(handles=[cc0, cc1, cc2], loc='best')
## Plot G1 vs S/G2/M phase scatter
# plt.figure()
# shapes=['o','s','v']
# colors=['m','r','k']
# for cc in range(max(undf['cc_num'].values)+1):
#     g1s=undf.loc[undf.cc_num==cc,'g1_ratio'].values
#     ss=undf.loc[undf.cc_num==cc,'s_ratio'].values
#     plt.scatter(g1s, ss, c=colors[cc], s=10, marker=shapes[cc], label=f"cc {cc}, n ={len(ss)}")
# plt.ylabel('S/G2/M Phase Ratio')
# plt.xlabel('G1 Phase Ratio')
# plt.title(f"{undf['date'].values[0]}")
# plt.legend(loc='best')
# plt.xlim(0,1)
# plt.ylim(0,1)
## Plot Individual cc phases histogram 
# plt.figure()
# shapes=['o','s','v']
# colors=['m','y','k']
# for cc in range(max(undf['cc_num'].values)+1):
#     g1s=undf.loc[undf.cc_num==cc,'g1_ratio'].values
#     ss=undf.loc[undf.cc_num==cc,'s_ratio'].values
#     n,b,p=plt.hist(g1s, bins=np.arange(0,1,0.05), density=True, alpha=0.75, color=colors[cc],edgecolor='black', label=f"cc {cc} G1 n={len(g1s)}")
#     # n,b,p=plt.hist(ss, bins='auto', alpha=0.75, color=colors[cc], edgecolor='black', label=f"cc {cc} S/G2/M n={len(ss)}")
# plt.ylabel('Prob Density')
# plt.xlabel('Ratio to CC Length')
# plt.title(f"{undf['date'].values[0]}")
# plt.legend(loc='best')
# %% Run cc lengths
# cpaths=[list_cell[0], r'C:\Users\jglic\Downloads\07152026 3c-ey2795\cell_measure\BF-timelapse_acdc_output_07152026_cpsam-d25-ms20.csv']
# g_lens,s_lens,lens_l=[],[],[]
# for i in range(len(cpaths)):
#     gl,sl,lens_dist, _=cc_lengths(cpaths[i])
#     g_lens.append(gl)
#     s_lens.append(sl)
#     lens_l.append(lens_dist)
# g_arr=np.concatenate(g_lens)
# s_arr=np.concatenate(s_lens)
# l_arr=np.concatenate(lens_l)

# %% Plot cc phase and total lengths 
# plt.figure()
# plt.xlabel("Time (min)")
# plt.title('Mother Cell CC Phase Annotation Lengths')
## bins=np.arange(0,max(cclen)+20,5)
# n0,b0,p0=plt.hist(s_arr,bins='auto', alpha=0.75, edgecolor='black', color='b', label='S/G2/M Phase n=('+f'{len(s_arr)})')
# n1,b1,p1=plt.hist(g_arr, bins='auto', alpha=0.75, edgecolor='black', color='y', label='G1 Phase n=('+f'{len(g_arr)})')
# n2,b2,p2=plt.hist(l_arr, bins='auto', alpha=0.75, edgecolor='black', color='k', label='Total CC n=('+f'{len(l_arr)})')
# plt.legend()

# %% plot initial bud size distribution
def plot_buddist(cdf: pd.DataFrame):
    """
    Args: Processed Cell ACDC dataframe
    Outputs: Numpy array of bud sizes at first detection
    """
    buds=cdf.loc[cdf.cc_group=='gbud']
    # buds=cdf.loc[cdf.relationship=='bud']
    bud_dist=[]
    for bud_id in np.unique(buds['Cell_ID'].values):
        bud_dist.append(buds.loc[buds.Cell_ID==bud_id, "cell_area_pxl"].values[0])
    fig,axes=plt.subplots(nrows=1,ncols=1)
    axes.hist(bud_dist, bins='auto', edgecolor='black', label='Buds')
    axes.set_title('Initial Bud Size Distribution')
    axes.set_xlabel('Cell Area (px^2)')
    fig.tight_layout()
    return bud_dist

# gg,ss=cc_sort(cdfpath)
# cdf=pd.concat([gg,ss])
# bd=plot_buddist(cdf)
# %% Run and Plot pairwise org parametric curves
# vo_m, vo_b = fullccbin(odf, 'vo', binwidth=0.05)
# ld_m, ld_b = fullccbin(odf, 'ld', binwidth=0.05)
# gl_m, gl_b = fullccbin(odf, 'gl', binwidth=0.05)
# fig,axes=plt.subplots(figsize=(6,6),nrows=3,ncols=1)
# im00=axes[0].scatter(vo_m[:,1], ld_m[:,1], c=vo_m[:,0], edgecolor='k', cmap='Reds', label='mothers')
# im01=axes[0].scatter(vo_b[:,1], ld_b[:,1], c=vo_b[:,0], marker='s', edgecolor='k', cmap='Reds', label='buds')
# axes[0].set_xlabel('Vacuole Volume (vox)')
# axes[0].set_ylabel('LD Volume (vox)')
# axes[0].legend(loc='upper left')
# fig.colorbar(im00, ax=axes[0], label='Normalized CC Position \n Div to Div')
# im10=axes[1].scatter(vo_m[:,1], gl_m[:,1], c=vo_m[:,0], edgecolor='k', cmap='Reds', label='mothers')
# im11=axes[1].scatter(vo_b[:,1], gl_b[:,1], c=vo_b[:,0], marker='s', edgecolor='k', cmap='Reds', label='buds')
# axes[1].set_xlabel('Vacuole Volume (vox)')
# axes[1].set_ylabel('Golgi Volume (vox)')
# axes[1].legend(loc='upper left')
# fig.colorbar(im10, ax=axes[1], label='Normalized CC Position \n Div to Div')
# im20=axes[2].scatter(gl_m[:,1], ld_m[:,1], c=vo_m[:,0], edgecolor='k', cmap='Reds', label='mothers')
# im21=axes[2].scatter(gl_b[:,1], ld_b[:,1], c=vo_b[:,0], marker='s', edgecolor='k', cmap='Reds', label='buds')
# axes[2].set_xlabel('Golgi Volume (vox)')
# axes[2].set_ylabel('LD Volume (vox)')
# axes[2].legend(loc='upper left')
# fig.colorbar(im20, ax=axes[2], label='Normalized CC Position \n Div to Div')
# fig.tight_layout()

# %% Organelle regression    
import itertools
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import PolynomialFeatures

def regr_pipeline(dataset: pd.DataFrame, deg: int, indep_vars: list):
    """
    Args: Dataframe including input and target data, other keyword args.
    Outputs: Arrays of R^2 score distribution mean & std for given regressors and outcomes.
    """
    scores=[]
    for n in range(1000):
        train, test = train_test_split(dataset,test_size=0.3)
        X_train,X_test = train[indep_vars], test[indep_vars]
        Y_train,Y_test = group_odf.loc[group_odf.index[train.index], outcomes[0]], group_odf.loc[group_odf.index[test.index], outcomes[0]]

        # regr = LinearRegression().fit(X_train, Y_train)
        polynomial_features = PolynomialFeatures(degree=deg, include_bias=False)
        linear_regression = LinearRegression()
        pipeline = Pipeline(
            [("polynomial_features", polynomial_features),
                ("linear_regression", linear_regression)])
        regr=pipeline.fit(X_train, Y_train)

        score=regr.score(X_test,Y_test)
        scores.append(score)
    output.append([np.mean(scores), np.std(scores)])
    return output
        # print(f"Independent Vars: {indep_vars}")
        # print(f"Score: {score}")
        # print('Coefficients: \n', regr.coef_)
        # print('Intercept: \n', regr.intercept_)
        # print('Coefficients: \n', regr.coef_)

###Need to add paramter output, especially for those with multiple regressor inputs.
#  Also org-org predictions. Also figure out whether the Pipeline with polynomial is doing what I would like it to be doing.
outcomes = ['frac_time']
# outcomes = ['real_time']
regressors = ['vo_vol', 'ld_vol', 'gl_vol']

def linregr(indiv_regressors: list, outcomes: list, deg: int=1, plot_graph=True):
    """
    Args: Lists of regressors and target outcomes, degree of polynomial features
    Outputs: Arrays of regression score stats after n runs and a graph.
    """
    combos = []
    for i in range(1, len(indiv_regressors)+1):
        els = [list(x) for x in itertools.combinations(indiv_regressors, i)]
        combos.extend(els)

    groups=['mother', 'bud']
    for group in groups: 
        output=[]
        group_odf=odf.loc[odf.relationship==group]
        group_odf.reset_index(drop=True, inplace=True)
        
        # ind=0
        for indep_vars in combos:
            dataset=group_odf[indep_vars]
            # output = regr_pipeline(dataset, deg, indep_vars)
            scores=[]
            for n in range(100):
                train, test = train_test_split(dataset,test_size=0.3)
                X_train,X_test = train[indep_vars], test[indep_vars]
                Y_train,Y_test = group_odf.loc[group_odf.index[train.index], outcomes[0]], group_odf.loc[group_odf.index[test.index], outcomes[0]]

                # regr = LinearRegression().fit(X_train, Y_train)
                polynomial_features = PolynomialFeatures(degree=deg, include_bias=False)
                linear_regression = LinearRegression()
                pipeline = Pipeline(
                    [("polynomial_features", polynomial_features),
                        ("linear_regression", linear_regression)])
                regr=pipeline.fit(X_train, Y_train)

                score=regr.score(X_test,Y_test)
                scores.append(score)
            output.append([np.mean(scores), np.std(scores)])

        if group=='mother':
            m_outarr=np.array(output)
        else:
            b_outarr=np.array(output)

    if plot_graph==True:
        fig, axes = plt.subplots(1, 1)
        xticks = np.linspace(0, len(m_outarr), len(m_outarr))
        axes.errorbar(xticks, m_outarr[:,0], m_outarr[:,1], fmt="o", label='mothers')
        axes.errorbar(xticks, b_outarr[:,0], b_outarr[:,1], fmt="o", c='r', label='buds')
        axes.set_title(f"Organelle Volumes as Regressors \n for Normalized CC Position deg={deg}")
        axes.set_ylabel('Average R^2 Value')
        axes.legend()
        # axes.set_ylim([0,1])
        xlabels = [" ".join(label) for label in combos]
        axes.set_xticks(xticks, labels=xlabels)
        axes.tick_params("x", rotation=-45)
    
    return m_outarr, b_outarr, combos

def plot_regrcomp(indiv_regressors: list, outcomes: list, deg: int=1):
    """
    Args:  Lists of regressors and target outcomes, degree of polynomial features
    Outputs: Plot overlaying the R^2 scores of the regressor combos on each of the target outcomes. 
    """
    fig, axes = plt.subplots(1, 1)
    markers=['o', 's', 'v', '<']
    colors=['m', 'k', 'r', 'y']
    for i in range(len(outcomes)):
        m_outarr, b_outarr, combos = linregr(indiv_regressors, [outcomes[i]], deg, plot_graph=False)
        xticks = np.linspace(0, len(m_outarr), len(m_outarr))
        axes.errorbar(xticks, m_outarr[:,0], m_outarr[:,1], fmt=markers[i], c=colors[i], label=f"{outcomes[i]} mothers")
        axes.errorbar(xticks, b_outarr[:,0], b_outarr[:,1], fmt=markers[i], c=colors[i], label=f"{outcomes[i]} buds")
    axes.set_title(f"Organelle Volumes as Regressors \n for Normalized CC Position deg={deg}")
    axes.set_ylabel('Average R^2 Value')
    axes.legend(fontsize='x-small')
    xlabels = [" ".join(label) for label in combos]
    axes.set_xticks(xticks, labels=xlabels)
    axes.tick_params("x", rotation=-45)
    return
# %% #including cross validation or some bs idk
# ## https://scikit-learn.org/stable/modules/cross_validation.html
# from sklearn.linear_model import LinearRegression
# from sklearn.model_selection import cross_val_score
# from sklearn.pipeline import Pipeline
# from sklearn.preprocessing import PolynomialFeatures

# polynomial_features = PolynomialFeatures(degree=2, include_bias=False)
# linear_regression = LinearRegression()
# pipeline = Pipeline(
#     [("polynomial_features", polynomial_features),
#         ("linear_regression", linear_regression)])
# pipeline.fit(X_train, Y_train)

# # Evaluate the models using crossvalidation
# scores = cross_val_score(
#     pipeline, X_test, Y_test, scoring="neg_mean_squared_error", cv=10
# )

# plt.figure(figsize=(14, 5))
# for i in range(len(degrees)):
# ax = plt.subplot(1, len(degrees), i + 1)
# plt.setp(ax, xticks=(), yticks=())
# X_test = np.linspace(0, 1, 100)
# plt.plot(X_test, pipeline.predict(X_test), label="Model")
# plt.scatter(X_test, Y_test, edgecolor="b", s=20, label="Samples")
# plt.xlabel("x")
# plt.ylabel("y")
# plt.xlim((0, 1))
# plt.ylim((-2, 2))
# plt.legend(loc="best")
# plt.title(
#     "Degree {}\nMSE = {:.2e}(+/- {:.2e})".format(
#         degrees[i], -scores.mean(), scores.std()
#     )
# )
# plt.show()
# %%
