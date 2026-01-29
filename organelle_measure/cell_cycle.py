# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# %%
## Separate cells within dataframe into two groups: those that start in G1 and those in S/G2/M. For each, the next major cell cycle checkpoint
    # serves as the reference point. For example, cells that start in G1 will have the G1/S checkpoint as their reference point for binning
    # organelle measurements in time, and as an alignment point if/when applicable. 

## So split into two groups, then within each determine frames between which cell cycle checkpoint occurs, use this to determine waiting times,
    # and then use those for plots & analysis.

## Find and replace g1_ and s_ with G1_ and S_ respectively after done typing, if still want to. 
# %%
def pp_acdc_output(df: pd.DataFrame) -> pd.DataFrame:
    """
    Args: Dataframe containing Cell-ACDC output metrics and analysis information.
    Outputs: Dataframe containing the same information sans dead or manually-excluded cells.
    """
    # dead_indices=df[df['is_cell_dead']==1].index
    excl_indices=df[df['is_cell_excluded']==1].index
    pp_df=df.drop(excl_indices)
    pp_df.columns=pp_df.columns.str.strip() #remove leading and trailing spaces
    return pp_df

def cc_sort(cell_path: str) -> (pd.DataFrame, pd.DataFrame):
    """
    Args: Path to csv containing cell cycle, cell size, and organelle metric information. 
    Outputs: Two dataframes produced by splitting the input according to initial cell cycle position.
    """
    #read in cell dataframe
    raw_df=pd.read_csv(cell_path)
    #preprocess raw outputs (remove unwanted cells)
    cell_df=pp_acdc_output(raw_df)
    #group cells by starting cell cycle position. Those mothers starting in G1 or as mother go into their respective
    #categories. Buds are split between the two based upon whether detection occured in frame_0 or at any later time.
    g1_cells=[]
    s_cells=[]
    for cell_id in np.unique(cell_df['Cell_ID'].values):
        cell_rows=cell_df[cell_df['Cell_ID']==cell_id]

        if cell_rows['cell_cycle_stage'].values[0]=='G1':
            g1_cells.append(cell_rows)

        elif cell_rows['relationship'].values[0]=='mother':
            excl_indices=cell_rows[cell_rows['will_divide']==0].index #exclude first-frame mothers that will not divide 
            cell_rows=cell_rows.drop(excl_indices) 
            s_cells.append(cell_rows)

        elif cell_rows['relationship'].values[0]=='bud':
            if cell_rows['frame_i'].values[0]==0:
                s_cells.append(cell_rows)
            elif cell_rows['frame_i'].values[0]>0:
                g1_cells.append(cell_rows)
            else:
                Print('Exception: Bud first registered outside of recorded frames?')
        
        else:
            print('Error: Cell state not recognized for cell ' + str(cell_id))

    g1_output=pd.concat(g1_cells, ignore_index=True)
    s_output=pd.concat(s_cells, ignore_index=True)

    return g1_output, s_output

def cc_time_analysis(cell_df_path: str, org_df_path: str, plot_graph: bool=True) -> (pd.DataFrame,pd.DataFrame):
    """
    Args: Paths to csv files containing cell and organelle information respectively.
    Outputs: Two dataframes containing input info along with metric for alignment according to their respective cc checkpoint (G1/S 
            for G1 df, M/G1 for S df). Optional plot(s)
    """
    g1_df,s_df=cc_sort(cell_df_path)
    org_df=pd.read_csv(org_df_path)
    org_label=str(org_df['organelle'][0])

    #Find the cell cycle transition frame for each group of cells, and append information to dataframe.
    g1_dfs=[]
    for cell_id in g1_df.loc[g1_df.frame_i==0,'Cell_ID'].values: #for each G1 cell (mothers) present in frame 0...
        cell_rows=g1_df.loc[g1_df.Cell_ID==cell_id].copy()
        if len(cell_rows.loc[cell_rows.cell_cycle_stage=='S', 'frame_i'].values>0): #if G1/S transition happens...
            alignment_frame=cell_rows.loc[cell_rows.cell_cycle_stage=='S', 'frame_i'].values[0]
            relative_index=[(cell_rows['frame_i'].values[i]-alignment_frame) for i in range(len(cell_rows))]
            cell_rows.loc[cell_rows.Cell_ID==cell_id, 'Relative Index']=relative_index
            g1_dfs.append(cell_rows)
    g1_mothers=pd.concat(g1_dfs)

    s_dfs=[]
    for cell_id in s_df.loc[s_df.frame_i==0,'Cell_ID'].values: #for each S cell (mothers) present in frame 0...
        if s_df.loc[s_df.Cell_ID==cell_id, 'relationship'].values[0]=='mother':
            cell_rows=s_df.loc[s_df.Cell_ID==cell_id].copy()
            if len(cell_rows.loc[cell_rows.cell_cycle_stage=='G1', 'frame_i'].values>0): #if M/G1 transition happens...
                alignment_frame=cell_rows.loc[cell_rows.cell_cycle_stage=='G1', 'frame_i'].values[0]
                relative_index=[(cell_rows['frame_i'].values[i]-alignment_frame) for i in range(len(cell_rows))]
                cell_rows.loc[cell_rows.Cell_ID==cell_id, 'Relative Index']=relative_index
                s_dfs.append(cell_rows)
    s_mothers=pd.concat(s_dfs)
    
    #Plot desired organelle metric as function of cc position wrt reference point for both groups.
    if plot_graph==True:
        fig,axes=plt.subplots(nrows=1,ncols=2)

        for cell_id in np.unique(g1_mothers['Cell_ID'].values):
            start_offset=g1_mothers.loc[g1_mothers.Cell_ID==cell_id, 'Relative Index'].values[0]
            org_vol=org_df.loc[org_df.idx_cell==cell_id, 'volume-pixel']
            axes[0].scatter(start_offset,org_vol,c='m',marker='*') #scatter of organelle volume fraction versus relative frame index
            
        for cell_id in np.unique(s_mothers['Cell_ID'].values):
            division_offset=s_mothers.loc[s_mothers.Cell_ID==cell_id, 'Relative Index'].values[0]
            org_vol=org_df.loc[org_df.idx_cell==cell_id, 'volume-pixel']
            axes[1].scatter(division_offset, org_vol,c='m',marker='*') #scatter of organelle volume fraction versus relative frame index
        
        axes[0].set_title('Pre-Start '+org_label+' vs CC position')
        axes[0].set_xlabel('Frames before Start (5 min interval)')
        axes[0].set_ylabel(org_label+' volume (voxels)')
        axes[1].set_title('Pre-division '+org_label+' vs CC position')
        axes[1].set_xlabel('Frames before Division (5 min interval)')
        axes[1].set_ylabel(org_label+' volume (voxels)')
        fig.tight_layout()

    return g1_mothers, s_mothers


# %%
