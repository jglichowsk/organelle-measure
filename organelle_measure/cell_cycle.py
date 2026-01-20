# %%
import numpy as np
import pandas as pd

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
        # print('cell_id = '+str(cell_id))
        cell_rows=cell_df[cell_df['Cell_ID']==cell_id]

        if cell_rows['cell_cycle_stage'].values[0]=='G1':
            g1_cells.append(cell_rows)

        elif cell_rows['relationship'].values[0]=='mother':
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

    g1_df=pd.concat(g1_cells, ignore_index=True)
    s_df=pd.concat(s_cells, ignore_index=True)

    return g1_df, s_df

def cc_time_analysis(cell_df: pd.DataFrame, org_df: pd.DataFrame) -> pd.DataFrame:
    """
    Args: Dataframes containing sorted cell info (via cc_sort) and their corresponding organelle metrics.
    Outputs: Single dataframe containing cell info aligned by cell. Optional plot(s)
    """
    
    return


# %%
