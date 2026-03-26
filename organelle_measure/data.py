import numpy as np
import pandas as pd
from pathlib import Path

def read_results(folder_i,subfolders,pixel_sizes,path_rate=None):
    px_x,px_y,px_z = pixel_sizes

    organelles = [
        "px",
        "vo",
        "er",
        "gl",
        "mt",
        "ld"
    ]


    dfs_cell = []
    dfs_orga = []
    for folder in subfolders:
        df_folder_cell = pd.concat((pd.read_csv(fcell) for fcell in (Path(str(folder_i+'/'+folder+'/cell_measure'))).glob("*.csv")))
        df_folder_orga = pd.concat((pd.read_csv(fcell) for fcell in (Path(str(folder_i+'/'+folder+'/org_measure'))).glob("*.csv")))
        # where folder_i is imgs_path and folder is experiment path
        
        df_folder_cell["folder"] = folder
        df_folder_orga["folder"] = folder
        
        dfs_cell.append(df_folder_cell)
        dfs_orga.append(df_folder_orga)
    df_cell_all = pd.concat(dfs_cell)
    df_orga_all = pd.concat(dfs_orga)

    df_cell_all["condition"] = df_cell_all["condition"].apply(lambda x:float(str(x).replace('-',".")))
    df_orga_all["condition"] = df_orga_all["condition"].apply(lambda x:float(str(x).replace('-',".")))

    type_cell = {
        "folder":     "string",
        "experiment": "string",
        "condition":  "float",
        "hour":       "int8",
        "field":      "string",
        "idx_cell":   "int16",
        "area":       "int16",
        "bbox-0":     "int16",
        "bbox-1":     "int16",
        "bbox-2":     "int16",
        "bbox-3":     "int16",
    }
    type_orga = {
        "folder":       "string",
        "experiment":   "string",
        "condition":    "float",
        "hour":         "int8",
        "field":        "string",
        "organelle":    "string",
        "idx_cell":     "int16",
        "idx-orga":     "int16",
        "volume-pixel": "int16",
        "volume-bbox":  "int16",
        "bbox-0":       "int16",
        "bbox-1":       "int16",
        "bbox-2":       "int16",
        "bbox-3":       "int16",
        "bbox-4":       "int16",
        "bbox-5":       "int16"
    }

    df_cell_all = df_cell_all.astype(type_cell)
    df_orga_all = df_orga_all.astype(type_orga)


    # GROUP BY CELL
    df_cell_all.loc[:,"effective-volume"] = (px_x*px_y)*np.sqrt(px_x*px_y)*(2.)*df_cell_all.loc[:,"area"]*np.sqrt(df_cell_all.loc[:,"area"])/np.sqrt(np.pi) 
    pivot_cell_bycell = df_cell_all.set_index(["folder","condition","field","idx_cell"])

    # data (in unit of microns)
    df_orga_all["volume-micron"] = np.empty_like(df_orga_all.index)

    df_orga_all.loc[df_orga_all["organelle"].eq("vacuole"),"volume-micron"] = None #handles warning brought about by following line and filling an empty column/row in pandas (I believe).
    df_orga_all.loc[df_orga_all["organelle"].eq("vacuole"),"volume-micron"] = (px_x*px_y)*np.sqrt(px_x*px_y)*(2.)*df_orga_all.loc[df_orga_all["organelle"].eq("vacuole"),"volume-pixel"]*np.sqrt(df_orga_all.loc[df_orga_all["organelle"].eq("vacuole"),"volume-pixel"])/np.sqrt(np.pi) 
    df_orga_all.loc[df_orga_all["organelle"].ne("vacuole"),"volume-micron"] = px_x*px_y*px_z*df_orga_all.loc[df_orga_all["organelle"].ne("vacuole"),"volume-pixel"]

    pivot_orga_bycell_mean = df_orga_all.loc[:,["folder","condition","field","organelle","idx_cell","volume-micron"]].groupby(["folder","condition","field","idx_cell","organelle"]).mean()["volume-micron"]
    pivot_orga_bycell_nums = df_orga_all.loc[:,["folder","condition","field","organelle","idx_cell","volume-micron"]].groupby(["folder","condition","field","idx_cell","organelle"]).count()["volume-micron"]
    pivot_orga_bycell_totl = df_orga_all.loc[:,["folder","condition","field","organelle","idx_cell","volume-micron"]].groupby(["folder","condition","field","idx_cell","organelle"]).sum()["volume-micron"]
    # index
    index_bycell = pd.MultiIndex.from_tuples(
        [(*index,orga) for index in pivot_cell_bycell.index.to_list() for orga in [*organelles,'non-organelle']],
        names=['folder','condition','field','idx_cell','organelle']
    )
    pivot_bycell = pd.DataFrame(index=index_bycell)
    # combine data with index
    pivot_bycell.loc[pivot_orga_bycell_mean.index,"mean"] = pivot_orga_bycell_mean
    pivot_bycell.loc[pivot_orga_bycell_mean.index,"count"] = pivot_orga_bycell_nums
    pivot_bycell.loc[pivot_orga_bycell_mean.index,"total"] = pivot_orga_bycell_totl
    pivot_bycell.fillna(0.,inplace=True)

    # include cell data
    pivot_bycell.reset_index("organelle",inplace=True) # comment out after 1st run 
    pivot_bycell.loc[:,"cell-area"] = pivot_cell_bycell.loc[:,"area"]
    pivot_bycell.loc[:,"cell-volume"] = pivot_cell_bycell.loc[:,"effective-volume"]
    pivot_bycell.loc[:,"total-fraction"] = pivot_bycell.loc[:,"total"]/pivot_bycell.loc[:,"cell-volume"]

    # calculate properties of regions that are not organelles
    df_bycell = pivot_bycell.reset_index()

    df_none = df_bycell[df_bycell['organelle'].ne("non-organelle")].groupby(['folder','condition','field','idx_cell'])[['total','cell-volume','total-fraction']].agg({'total':'sum','cell-volume':'first','total-fraction':'sum'})
    pivot_bycell.loc[pivot_bycell['organelle'].eq("non-organelle"),"count"] = 1
    pivot_bycell.loc[pivot_bycell['organelle'].eq("non-organelle"),"mean"] = df_none["cell-volume"] - df_none["total"]
    pivot_bycell.loc[pivot_bycell['organelle'].eq("non-organelle"),"total"] = df_none["cell-volume"] - df_none["total"]
    pivot_bycell.loc[pivot_bycell['organelle'].eq("non-organelle"),"total-fraction"] = 1 - df_none["total-fraction"]

    df_bycell = pivot_bycell.reset_index()

    if path_rate is not None:
        df_rates = pd.read_csv(str(path_rate))
        df_rates.rename(columns={"experiment":"folder"},inplace=True)
        df_bycell.set_index(["folder","condition"],inplace=True)
        df_rates.set_index(["folder","condition"],inplace=True)
        df_bycell.loc[:,"growth_rate"] = df_rates["growth_rate"]
        df_bycell.reset_index(inplace=True)

    return df_bycell

