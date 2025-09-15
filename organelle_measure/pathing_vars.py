# This file contains all of the pathing variables which will be called in the processing scripts to:
#     1. Locate & read-in data.
#     2. Export data to the appropriate folders.
#     3. Create folders when needed.
#     4. Parse experiment details for Excel labelling. 

## Master path to experiment images
master_path = 
## Extension to folder of desired experiment. At the end of this path should be folders of different conditions for the given species.
experiment_path  = 

## List of experiment folders to operate on. Can take all folders in designated path or specify a subset.
import os
# folders_list = [f for f in os.listdir(str(master_path+'/'+experiment_path))]
# folders_list=['6-26_DNM1', '6-28_FZO1', '6-28_RTN1']
# folders_list=['insert folder names here']