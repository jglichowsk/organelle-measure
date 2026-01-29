import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from torch import nn
from torch import optim
from torch.utils.data import TensorDataset,DataLoader,random_split
from torch.utils.tensorboard import SummaryWriter
from pathlib import Path
from organelle_measure.data import read_results


# Chores: Using GPU, TensorBoard
print("Is CUDA available?: ",torch.cuda.is_available())
dev = torch.device("cuda") if torch.cuda.is_available() else torch.device("cpu")

writer = SummaryWriter('data/tensorboard/draft')


# Data
def get_data(train_ds, valid_ds, bs):
    return (
        DataLoader(train_ds, batch_size=bs, shuffle=True),
        DataLoader(valid_ds, batch_size=bs * 2),
    )


# Model
def quadratic_form(vec,matrix):
    return vec @ matrix @ vec

class CellSizePredictor(nn.Module):
    """ 
    input should have (2*num_organelle) dimensions,
    first half to be average volume, second half numbers 
    """
    def __init__(self,n_input):
        """
        output = Sum_i(A_i*x_i)+Sum_ij(x_i*B_ij*x_j)+C 
        """
        super().__init__()
        self.n_input = n_input
        self.A = nn.Parameter(torch.randn(self.n_input))
        self.B = nn.Parameter(torch.randn((self.n_input,self.n_input)))
        self.C = nn.Parameter(torch.zeros(1))
        self.alpha = nn.Parameter(torch.ones(self.n_input))

    def forward(self,x):
        # assert x.shape[0]==2*self.n_input, "Input Dimension Error!"
        average = x[:,:self.n_input]
        numbers = x[:,-self.n_input:]
        organelle = numbers * torch.pow(average,self.alpha)
        return (
            organelle @ self.A + 
            torch.stack(
                [quadratic_form(vec,self.B) for vec in torch.unbind(organelle,dim=0)],
                dim=0
            ) + 
            self.C
        )

class OrganelleClassifier(nn.Module):
    def __init__(self,n_input,n_output) -> None:
        super().__init__()
        self.n_input  = n_input
        self.n_output = n_output
        self.A = nn.Parameter((torch.randn(self.n_output,self.n_input)))
        self.B = nn.Parameter(torch.randn((self.n_input,self.n_output,self.n_input)))
        self.C = nn.Parameter(torch.zeros(self.n_output))
        self.alpha = nn.Parameter(torch.ones(self.n_input))

    def forward(self,x):
        # assert x.shape[0]==2*self.n_input, "Input Dimension Error!"
        average = x[:,:self.n_input]
        numbers = x[:,-self.n_input-1:-1]
        cellvol = x[-1]
        organelle = numbers * torch.pow(average,self.alpha) / cellvol
        # raw = (
        #     torch.tensordot(self.A,organelle,dims=1) + 
        #     torch.tensordot(
        #         torch.tensordot(organelle,self.B,dims=1),
        #         organelle,
        #         dims=1
        #     ) + 
        #     self.C
        # )
        return raw/torch.sum(raw)
# model = OrganelleClassifier(6,4)
# writer.add_graph(model=model,input_to_model=torch.zeros(13))

# Pipeline

def loss_batch(model, loss_func, xb, yb, opt=None):
    loss = loss_func(model(xb), yb)

    if opt is not None:
        loss.backward()
        opt.step()
        opt.zero_grad()

    return loss.item(), len(xb)

# def accuracy_classify(out,yb):
#     predictions = torch.argmax(out,axis=1)
#     targets = torch.argmax(out,axis=1)
#     return (predictions==targets).float()/predictions.shape[0]

def fit(epochs, model, loss_func, opt, train_dl, valid_dl):
    for epoch in range(epochs):
        model.train()
        for xb, yb in train_dl:
            loss_train,_  = loss_batch(model, loss_func, xb, yb, opt)
            writer.add_scalar(tag="loss in a training batch",scalar_value=loss_train)
            print(f"training batch loss={loss_train}")

        model.eval()
        with torch.no_grad():
            losses, nums = zip(
                *[loss_batch(model, loss_func, xb, yb) for xb, yb in valid_dl]
            )
        val_loss = np.sum(np.multiply(losses, nums)) / np.sum(nums)
        writer.add_scalar(tag="loss in a validation batch",scalar_value=val_loss)

        print(f"{epoch=}", f"{val_loss=}")
    return None


# def __main__():

# hyperparameters
learning_rate = 0.1
epochs = 500
batch_size = 256
training_ratio = 0.80
training_type = "regress"
# training_type = "classify"

# data
px_x,px_y,px_z = 0.41,0.41,0.20
organelles = [
    "peroxisome",
    "vacuole",
    "ER",
    "golgi",
    "mitochondria",
    "LD"
]
subfolders = [
    "EYrainbow_glucose",
    "EYrainbow_glucose_largerBF",
    "EYrainbow_rapamycin_1stTry",
    "EYrainbow_rapamycin_CheckBistability",
    "EYrainbow_1nmpp1_1st",
    "EYrainbow_leucine_large",
    "EYrainbow_leucine",
    "EYrainbowWhi5Up_betaEstrodiol"
]
folder_i = Path("./data/results")
path_rate = Path("./data/growthrate/growth_rate.csv")
df_bycell = read_results(folder_i,subfolders,(px_x,px_y,px_z),path_rate=path_rate)

df_bycell.set_index(["folder","condition","field","idx_cell"],inplace=True)
idx2learn = df_bycell.loc[df_bycell["organelle"].eq("ER")].index
df2learn = pd.DataFrame(index=idx2learn)
col_x = []
for stat in ["mean","count"]:
    for organelle in organelles:
        col = f"{stat}-{organelle}"
        df2learn[col] = df_bycell.loc[df_bycell["organelle"].eq(organelle),stat]
        col_x.append(col)

df2learn["cell-volume"] = df_bycell.loc[df_bycell["organelle"].eq("ER"),"cell-volume"]
# col_x.append("cell-volume")

data_x = df2learn[col_x].to_numpy()
data_y = df2learn["cell-volume"].to_numpy()
data_x = torch.Tensor(data_x)
data_y = torch.Tensor(data_y)
# data_x.to(dev)
# data_y.to(dev)

dataset_all = TensorDataset(data_x,data_y)
len_all   = len(dataset_all)
len_train = int(len_all*training_ratio)
dataset_train,dataset_valid = random_split(
    dataset_all,
    lengths=[len_train,len_all-len_train]
)
train_dataloader, valid_dataloader = get_data(dataset_train,dataset_valid,batch_size)

# model
model = CellSizePredictor(n_input=6)
# model.to(dev)
loss_function = F.cross_entropy if training_type=="classify" else F.mse_loss
optimizer = optim.SGD(model.parameters(),lr=learning_rate,momentum=0.9)

# training
fit(epochs,model,loss_function,optimizer,train_dataloader,valid_dataloader)
