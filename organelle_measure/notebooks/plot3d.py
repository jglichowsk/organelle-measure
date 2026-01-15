# This script/notebook test how to tune the 3d Scatter plot, especially how to 
# set the background to white and add bounding box.

from matplotlib import projections
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path

fig = go.Figure()
for i in range(3):
    fig.add_trace(
        go.Scatter3d(
            x=np.random.random(10),
            y=np.random.random(10),
            z=np.random.random(10),
            name=f"{i}", mode='markers',
            marker=dict(size=2,opacity=0.8)
        )
        
    )
fig.update_layout(
    scene=dict(
        xaxis=dict(
            title=f"test_{i}",
            backgroundcolor='rgba(1,1,1,0)',
            gridcolor="grey",
            showline=True
        ),
        yaxis=dict(
            title=f"test_{i}",
            backgroundcolor='rgba(1,1,1,0)',
            gridcolor="grey",
            showline=True
        ),
        zaxis=dict(
            title=f"test_{i}",
            backgroundcolor='rgba(1,1,1,0)',
            gridcolor="grey",
            showline=True
        )
    )
)
fig.show()

import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(projection="3d")
for i in range(3):
    x=np.random.random(10),
    y=np.random.random(10),
    z=np.random.random(10),
    ax.scatter(x,y,z)

ax.set_xlabel("proj 0")
ax.set_ylabel("proj 1")
ax.set_zlabel("proj 2")

ax.xaxis.pane.set_edgecolor('black')
ax.yaxis.pane.set_edgecolor('black')
ax.zaxis.pane.set_edgecolor('black')

ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False

ax.xaxis.set_major_locator(plt.MultipleLocator(0.2))
ax.yaxis.set_major_locator(plt.MultipleLocator(0.2))
ax.zaxis.set_major_locator(plt.MultipleLocator(0.2))

ax.grid(which='major', axis='x', linewidth=0.75, linestyle=':', color='0.25')
ax.grid(which='major', axis='y', linewidth=0.75, linestyle=':', color='0.25')
ax.grid(which='major', axis='z', linewidth=0.75, linestyle=':', color='0.25')

fig.show()
