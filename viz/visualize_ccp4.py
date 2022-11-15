# %%
import numpy as np
import plotly.graph_objects as go
from plotly.offline import plot
import gemmi
import os

file_dir = os.path.dirname(__file__)
# s = input("give the spacing of desired grid:\n")
# t = input("give the energy thershold used in the grid calc:\n")
s, t = "0.1","100"
structure = "PSI"
structure = input("Structure name:\n")
grid_path = "../grid/%s_%s_%s.ccp4"%(structure,s,t)

# map generated by the binary of the code src/grid.cpp
ccp4 = gemmi.read_ccp4_map(os.path.join(file_dir,grid_path))
ccp4.setup(2)
arr = np.array(ccp4.grid, copy=False)

# %%
# Cubic box
# orth_mat = np.array(ccp4.grid.unit_cell.orthogonalization_matrix)
# cube_len = max(20,orth_mat[0,0],orth_mat[1,1],orth_mat[2,2])
cube_len = 20 # angstrom
spacing = 0.3
num = int(cube_len//spacing)
x = np.linspace(-15, cube_len, num=num)
y = np.linspace(-15, cube_len, num=num)
z = np.linspace(-15, cube_len, num=num)
X, Y ,Z = np.meshgrid(x, y, z,indexing='ij')
pts = np.column_stack( (X.flatten(), Y.flatten(), Z.flatten()))

values = []
for u in x:
    for v in y:
        for w in z:
            values.append(ccp4.grid.interpolate_value(gemmi.Position(u,v,w)))

# Plotly plot
fig = go.Figure(data=go.Volume(
    x=X.flatten(),
    y=Y.flatten(),
    z=Z.flatten(),
    value=np.array(values),
    isomin=np.min(arr),
    isomax=100,
    opacity=0.1, # needs to be small to see through all surfaces
    surface_count=17, # needs to be a large number for good volume rendering
    colorscale='rdylbu',reversescale=False # Choose color
    ))
plot(fig, filename=os.path.join(file_dir,"%s.html"%structure))
# fig.show()

