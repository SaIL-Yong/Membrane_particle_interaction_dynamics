# mesh_modification.py
# Upsample (refine) and coarsen (decimate) a mesh using libigl’s Python bindings

# Install libigl bindings:
#    pip install pyg

import igl
import numpy as np

def main():
    # 1) Load your OFF mesh into NumPy arrays
    #    V: (n,3) float64, F: (m,3) int32
    V, F = igl.read_triangle_mesh("equilibrated.off")
    print(f"Original mesh: {V.shape[0]} vertices, {F.shape[0]} faces")

    # 2) Upsample (midpoint subdivision)
    n_subdivs = 1
    V_up, F_up = igl.upsample(V, F, n_subdivs)
    print(f"Upsampled mesh: {V_up.shape[0]} vertices, {F_up.shape[0]} faces")
    igl.write_triangle_mesh("equilibrated_upsampled.off", V_up, F_up)

    # # 3) Coarsen (quadric decimation) to ~1000 faces
    # target_faces = 1000
    # ratio = target_faces / F.shape[0]
    # # Estimate target vertex count proportionally
    # target_vertices = max(4, int(V.shape[0] * ratio))

    # # Ensure appropriate dtypes and memory layouts
    # Vf = np.array(V, dtype=np.float64, order='F')
    # Ff = np.array(F, dtype=np.int32,   order='F')

    # # decimate(V, F, max_m) → (V_coarse, F_coarse, J, I)
    # V_co, F_co, J, I = igl.decimate(Vf, Ff, target_vertices)
    # print(f"Coarsened mesh: {V_co.shape[0]} vertices, {F_co.shape[0]} faces")
    # igl.write_triangle_mesh("hauser_cube_coarse.off", V_co, F_co)

if __name__ == "__main__":
    main()
