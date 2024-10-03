import argparse
from ase.io import read, write
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
import numpy as np

def remove_surface_atoms(traj):
    """Remove Cu, Ag, Au atoms from the trajectory frames."""
    metals = ['Cu', 'Ag', 'Au']
    filtered_traj = []

    for atoms in traj:
        # Create a new atoms object with only non-metal atoms
        filtered_atoms = atoms[[atom.index for atom in atoms if atom.symbol not in metals]]
        filtered_traj.append(filtered_atoms)
    
    return filtered_traj

def select_anchors(traj_path, anchor_path, no_surface):
    # Load ASE trajectory (e.g., extxyz format)
    traj = read(traj_path, index=':')

    # If the --no_surface flag is present, remove metal atoms from the trajectory
    if no_surface:
        print("Removing metal surface atoms (Cu, Ag, Au)...")
        traj_for_pca = remove_surface_atoms(traj)
    else:
        traj_for_pca = traj

    # Extract Cartesian coordinates of all frames
    n_frames = len(traj_for_pca)
    n_atoms = len(traj_for_pca[0])
    coordinates = np.zeros((n_frames, n_atoms * 3))

    for i, atoms in enumerate(traj_for_pca):
        coordinates[i, :] = atoms.get_positions().flatten()

    # Apply PCA to the coordinates
    pca = PCA(n_components=6)  # Reduce dimensionality to 6 principal components
    reduced_coords = pca.fit_transform(coordinates)

    # Select the number of clusters (anchors)
    n_clusters = int(np.ceil(np.min([50, np.sqrt(len(traj_for_pca))])))
    print(f'Selecting {n_clusters} anchors')

    # Use KMeans to select representative frames
    kmeans = KMeans(n_clusters=n_clusters).fit(reduced_coords)

    # Get indices of the representative frames (closest to the cluster centers)
    representative_indices = []
    for cluster in range(n_clusters):
        indices = np.where(kmeans.labels_ == cluster)[0]
        centroid = kmeans.cluster_centers_[cluster]
        closest_frame = indices[np.argmin(np.linalg.norm(reduced_coords[indices] - centroid, axis=1))]
        representative_indices.append(closest_frame)

    # Output the representative frames (anchors)
    print("Representative frames (anchor points):", representative_indices)

    # Extract and save the anchor frames (with surface included)
    anchors = [traj[i] for i in representative_indices]
    write(anchor_path, anchors)

if __name__ == "__main__":
    # Argument parser setup
    parser = argparse.ArgumentParser(description="Select anchor frames from an ASE trajectory")
    parser.add_argument('traj_path', type=str, help='Path to the input trajectory file (e.g., extxyz)')
    parser.add_argument('anchor_path', type=str, help='Path to save the anchor frames (e.g., anchors.xyz)')
    parser.add_argument('--no_surface', action='store_true', help='Remove Cu, Ag, Au atoms for PCA but include them in output')

    # Parse arguments
    args = parser.parse_args()

    # Select anchors based on the given trajectory path and save to anchor_path
    select_anchors(args.traj_path, args.anchor_path, args.no_surface)

