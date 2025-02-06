import sys
import random
import numpy as np
import os
import subprocess
import math
import shutil
import time
import glob
#from mpi4py import MPI

#All input functions for extracting random resids
# Function to read a GRO file
def readGRO(input_gro_file):
    with open(input_gro_file, 'r') as f:
        lines = f.readlines()
    return lines[2:-1]  # Skip header lines and the last line

# Function defining the atomic masses of elements
def atomicMass(lbl):
    mass = {
        "H": 1.008,
        "C": 12.011,
        "N": 14.007,
        "O": 15.999,
        "S": 32.065,
    }
    return mass.get(lbl, 0)

#Function to calculate the COM of a molecule.
def getCOM(mol, vel=True):
    mr_tot = np.zeros(3)
    m_tot = 0.0
    for line in mol:
        l = line.split()
        m = atomicMass(l[1][0])
        if vel:
            coords = np.array([float(i) for i in l[-6:-3]])
        else:
            coords = np.array([float(i) for i in l[-3:]])
                        
        mr_tot += coords * m
        m_tot += m
    return mr_tot / m_tot

#Function to select molecules from the input GRO file based on COM and limits.
def selectMolec(lines, natoms,  x_min, x_max, y_limits, z_min, z_max):
    mol = []
    selected_resids = {}

    for line in lines:
        mol.append(line)

        if len(mol) == natoms:
            com = getCOM(mol, vel=True)
            #print(f"COM: {com}")
            x_com = com[0]  # Get the X-coordinate of the COM 
            y_com = com[1]  # Get the Y-coordinate of the COM
            z_com = com[2]  # Get the Z-coordinate of the COM

            if x_min <= x_com <= x_max and z_min <= z_com <= z_max:
                # check y limits
                if isinstance(y_limits, tuple):   # range 1
                    if y_limits[0] <= y_com <= y_limits[1]:
                        resid = ''.join(filter(str.isdigit, mol[0].split()[0]))
                        selected_resids[resid] = com
                        #print(f"Selected residue: {resid}")
                elif isinstance(y_limits, list): # range 2
                    for y_range in y_limits:
                        if y_range[0] <= y_com <= y_range[1]:
                            resid = ''.join(filter(str.isdigit, mol[0].split()[0]))
                            selected_resids[resid] = com
                            #print(f"Selected residue: {resid}")
                            break

            mol = []
    #print(f"Total selected residues: {len(selected_resids)}")
    return selected_resids

#Function to load resids from the input file.
def load_resids(traj_resid_file):
    if not os.path.exists(traj_resid_file):
        return set()
    with open(traj_resid_file, 'r') as f:
        return set(line.strip() for line in f)
    
# Function to save resids to the output file.
def save_resids(traj_resid_file, resids):
    with open(traj_resid_file, 'w') as f:
        for resid in resids:
            f.write(f"{resid}\n")
 
#Function to ask axes and limits from the user for extracting the input random resids
def ask_axes_and_limits(axes_count, axes, selection_choice, y_range_choice, num_to_select=None):
    limits = {"x_min": -np.inf, "x_max": np.inf, "z_min": -np.inf, "z_max": np.inf}
    predefined_limits = {
        "x": {"min": 2.0, "max": 14.0},  # Example limits for x-axis
        "z": {"min": 7.5, "max": 9.0}    # Example limits for z-axis
    }
    # Predefined y-limits ranges
    y_limits_range_1 = (4.0, 9.0)  # range 1
    y_limits_range_2 = [(4.5, 5.5), (9.0, 10.0)]  # range 2

    axes = axes.split(",")
    if len(axes) != axes_count:
        print(f"Error: You specified {axes_count} axes, but provided {len(axes)} axes.")
        sys.exit(1)

    for axis in axes:
        axis = axis.strip().lower()
        if axis == "y":
            if y_range_choice == 1:
                limits["y_limits"] = y_limits_range_1
            elif y_range_choice == 2:
                limits["y_limits"] = y_limits_range_2
            else:
                print("Invalid choice for y_range_choice. Please enter 1 or 2.")
                sys.exit(1)
        elif axis in predefined_limits:
            limits[f"{axis}_min"] = predefined_limits[axis]["min"]
            limits[f"{axis}_max"] = predefined_limits[axis]["max"]
        else:
            print(f"Invalid axis: {axis}")
            sys.exit(1)

    if selection_choice == "yes":
        limits["select_all"] = True
    elif selection_choice == "no":
        if num_to_select is None:
            print("Error: You must specify the number of residues to select when selection_choice is 'no'.")
            sys.exit(1)
        limits["select_all"] = False
        limits["num_to_select"] = num_to_select
    else:
        print("Invalid choice for selection. Please enter 'yes' or 'no'.")
        sys.exit(1)

    return limits

#Function to extract random resids based on the defined axes and limits.
def extract_and_update_resids(input_gro_file, traj_resid_file, natoms, limits):
    if os.path.exists(traj_resid_file):
        previous_resids = load_resids(traj_resid_file)
        os.rename(traj_resid_file, "old_random_resids.txt")

        lines = readGRO(input_gro_file)
        # Pass only the limits relevant to x, y, and z axes
        selected_resids = selectMolec(lines, natoms, limits["x_min"], limits["x_max"], limits["y_limits"], limits["z_min"], limits["z_max"])

        new_resids = list(set(selected_resids) - previous_resids)  # Find new residues
        # new_resids = list(set(selected_resids) - load_resids("old_random_resids.txt")) if os.path.exists("old_random_resids.txt") else list(selected_resids)
        if limits["select_all"]:
            selected_residues = new_resids  # Select all residues
        else:
            num_to_select = min(limits["num_to_select"], len(new_resids))  # Select the specified number
            selected_residues = random.sample(list(new_resids), num_to_select)  # Convert to list

        save_resids(traj_resid_file, selected_residues)
        print(f"Extracted {len(selected_residues)} residues and saved to {traj_resid_file}")

        # # Print the COM coordinates of the selected residues
        # print("\nCOM coordinates of selected residues:")
        # for resid in selected_residues:
        #     print(f"Resid {resid}: {selected_resids[resid]}")

    else:
        lines = readGRO(input_gro_file)
        # Pass only the limits relevant to x, y, and z axes
        selected_resids = selectMolec(lines, natoms, limits["x_min"], limits["x_max"], limits["y_limits"], limits["z_min"], limits["z_max"])
        
        if limits["select_all"]:
            selected_residues = selected_resids  # Select all residues
        else:
            num_to_select = min(limits["num_to_select"], len(selected_resids))  # Select the specified number
            selected_residues = random.sample(list(selected_resids), num_to_select)  # Convert to list

        save_resids(traj_resid_file, selected_residues)
        print(f"Extracted {len(selected_residues)} residues and saved to {traj_resid_file}")

        # # Print the COM coordinates of the selected residues
        # print("\nCOM coordinates of selected residues:")
        # for resid in selected_residues:
        #     print(f"Resid {resid}: {selected_resids[resid]}")
        
def extract_random_resids(input_gro_file, traj_resid_file, natoms, axes_count, axes, selection_choice, y_range_choice, num_to_select=None):
    limits = ask_axes_and_limits(axes_count, axes, selection_choice, y_range_choice, num_to_select)
    extract_and_update_resids(input_gro_file, traj_resid_file, natoms, limits)
    
#All functions for nearest neighbors calculation using ASANN algorithm
# Atomic weights for COM calculations
ATOMIC_WEIGHTS = {
    'C': 12.011,
    'H': 1.008
}
#Function to get residues based on the source_resid
def get_residues(input_gro_file):
    """
    Reads the GRO file and returns atom coordinates and residue information.
    """
    residues = {}
    with open(input_gro_file, 'r') as f:
        lines = f.readlines()[2:-1]  # Skip header and footer
        for line in lines:
            resid = int(line[0:5].strip())
            atom_name = line[10:15].strip()
            x, y, z = map(float, (line[20:28].strip(), line[28:36].strip(), line[36:44].strip()))
            if resid not in residues:
                residues[resid] = []
            residues[resid].append((atom_name, np.array([x, y, z])))
    return residues

#function to get COM for that particular residues
def compute_com(residues):
    """
    Computes the center of mass (COM) for each molecule.
    """
    coms = {}
    for resid, atoms in residues.items():
        total_mass = 0.0
        weighted_positions = np.zeros(3)
        for atom_name, position in atoms:
            # Get atomic weight (default to 0 if not found)
            element = atom_name[0]  # Assume the first letter indicates the element
            mass = ATOMIC_WEIGHTS.get(element, 0.0)
            total_mass += mass
            weighted_positions += mass * position
        com = weighted_positions / total_mass
        coms[resid] = com
    return coms

#Function to get PBC vectors
def get_pbc_vectors(coords, nb_atoms=None):
    """
    Compute pairwise vectors without periodic boundaries conditions.
    
    Parameters:
        coords (numpy 2D array): List of atom coordinates.
        pbc (bool): Ignored; function assumes no PBC.
        nb_atoms (int, optional): Number of real atom coordinates.
    
    Returns:
        numpy 2D array of pairwise vectors (cartesian coordinates).
    """
    # Retrieve number of real atoms
    if nb_atoms is None:
        nb_atoms = coords.shape[0]
    
    # Compute pairwise vectors
    vectors = coords[np.newaxis, :, :] - coords[:nb_atoms, np.newaxis, :]
    return (vectors)

#Function to get sorted distances for SANN
def get_sorted_distances(coords, nb_atoms=None):
    """
    Compute sorted pairwise distances and vectors.
    """
     # Retrieve number of atoms if not given
    if nb_atoms is None:
        nb_atoms = coords.shape[0]
    
    # Compute pairwise vectors
    vectors = get_pbc_vectors(coords, nb_atoms=nb_atoms)  # No PBC applied
    
    # Compute pairwise distances
    distances = np.sqrt(np.sum(vectors**2, axis=-1))
    
    # Sort distances and vectors
    sorted_index_axis1 = np.argsort(distances, axis=-1)
    sorted_index_axis0 = np.arange(nb_atoms)[:, None]
    distances = distances[sorted_index_axis0, sorted_index_axis1]
    vectors = vectors[sorted_index_axis0, sorted_index_axis1]
    
    return distances, vectors, sorted_index_axis1

#Function to get SANN outputs
def get_SANN(all_distances):
    """
    Computes coordination numbers using the SANN algorithm.
    """
    # Retrieve number of atoms
    nb_coords = all_distances.shape[1]
 
    list_sann_CN = []
    list_sann_radius = []
    list_dist_sum = all_distances[:, 1:4].sum(axis=1)
    for (dist_sum, atom_distances) in zip(list_dist_sum, all_distances):
        sann_CN = 3
        while (sann_CN + 1 < nb_coords) and (dist_sum / (sann_CN - 2) >= atom_distances[sann_CN + 1]):
            dist_sum += atom_distances[sann_CN + 1]
            sann_CN += 1
        list_sann_CN.append(sann_CN)
        list_sann_radius.append(dist_sum / (sann_CN - 2))
    return (np.array(list_sann_CN), np.array(list_sann_radius))

def dist_to_barycenter(nearest_neighbors, nearest_distances, radius):
    # Compute solid angles of neighbors
    list_SA = 1 - nearest_distances / radius

    # Compute solid-angle-weighted barycenter
    bary_vector = np.sum(nearest_neighbors * list_SA[:, np.newaxis], axis=0) / np.sum(list_SA)

    # Return distance from the central atom to the barycenter
    return(math.sqrt(np.sum(bary_vector**2)))

def angular_correction(nearest_neighbors, nearest_distances, radius):
    """
    Computes the angular correction for ASANN.
    """
    alpha = dist_to_barycenter(nearest_neighbors, nearest_distances, radius) / radius
    return ((alpha + math.sqrt(alpha**2 + 3 * alpha)) / 3)

def get_ASANN(sorted_distances, sorted_vectors, sann_CNs, sann_radii):
    """
    Computes ASANN-based coordination numbers.
    """
    # Retrieve number of atoms
    nb_coords = sorted_distances.shape[1]
    
    list_asann_CN = []
    list_asann_radius = []
    for (atom_distances, atom_neighbours, sann_CN, sann_radius) in zip(sorted_distances, sorted_vectors, sann_CNs, sann_radii):
        # computes angular correction
        nearest_distances = atom_distances[1:sann_CN+1]
        nearest_neighbours = atom_neighbours[1:sann_CN+1]
        ang_corr = angular_correction(nearest_neighbours, nearest_distances, sann_radius)
        beta = 2*(1 - ang_corr)
        
        # ASANN algorithm (i.e. while ASANN radius sum(r_ij, j=1..m)/(m-2*(1-ang_corr)) >= r_i(m+1), increase m by 1)
        asann_CN = int(beta) + 1 # set CN to floor(2*(1-ang_corr)) + 1 (i.e. the minimum CN value for the ASANN algorithm)
        dist_sum = atom_distances[1:asann_CN+1].sum()
        while (asann_CN + 1 < nb_coords) and (dist_sum / (asann_CN - beta) >= atom_distances[asann_CN + 1]):
            dist_sum += atom_distances[asann_CN + 1]
            asann_CN += 1
        # store the ASANN CN and radius found
        list_asann_CN.append(asann_CN)
        list_asann_radius.append(dist_sum / (asann_CN - beta))
    return (np.array(list_asann_CN), np.array(list_asann_radius))

def nearest_neighbors_asann(input_gro_file, source_resid):
    residues = get_residues(input_gro_file)
    coms = compute_com(residues)
    coords = np.array(list(coms.values()))
    resid_indices = np.array(list(coms.keys()))
    if source_resid not in resid_indices:
        print(f"Residue ID {source_resid} not found in the file.")
        return {}
    sorted_distances, sorted_vectors, sorted_idx = get_sorted_distances(coords)
    sann_CNs, sann_radii = get_SANN(sorted_distances)
    asann_CNs, asann_radii = get_ASANN(sorted_distances, sorted_vectors, sann_CNs, sann_radii)
    source_index = np.where(resid_indices == source_resid)[0][0]
    sann_neighbors = sorted_idx[source_index, 1:sann_CNs[source_index] + 1] + 1
    asann_neighbors = sorted_idx[source_index, 1:asann_CNs[source_index] + 1] + 1
    extracted_nn = {
        'source_resid': source_resid,
        #'sann_neighbors': sann_neighbors.tolist(),
        'asann_neighbors': asann_neighbors.tolist()
    }
    
    return extracted_nn

# All functions for calculating average coupling values
# Function to generate Charge transfer.dat file
def generate_charge_transfer_dat(source_resid, target_site, spec_file):
    charge_transfer_dat = f"""seed = 1
hamiltonian = dftb
tfermi = 300
slkopath = /home/ka/ka_ipc/ka_rv4081/Sonali/dftb-parameters/slko-files-mio/
wavefunctionreal = 1.0 0.0
chargecarrier = hole
offdiagscaling = yes
extchrmode = vacuo
espscaling = 1
nsitetypes = 1
sic = 0.0
typefiles = {spec_file}
nsites = 2
zonesize = 2
sites = {source_resid} {target_site}
atomindex = 1 20 29
sitetypes = 1 1
foshift = 0.0 0.0
sitescc = 0 0
jobtype = NOM
deltaqmode = mulliken
internalrelax = onsite"""
    return charge_transfer_dat

# Function to calculate average coupling from TB_HAMILTONIAN.xvg
def calculate_average_coupling(directory, hamiltonian_file):
    hamiltonian_file_path = os.path.join(directory, hamiltonian_file)
    if os.path.exists(hamiltonian_file_path):
        with open(hamiltonian_file_path, 'r') as f:
            lines = f.readlines()
            coupling_values = [float(line.split()[2]) for line in lines if line.strip()]
            avg_coupling = sum(coupling_values) / len(coupling_values) if coupling_values else 0.0
            return abs(avg_coupling)
    else:
        print(f"TB_HAMILTONIAN.xvg file not found in directory: {directory}")
        return None

# Function to check if a job is finished by looking for 'Finished' in ham.log
def is_job_finished(directory, ham_log_file):
    log_file = os.path.join(directory, ham_log_file)
    if os.path.exists(log_file):
        with open(log_file, 'r') as f:
            return 'Finished' in f.read()
    return False

# Function to clean up the directory, keeping only the specified file
def cleanup_directory(directory, file_to_keep):
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
        if filename != file_to_keep and os.path.isfile(file_path):
            os.remove(file_path)

def average_coupling_main(source_resid, input_gro_file, charge_transfer_submission_script, spec_file, mdp_file, top_file, hamiltonian_file, ham_log_file, sub_dir):
    extracted_coupling_values = {}
    dict_key = f"average_coupling_values_{source_resid}"
    extracted_coupling_values[dict_key] = {}
    
    # Extract nearest_neighbors around the source site using asann algorithm
    extracted_nn = nearest_neighbors_asann(input_gro_file, source_resid)
    print(extracted_nn)
    output_resids = extracted_nn.get('asann_neighbors', set())
    
    # Loop through each unique residue and calculate coupling
    for target_site in output_resids:
        print(f"Calculating coupling between sites {source_resid} and {target_site}...")

        # Create a directory for the pair
        pair_directory = os.path.join(sub_dir, f"{source_resid}_{target_site}")
        os.makedirs(pair_directory, exist_ok=True)
  
        # Generate charge_transfer.dat for each pair
        charge_transfer_dat = generate_charge_transfer_dat(source_resid, target_site, spec_file)
        charge_transfer_dat_file = os.path.join(pair_directory, "charge-transfer.dat")
        with open(charge_transfer_dat_file, 'w') as f:
            f.write(charge_transfer_dat)
    
        # Copy necessary files to the pair directory
        shutil.copy(charge_transfer_submission_script, pair_directory)
        shutil.copy(spec_file, pair_directory)
        shutil.copy(mdp_file, pair_directory)
    
        # Run the GROMACS command to generate ham.tpr in the current directory
        print("Running GROMACS grompp...")
        subprocess.run(["/lustre/home/ka/ka_ipc/ka_ki1980/gromacs-sh-old_0/build/src/kernel/grompp", "-f", mdp_file, "-c", input_gro_file, "-p", top_file, "-o", 'ham.tpr', "-maxwarn", '1'], cwd=pair_directory)
        
        # Check if the job has already finished
        if is_job_finished(pair_directory, ham_log_file):
            print(f"Job already completed for pair {source_resid} and {target_site}. Skipping...")
        else:
            # Run the charge transfer calculation script
            print("Running charge transfer calculation...")
            subprocess.run(['sbatch', charge_transfer_submission_script], cwd=pair_directory)
            print(f"Calculation initiated for pair {source_resid} and {target_site} in directory: {pair_directory}.")
    
        # Wait for the job to finish
        while not is_job_finished(pair_directory, ham_log_file):
            print(f"Job still running for pair {source_resid} and {target_site}. Checking again in 10 seconds...")
            time.sleep(10)
    
        print(f"Job done for pair {source_resid} and {target_site}.")
    
        # Calculate average coupling for the pair and store it dynamically
        avg_coupling = calculate_average_coupling(pair_directory, hamiltonian_file)
        if avg_coupling is not None:
            extracted_coupling_values[dict_key][target_site] = avg_coupling * 1000
        else:
            extracted_coupling_values[dict_key][target_site] = None
        
        # Remove test.sh to ensure it's not run again in subsequent iterations
        os.remove(os.path.join(pair_directory, charge_transfer_submission_script))
    
        # Clean up the directory, keeping only TB_HAMILTONIAN.xvg
        cleanup_directory(pair_directory, hamiltonian_file)
    
    print("All calculations completed.")
    return extracted_coupling_values

# Extract residue lines from .gro file
def extractResidueLines(resid, input_gro_file):
    residue_lines = []
    with open(input_gro_file, "r") as f:
        lines = f.readlines()
        for line in lines[2:-1]:
            try:
                current_resid = int(line[:5].strip())
                if current_resid == resid:
                    residue_lines.append(line)
            except ValueError:
                continue
    return residue_lines

# Calculate COM for a specific residue
def calculateCOM(resid, input_gro_file):
    residue_lines = extractResidueLines(resid, input_gro_file)
    if residue_lines:
        return getCOM(residue_lines)
    return None

# Function to read source_resid and target_resid from a file
def read_resid_data(average_coupling_values):
    source_resids = []
    target_resids = []
    coupling_data = []

    for source_resid_str, targets in average_coupling_values.items():
        source_resid = int(source_resid_str.replace("average_coupling_values_", ""))
        for target_resid, avg_cpl in targets.items():
                if avg_cpl is not None:
                    source_resids.append(source_resid)
                    target_resids.append(target_resid)
                    coupling_data.append(avg_cpl)
    
    return source_resids, target_resids, coupling_data

#Function to calculate probabilities
# def calculate_probabilities(target_resids, coupling_data):
#     probabilities = {}
#     weights = {}
#     total_weight = 0.0
#     kT = 25.7  # Boltzmann constant times temperature (298 K)

#     # Compute weights for each neighbor
#     for target_resid, coupling in zip(target_resids, coupling_data):
#         weight = math.exp(-coupling / kT)
#         weights[target_resid] = weight
#         total_weight += weight
#     else:
#         print(f"Warning: No dot product found for neighbor {target_resid}")

#     # Compute probabilities
#     for target_resid, weight in weights.items():
#         probabilities[target_resid] = 1- (weight / total_weight) if total_weight > 0 else 0.0

#     return probabilities

# probabilities calculated using softmax activation function
def calculate_probabilities(target_resids, coupling_data, beta=(1/25.7)):
    probabilities = {}
    exp_weight = {}
    total_weight = 0.0

    # Compute softmax weights for each neighbor
    max_coupling = max(coupling_data)  # Normalize to avoid overflow issues with large numbers
    exp_weights = []

    for coupling in coupling_data:
        exp_weight = math.exp(beta * (coupling - max_coupling))  # Shifted for numerical stability
        exp_weights.append(exp_weight)
        total_weight += exp_weight

    # Calculate the probabilities for each target_resid
    for target_resid, exp_weight in zip(target_resids, exp_weights):
        probabilities[target_resid] = exp_weight / total_weight if total_weight > 0 else 0.0

    return probabilities

def check_too_far(source_com, selected_neighbor_com):
    """Checks if the neighbor is too far from the source."""
    distance = math.sqrt((selected_neighbor_com[0] - source_com[0])**2 + (selected_neighbor_com[1] - source_com[1])**2 + (selected_neighbor_com[2] - source_com[2])**2)
    print(f"Distance between source_resid and neighbor_resid: {distance:.2f} nm")
    return distance > 0.70

#Function to check if the selected neighbor has already been visited
def check_visited(selected_neighbor, visited_resids):
    """Checks if the neighbor has already been visited."""
    return selected_neighbor in visited_resids

def check_range(resid, input_gro_file):
    """
    Checks if a given resid's center of mass (y_COM) is within predefined ranges.
    Logs the result and returns the range classification.
    """
    com = calculateCOM(resid, input_gro_file)
    if com is not None:
        y_COM = com[1]  # Extract y-coordinate
        
        # Define the ranges
        y0, y1 = 4.5, 5.5  # First range
        y2, y3 = 9.0, 10.0  # Second range
        
        # Check if y_COM is in either range
        if y0 <= y_COM <= y1:
            print(f"Resid {resid} satisfies y_COM {y_COM} in range ({y0}, {y1})")
            return "in_range_1"
        elif y2 <= y_COM <= y3:
            print(f"Resid {resid} satisfies y_COM {y_COM} in range ({y2}, {y3})")
            return "in_range_2"
        else:
            print(f"Resid {resid} does NOT satisfy y_COM {y_COM} in any range")
            return "out_of_range"
    return None

#original select_neighbor function without fallback mechanism
# def select_neighbor(probabilities, source_resid, input_gro_file, visited_resids):
#     while True: #Loop allows for backtracking if needed
#         random_value = random.uniform(0, 1)  # Generate a random  number between 0 and 1
#         print(f"Random number for source_resid {source_resid}: {random_value}")
    
#         # Compute source center of mass
#         source_com = calculateCOM(source_resid, input_gro_file)
#         source_range = check_range(source_resid, input_gro_file)
    
#         # Sort neighbors by probability closeness to `random_value`
#         sorted_neighbors = sorted(probabilities.items(), key=lambda x: abs(x[1] - random_value))
    
#         for target_resid, _ in sorted_neighbors:
#             if target_resid in visited_resids:
#                 print(f"Neighbor {target_resid} has already been visited. Skipping...")
#                 continue
    
#             # Compute neighbor COM
#             target_com = calculateCOM(target_resid, input_gro_file)
    
#             # Check movement in the correct y-direction based on first source_resid range
#             if source_range == "in_range_1" and target_com[1] <= source_com[1]:
#                 print(f"Neighbor {target_resid} does not move forward in range 1. Skipping...")
#                 continue
#             elif source_range == "in_range_2" and target_com[1] >= source_com[1]:
#                 print(f"Neighbor {target_resid} does not move forward in range 2. Skipping...")
#                 continue
    
#             # Check if the selected neighbor is too far
#             if check_too_far(source_com, target_com):
#                 print(f"Neighbor {target_resid} is too far. Skipping...")
#                 continue
    
#             # Assign selected neighbor if valid
#             selected_neighbor = target_resid
#             print(f"Selected valid neighbor: {target_resid} (probability: {probabilities[target_resid]})")
#             return selected_neighbor  # Return the first valid neighbor found
    
#         # **Backtracking if no valid neighbor found**
#         if visited_resids:
#             last_visited = visited_resids.pop()  # Remove and get the last visited residue
#             print(f"No valid neighbor found. Backtracking to {last_visited}...")
#             source_resid = last_visited
#             return select_neighbor(probabilities, source_resid, input_gro_file, visited_resids)
    
#         else:
#             # If no previous residues left, return None
#             print("No previous residues to backtrack to. Stopping.")
#             return None


#modified select_neighbor function according to fallback mechanism
def select_neighbor(probabilities, source_resid, first_source_resid, input_gro_file, visited_resids, sampled_paths, charge_transfer_submission_script, spec_file, mdp_file, top_file, hamiltonian_file, ham_log_file, sub_dir, n=1):
    """
    Selects a valid neighbor for the given source_resid based on probabilities.
    if no valid neighbor is found, backtracks and calls fallback_select_neighbor.
    """
    while True: #Loop allows for backtracking if needed
        random_value = random.uniform(0, 1)  # Generate a random  number between 0 and 1
        print(f"Random number for source_resid {source_resid}: {random_value}")
    
        # Compute source center of mass
        source_com = calculateCOM(source_resid, input_gro_file)
        source_range = check_range(first_source_resid, input_gro_file)
    
        # Sort neighbors by probability closeness to `random_value`
        sorted_neighbors = sorted(probabilities.items(), key=lambda x: abs(x[1] - random_value))
    
        for target_resid, _ in sorted_neighbors:
            if target_resid in visited_resids:
                print(f"Neighbor {target_resid} has already been visited. Skipping...")
                continue
    
            # Compute neighbor COM
            target_com = calculateCOM(target_resid, input_gro_file)
    
            # Check movement in the correct y-direction based on first source_resid range
            if source_range == "in_range_1" and target_com[1] <= source_com[1]:
                print(f"Neighbor {target_resid} does not move forward in range 1. Skipping...")
                continue
            elif source_range == "in_range_2" and target_com[1] >= source_com[1]:
                print(f"Neighbor {target_resid} does not move forward in range 2. Skipping...")
                continue
    
            # Check if the selected neighbor is too far
            if check_too_far(source_com, target_com):
                print(f"Neighbor {target_resid} is too far. Skipping...")
                continue
    
            # Assign selected neighbor if valid
            selected_neighbor = target_resid
            print(f"Selected valid neighbor: {target_resid} (probability: {probabilities[target_resid]})")
            return selected_neighbor  # Return the first valid neighbor found
    
        # **Backtracking if no valid neighbor found**
        if visited_resids:
            sampled_paths.remove(source_resid)  # Remove the source_resid from the sampled paths
            last_visited = visited_resids[-n]  # Remove and get the last visited residue
            print(f"No valid neighbor found. Backtracking to {last_visited}...")

            #call the fallback_select_neighbor function
            fallback_resid = last_visited
            selected_neighbor = fallback_select_neighbor(fallback_resid, input_gro_file, charge_transfer_submission_script, spec_file, mdp_file, top_file, hamiltonian_file, ham_log_file, sub_dir, visited_resids, sampled_paths, first_source_resid)
            return selected_neighbor            
    
        else:
            # If no previous residues left, return None
            print("No previous residues to backtrack to. Stopping.")
            return None

# functio defining the fallback mechanism
def fallback_select_neighbor(fallback_resid, input_gro_file, charge_transfer_submission_script, spec_file, mdp_file, top_file, hamiltonian_file, ham_log_file, sub_dir, visited_resids, sampled_paths, first_source_resid):
    """
    A fallback mechanism to find a neighboring residue if the primary attempt fails.
    """
    print(f"Fallback mechanism activated for fallback_resid {fallback_resid}...")
    
    #Compute avearge coupling values for the fallback_resid
    average_coupling_values = average_coupling_main(fallback_resid, input_gro_file, charge_transfer_submission_script, spec_file, mdp_file, top_file, hamiltonian_file, ham_log_file, sub_dir, visited_resids, sampled_paths, first_source_resid)
    print(f"Average coupling values for fallback_resid {fallback_resid}: {average_coupling_values}")

    # read and process the coupling data
    _, target_resids, coupling_data = read_resid_data(average_coupling_values)

    # calculate probabilities for the fallback_resid
    probabilities = calculate_probabilities(target_resids, coupling_data)
    print(f"Final probabilities for fallback_resid {fallback_resid}: {probabilities}")

    #call select_neighbor with the fallback resid
    selected_neighbor = select_neighbor(probabilities, fallback_resid, first_source_resid, input_gro_file, visited_resids, sampled_paths, charge_transfer_submission_script, spec_file, mdp_file, top_file, hamiltonian_file, ham_log_file, sub_dir, n=1)

    if selected_neighbor is None:
        print(f"No selected neighbor found for fallback_resid {fallback_resid}. Exiting...")

    return selected_neighbor

# main function updated and testing
def main():
    #MPI setup
    #comm = MPI.COMM_WORLD
    #rank = comm.Get_rank()  # Process ID
    #size = comm.Get_size()  # Total numbe rof processes

    #set parameters
    traj_path = "/lustre/work/ws/ws1/ka_rv4081-sfb_cm/pentacene/gb_pen-schellhammer/y45/transfer_integral_distr/snapshots_asann/input_files/trajectories"
    traj_file_prefix = "TRAJSS"
    output_extract_resid_file = "random_resids.txt"
    natoms = 36
    axes_count = 3
    axes = "x,y,z"
    selection_choice = "no"
    y_range_choice = 2
    num_to_select = 1
    mdp_file = "/lustre/work/ws/ws1/ka_rv4081-sfb_cm/pentacene/gb_pen-schellhammer/y45/parallel_workflow_test/input_files/namd-qmmm.mdp"
    top_file = "/lustre/work/ws/ws1/ka_rv4081-sfb_cm/pentacene/gb_pen-schellhammer/y45/parallel_workflow_test/input_files/pen-esp.top"
    charge_transfer_submission_script = "test.sh"
    spec_file = "pen-hole.spec"
    hamiltonian_file = "TB_HAMILTONIAN.xvg"
    ham_log_file = "ham.log"

    # define trajectory range
    itraj = 1
    ftraj = 1

    #Distribute workloads among MPI ranks
    #all_trajs = list(range(itraj, ftraj+1))
    #trajs_per_rank = [all_trajs[i::size] for i in range(size)]  #split across ranks
    #my_trajs = trajs_per_rank[rank]  #get the trajs for this rank

    #for traj in my_trajs:
    for traj in range(itraj, ftraj+1):
        traj_dir = f"TRAJ{traj}"
        os.makedirs(traj_dir, exist_ok=True)
        input_gro_file = os.path.join(traj_path, f"{traj_file_prefix}{traj}.gro")
        print(f"Processing trajectory: {traj}, input file: {input_gro_file}")
        
        if not os.path.exists(input_gro_file):
            print(f"Error: File {input_gro_file} does not exist. Skipping this trajectory.")
            break

        #step 1: calling extract_random_resids function to extract random start resids within the specified parameters
        traj_resid_file = os.path.join(traj_dir, output_extract_resid_file)
        extract_random_resids(input_gro_file, traj_resid_file, natoms, axes_count, axes, selection_choice, y_range_choice, num_to_select)

        # step 2: Looping over the source_resid present in random_resids.txt
        with open(traj_resid_file, 'r') as f:
            source_resids = [line.strip() for line in f.readlines()]
            
        for source_resid in source_resids:
            source_resid = int(source_resid)
            
            # Create a subdirectory for each source_resid inside respective trajectory directory
            sub_dir = os.path.join(traj_dir, f"subdir_{source_resid}")
            os.makedirs(sub_dir, exist_ok=True)
            print(f"Processing source_resid: {source_resid}, created subdir: {sub_dir}")
    
            #Assign first source_resid to the current source_resid
            first_source_resid = source_resid
            sampled_paths = [first_source_resid]
            print(f"First source_resid: {first_source_resid}")

            visited_resids = []
            visited_resids.append(first_source_resid)
              
            first_range = check_range(first_source_resid, input_gro_file)
            if first_range == "out_of_range":
                print(f"Source_resid {first_source_resid} is out of range. Skipping this source_resid.")
                continue

            while True:
                average_coupling_values = average_coupling_main(source_resid, input_gro_file, charge_transfer_submission_script, spec_file, mdp_file, top_file, hamiltonian_file, ham_log_file, sub_dir)
                print(f"Average coupling values for source_resid {source_resid}: {average_coupling_values}")

                source_resids, target_resids, coupling_data = read_resid_data(average_coupling_values)
                
                probabilities = calculate_probabilities(target_resids, coupling_data)
                print(f"Final probabilities for source_resid {source_resid}: {probabilities}")
                
                selected_neighbor = select_neighbor(probabilities, source_resid, first_source_resid, input_gro_file, visited_resids, sampled_paths, charge_transfer_submission_script, spec_file, mdp_file, top_file, hamiltonian_file, ham_log_file, sub_dir, n=1)
                print(f"Selected neighbor for source_resid {source_resid}: {selected_neighbor}")
                        
                if selected_neighbor is None:
                    print(f"No selected neighbor found for source_resid {source_resid}. Exiting...")
                    break

                visited_resids.append(selected_neighbor)

                # Add selected_neighbor to the sampled paths list
                sampled_paths.append(selected_neighbor)

                # Check the range of selected_neighbor
                selected_range = check_range(selected_neighbor, input_gro_file)

                # stop if selected_neighbor is in opposite range
                if (first_range == "in_range_1" and selected_range == "in_range_2") or \
                    (first_range == "in_range_2" and selected_range == "in_range_1"):
                    print(f"Selected neighbor {selected_neighbor} is in opposite range. Exiting...")

                    # Identify the current TRAJ directory from sub_dir
                    traj_dir = os.path.dirname(sub_dir)  # Get the parent TRAJ* directory

                    # Remove all subdir_* inside TRAJ* directories
                    sub_dirs = glob.glob(os.path.join(traj_dir, "subdir_*"))
                    for sub_dir in sub_dirs:
                        try:
                            os.system(f"rm -r {sub_dir}")  # Remove the subdir safely
                        except Exception as e:
                            print(f"Error removing {sub_dir}: {e}")

                     # Write final sampled_paths to a text file
                     output_file = os.path.join(traj_dir, "final_sampled_paths.txt")
                     with open(output_file, "w") as f:
                          for resid in sampled_paths:
                              f.write(f"{resid}\n")

                     print(f"Final sampled paths saved to {output_file}")
                     break    
                
                # Continue iteration with selected_neighbor as the new source_resid
                print(f"Continuing with selected_neighbor {selected_neighbor} as the new source_resid.")
                source_resid = selected_neighbor
                
            print(f"Sampled paths for source_resid {first_source_resid}: {sampled_paths}")

if __name__ == "__main__":
    main()
