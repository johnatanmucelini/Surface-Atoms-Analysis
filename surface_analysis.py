import sys
import numpy as np
import ase.io
import argparse
from scipy.stats import describe
from scipy.spatial.distance import cdist
from sklearn.cluster import DBSCAN

def RegRDS_set(sampling_distance, N):
    """Return a set of N R3 dots (almost) regular  distributed in the surface of
    a sphere of radius 'sampling_distance'.
    More deatiail of the implementation in the article "How to generate
    equidistributed points on the surface of a sphere" by Markus Deserno:
    https://www.cmu.edu/biolphys/deserno/pdf/sphere_equi.pdf
    samplind_distance: a float or a int grater than zero.
    N: intiger grater than zero."""

    if not sampling_distance > 0:
        sys.exit("sampling_distance must be higher than zero! Aborting...")

    if (type(N) != int) or (N <= 0):
        sys.exit("N must be an intiger grater than zero! Aborting...")

    cart_coordinates=[]
    r = 1
    Ncount = 0
    a = 4. * np.pi * r**2 / N
    d = np.sqrt(a)
    Mtheta = int(round(np.pi/d))
    dtheta = np.pi / Mtheta
    dphi = a / dtheta
    for m in range(0, Mtheta ) :
        theta = np.pi *( m + 0.5 ) / Mtheta
        Mphi = int(round(2 *np.pi * np.sin(theta) / dphi ))
        for n in range( 0 , Mphi ) :
            phi = 2* np.pi * n / Mphi
            Ncount += 1
            y = sampling_distance * np.sin(theta) * np.cos(phi)
            x = sampling_distance * np.sin(theta) * np.sin(phi)
            z = sampling_distance * np.cos(theta)
            cart_coordinates.append([x,y,z])
    cart_coordinates = np.array(cart_coordinates)

    return cart_coordinates

def writing_points_xyz(file_name, positions):
    """Write positions, a list or array of R3 points, in a xyz file file_named.
    Several softwares open xyz files, such as Avogadro and VESTA
    In the xyz file all the atoms are H.
    file_name: a string with the path of the xyz document which will be writed.
    positions: a list or numpy array with the atoms positions."""

    if type(file_name) != str :
        sys.exit("file_name must be a string")

    for index, element in enumerate(positions) :
        if len(element) != 3 :
            sys.exit("Element " + str(index) + " of positions does not present three elements." )

    if type(positions) != list : positions = np.array(positions)

    ase.io.write( file_name , ase.Atoms('H'+str(len(positions)) , list(map(tuple,positions)) ))

def large_surfaces_index(surface_dots_positions, eps):
    """Seach if there are more than one surfaces in surface_dots_positions, than
    return a np.array of booleans with True for index of atoms for the surfaces
    more points.
    surface_dots_positions: np.array with R3 dots.
    eps: minimun distance between different surfaces, should be a float grater
         than zero."""

    if (not isinstance(surface_dots_positions, np.ndarray)
            or (np.shape(surface_dots_positions)[1] != 3)):
        sys.exit("surface_dots_positions must be a (n,3) shaped np.array. Aborting...")

    if not eps > 0:
        sys.exit("eps must be large than zero.")

    db = DBSCAN(eps=eps, min_samples=1).fit_predict(surface_dots_positions)
    labels, quanity = np.unique(db, return_counts=True)
    if len(labels) > 1:
        print('{} surfaces were found (sizes: {}), the bigger were selected as the external!'.format(
              labels, str(quanity).replace('[', '').replace(']', '')))
    result = db == labels[np.argmax(quanity)]

    return result

def surfcore_classyfier_adatom(positions, atomic_radii, adatom_radius,
                               remove_is=True, ssamples=1000, sp_file="surface_points.xyz"):
    """This algorithm classify atoms in surface and core atoms employing the
    concept of atoms as ridge spheres. Then the surface atoms are the ones that
    could be touched by an fictitious adatom that approach the cluster, while
    the core atoms are the remaining atoms.
    See more of my algorithms im GitHub page Johnatan.mucelini.
    Articles which employed thi analysis: Mendes P. XXX
    .
    Parameters
    ----------
    positions: numpy array of floats (n,3) shaped.
               Cartezian positions of the atoms, in angstroms.
    atomic_radii: numpy array of floats (n,) shaped.
                  Radius of the atoms, in the same order which they appear in
                  positions, in angstroms.
    adatom_radius: float (optional, default=1.1).
                   Radius of the dummy adatom, in angstroms.
    ssampling: intiger (optional, default=1000).
               Quantity of samplings over the touched sphere surface of each
               atom.
    write_sp: boolean (optional, default=True).
              Define if the xyz positions of the surface points will
              be writed in a xyz file (surface_points.xyz).
    sp_file : string (optional, default='surface_points.xyz').
              The name of the xyz file to write the surface points positions,
              in angstroms, if write_sp == True.
    Return
    ------
    surface_exposition: numpy array of floats (n,).
                        The percentual of surface exposition of each atom.
    Example
    ------
    >>> ...
    """

    print('Surface analysis:')

    # Centralizing atoms positions:
    positions = positions - np.average(positions, axis=0)
    touch_radii = atomic_radii + adatom_radius
    qtna = len(positions)

    dots_try = RegRDS_set(adatom_radius, ssamples)
    rssamples = len(dots_try)
    max_dots = rssamples * qtna
    print('    Number of dots per atom: {}'.format(len(dots_try)))
    print('    Number of investigated dots: {}'.format(max_dots))

    dots = positions[0] + RegRDS_set(touch_radii[0] + 0.001, ssamples)
    dot_origin = [[0] * len(dots_try)]
    for atom_id in range(1, qtna):
        dots = np.append(dots, positions[atom_id] + RegRDS_set(touch_radii[atom_id] + 0.001, ssamples), axis=0)
        dot_origin.append([atom_id] * len(dots_try))
    dot_origin = np.array(dot_origin).flatten()

    # removing dots inside other touch sphere
    dots_atoms_distance = cdist(positions, dots)
    atomos_radii_projected = np.array([touch_radii]*len(dots)).reshape(len(dots), qtna).T
    surface_dots = np.sum(dots_atoms_distance < atomos_radii_projected, axis=0) == 0
    dots_surface = dots[surface_dots]

    # removing internal surfaces
    if remove_is:
        print("    Removing_internal_surface")
        dots_for_find_eps = RegRDS_set(max(touch_radii), ssamples)
        dotdot_distances = cdist(dots_for_find_eps, dots_for_find_eps)
        min_dotdot_dist = np.min(dotdot_distances + np.eye(rssamples)*10, axis=0)
        eps = 2.1 * np.max(min_dotdot_dist)
        external_surface_dots = large_surfaces_index(dots_surface, eps)
        dots_surface = dots_surface[external_surface_dots]
        surface_dot_origin = dot_origin[surface_dots][external_surface_dots]
    else:
        surface_dot_origin = dot_origin[surface_dots]
    surface_atoms, dots_per_atom = np.unique(surface_dot_origin, return_counts=True)

    # Getting exposition
    atoms_area_per_dot = (4 * np.pi * atomic_radii**2)/(1.*rssamples)
    exposition = np.zeros(qtna)
    is_surface = np.zeros(qtna, dtype=bool)
    incidence = np.zeros(qtna)
    for atom_id, atom_incidence in zip(surface_atoms, dots_per_atom):
        #print(' found surface atom: ' + str(atom_id))
        is_surface[atom_id] = True
        incidence[atom_id] = atom_incidence / rssamples
        exposition[atom_id] = atom_incidence * atoms_area_per_dot[atom_id] 

    if sp_file:
        writing_points_xyz(sp_file, dots_surface)

    centered_dots_surface = dots_surface - np.average(dots_surface, axis=0)
    origin = np.array([[0., 0., 0.]])
    dots_origin_dist = cdist(origin, centered_dots_surface).flatten()
    d=describe(dots_origin_dist)

    print("    Surface atom: \n{}".format(is_surface))
    print("    Atomic exposition: \n{}.".format(exposition))
    print("    Surface points distance to mol geometric centre:")
    print("        Min, Max, Mean, Variance, Skewness, Kurtosis")
    print("        {:2.5f}, {:2.5f}, {:2.5f}, {:2.5f}, {:2.5f}, {:2.5f}".format(d.minmax[0], d.minmax[1], d.mean, d.variance, d.skewness, d.kurtosis))
    print("    Total surface area: {}\n".format(sum(exposition)))
    return is_surface, exposition

def compute_ecn_dav(positions, print_convergence=True):
    """Return the effective coordination number (ecn) and the average bound
    distance and the conective index matrix Pij.
    positions: np.array of (n,3) shape with atoms positions.
    pij_to_int: around the conective index matrix to intiger."""

    if (type(positions) != np.ndarray) or (np.shape(positions)[1] != 3):
        sys.exit("positions must be a (n,3) shaped numpy.ndarray! Aborting...")

    print("ECN analysis:")
    dij = cdist(positions,positions) + 100*np.eye(len(positions))
    dav = np.max(cdist(positions,positions), axis=0)
    ecn = np.zeros_like(dav)
    dav_pre = np.zeros_like(dav)
    ecn_pre = np.zeros_like(dij)
    v=0
    print("    MAE davs       MAE ECNs")
    while np.average(np.abs(dav_pre - dav)) > 10E-8 and np.average(np.abs(dav_pre - ecn)) or v < 2 :
        if v > 0 :
            dav_pre = dav * 1.
            ecn_pre = ecn * 1.
        Pij = np.exp(1-(2* dij / ( dav.reshape(-1,1) + dav.reshape(1,-1)   ))**6)
        ecn = np.sum(Pij, axis=1)
        dav = np.sum(Pij * dij,axis=1) / np.sum(Pij,axis=1)
        ecn = np.sum(Pij, axis=1)
        print('    {:1.10f}   {:1.10f}'.format(np.average(np.abs(dav_pre - dav )), np.average(np.abs(ecn - ecn_pre))))
        v+=1
    print('    ECNs: \n{}'.format(np.array(ecn)))
    print('    davs: \n{}\n'.format(np.array(dav)))
    return ecn , dav , Pij


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='This script calculate the atoms exposition to the vacuum.')
    parser.add_argument('--mol', required=True, 
            help='the molecules (xyz, geometry.in, etc) to analyze.')
    parser.add_argument('--r_adatom', action='store', default=1.1,
            help='the radius for the adatom radius. (Default=1.1)')
    parser.add_argument('--r_atoms', action='store', default='',
            help='the radius for the mol atoms. (Default=dav/2)')
    parser.add_argument('--ssamples', action='store', default=1000,
            help='the number of points to distribute over each atom. (Default=1000)')
    parser.add_argument('--sp_file', action='store', default='',
            help='if defined (ex: surf.xyz), the surface points found will be printed in this file.')
    parser.add_argument('--save_json', action='store', default='',
            help='if defined (ex: dados.json), all the data will be saved in a json file.')
    args = parser.parse_args()
 
    mols_names = args.mol.split()
    adatom_radius = float(args.r_adatom)
    ssamples = int(args.ssamples)
    sp_file = args.sp_file
    json_file = args.save_json

    use_ecn = True
    r_atoms = args.r_atoms.split()
    if len(r_atoms) >= 1:
        use_ecn = False

    # creating variables to save a json file
    if json_file:
        import pandas as pd
        meus_dados = pd.DataFrame()
        all_positions = []
        all_chemical_symbols = []
        all_ecn = []
        all_dav = []
        all_is_surf = []
        all_exposition = []

    for mol_name in mols_names:
        print("Investigation mol: {}\n".format(mol_name))
        mol = ase.io.read(mol_name)
        positions = mol.positions
        chemical_symbols = np.array(mol.get_chemical_symbols())
        size = mol.get_global_number_of_atoms()

        if use_ecn:
            # compute dav
            ecn, dav, pij = compute_ecn_dav(positions, print_convergence=True)
            atomic_radius = dav/2.
        else:
            # processing r_atoms
            if len(r_atoms) == 1:
                atomic_radius = np.array(r_atoms*size, dtype=float)
            else:
                atomic_radius = np.array(r_atoms, dtype=float)
                if len(atomic_radius) != size:
                    sys.exit('ERROR: The number radii provided ({}) deffer from the '\
                             'number of atoms ({}) in the molecule {}.'.format(
                                  len(atomic_radius), size, mol_name))

        # finding surface atoms
        is_surface, exposition = surfcore_classyfier_adatom(positions, atomic_radius, adatom_radius, 
                ssamples=ssamples, sp_file=sp_file)

        # printing:
        print('index chemical_element ecn dav is_surf exposition')
        for i, (chemi_i, ecn_i, dav_i, is_surf_i, expos_i) in enumerate(zip(chemical_symbols, ecn, dav, is_surface, exposition)):
            print(i+1, chemi_i, ecn_i, dav_i, is_surf_i, expos_i)

        # saving mol data
        if json_file:
            all_positions.append(positions)
            all_chemical_symbols.append(chemical_symbols)
            all_ecn.append(ecn)
            all_dav.append(dav)
            all_is_surf.append(is_surface)
            all_exposition.append(exposition)

    # writing the data
    if json_file:
        meus_dados['mols'] = mols_names
        meus_dados['chemical_symbols'] = all_chemical_symbols 
        meus_dados['ecn'] = all_ecn
        meus_dados['dav'] = all_dav
        meus_dados['is_surf'] = all_is_surf
        meus_dados['exposition'] = all_exposition
        print('Writing json_file: {}'.format(json_file))
        meus_dados.to_json(json_file)
