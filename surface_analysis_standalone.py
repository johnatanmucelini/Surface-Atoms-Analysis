import sys
import numpy as np
import ase.io
import logging

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

    logging.debug("Initializing RegRDS_set function!")
    logging.debug("sampling_distance: " + str(sampling_distance))
    logging.debug("N: " + str(N))

    if not sampling_distance > 0:
        logging.error("sampling_distance must be higher than zero! Aborting...")
        sys.exit(1)

    if (type(N) != int) or (N <= 0):
        logging.error("N must be an intiger grater than zero! Aborting...")
        sys.exit(1)

    cart_coordinates=[]
    r = 1
    Ncount = 0
    a = 4. * np.pi * r**2 / N
    d = np.sqrt(a)
    Mtheta = int(round(np.pi/d))
    dtheta = np.pi / Mtheta
    dphi = a / dtheta
    logging.debug("Mtheta: " + str(Mtheta))
    for m in range(0, Mtheta ) :
        theta = np.pi *( m + 0.5 ) / Mtheta
        Mphi = int(round(2 *np.pi * np.sin(theta) / dphi ))
        logging.debug("Mtheta: " + str(Mphi))
        for n in range( 0 , Mphi ) :
            phi = 2* np.pi * n / Mphi
            Ncount += 1
            y = sampling_distance * np.sin(theta) * np.cos(phi)
            x = sampling_distance * np.sin(theta) * np.sin(phi)
            z = sampling_distance * np.cos(theta)
            cart_coordinates.append([x,y,z])
    cart_coordinates = np.array(cart_coordinates)
    logging.info("Final quanity of points in the radial grid: "
                 + str(len(cart_coordinates)))

    logging.debug("RegRDS_set function finished sucsessfuly!")
    return cart_coordinates

def writing_points_xyz( file_name , positions ) :
    """Write positions, a list or array of R3 points, in a xyz file file_named.
    Several softwares open xyz files, such as Avogadro and VESTA
    In the xyz file all the atoms are H.

    file_name: a string with the path of the xyz document which will be writed.
    positions: a list or numpy array with the atoms positions."""

    logging.debug("Initializing writing_points_xyz function!")
    logging.debug("file_name: " + str(file_name))
    logging.debug("positions: " + str(positions))

    if type(file_name) != str :
        logging.error("file_name must be a string")
        sys.exit(1)

    for index, element in enumerate(positions) :
        if len(element) != 3 :
            logging.error("Element " + str(index) + " of positions does not" \
                          "present three elements." )
    if type(positions) != list : positions = np.array(positions)

    logging.debug("Writing points in the xyz file...")
    ase.io.write( file_name , ase.Atoms('H'+str(len(positions)) , list(map(tuple,positions)) ))
    logging.debug("Finished.")

    logging.debug("writing_points_xyz function finished!")

def large_surfaces_index(surface_dots_positions, eps):
    """Seach if there are more than one surfaces in surface_dots_positions, than
    return a np.array of booleans with True for index of atoms for the surfaces
    more points.
    surface_dots_positions: np.array with R3 dots.
    eps: minimun distance between different surfaces, should be a float grater
         than zero."""

    logging.debug("Initializing remove_pseudo_surfaces function!")

    if (not isinstance(surface_dots_positions, np.ndarray)
            or (np.shape(surface_dots_positions)[1] != 3)):
        logging.error("surface_dots_positions must be a (n,3) shaped np.array.\
                  Aborting...")
        sys.exit(0)

    if not eps > 0:
        logging.error("eps must be large than zero.")
        sys.exit(1)

    db = DBSCAN(eps=eps, min_samples=1).fit_predict(surface_dots_positions)
    labels, quanity = np.unique(db, return_counts=True)
    if len(labels) > 1:
        logging.warning(str(len(labels)) + ' surfaces were found, of sizes: '
                        + str(quanity).replace('[', '').replace(']', '')
                        + '. The bigger will be selected!')
    result = db == labels[np.argmax(quanity)]

    logging.debug("remove_pseudo_surfaces finished sucsessfuly!")
    return result

def surfcore_classyfier_adatom(positions, atomic_radii, adatom_radius,
                               remove_is=True, ssamples=1000, writw_sp=True,
                               return_expositions=True,
                               print_surf_properties=False,
                               sp_file="surface_points.xyz"):
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

    # Centralizing atoms positions:
    positions = positions - np.average(positions, axis=0)
    touch_radii = atomic_radii + adatom_radius
    qtna = len(positions)

    dots_try = RegRDS_set(adatom_radius, ssamples)
    rssamples = len(dots_try)
    max_dots = rssamples * qtna
    print('Quantity of dots per atom: ' + str(len(dots_try)))
    print('Quantity of investigated dots: ' + str(max_dots))

    dots = positions[0] + RegRDS_set(touch_radii[0] + 0.001, ssamples)
    dot_origin = [[0] * len(dots_try)]
    for atom_id in range(1, qtna):
        dots = np.append(dots, positions[atom_id] + RegRDS_set(
            touch_radii[atom_id] + 0.001, ssamples), axis=0)
        dot_origin.append([atom_id] * len(dots_try))
    dot_origin = np.array(dot_origin).flatten()

    # removing dots inside other touch sphere
    dots_atoms_distance = cdist(positions, dots)
    atomos_radii_projected = np.array(
        [touch_radii]*len(dots)
        ).reshape(len(dots), qtna).T
    surface_dots = np.sum(
        dots_atoms_distance < atomos_radii_projected,
        axis=0
        ) == 0
    dots_surface = dots[surface_dots]

    # removing internal surfaces
    if remove_is:
        print("removing_internal_surface")
        dots_for_find_eps = RegRDS_set(max(touch_radii), ssamples)
        dotdot_distances = cdist(dots_for_find_eps, dots_for_find_eps)
        min_dotdot_dist = np.min(dotdot_distances + np.eye(rssamples) * 10,
                                 axis=0)
        eps = 2.1 * np.max(min_dotdot_dist)
        external_surface_dots = large_surfaces_index(dots_surface, eps)
        dots_surface = dots_surface[external_surface_dots]
        surface_dot_origin = dot_origin[surface_dots][external_surface_dots]
    else:
        surface_dot_origin = dot_origin[surface_dots]
    surface_atoms, dots_per_atom = np.unique(surface_dot_origin,
                                             return_counts=True)

    # Getting exposition
    atoms_area_per_dot = (4 * np.pi * atomic_radii**2)/(1.*rssamples)
    exposition = np.zeros(qtna)
    is_surface = np.zeros(qtna, dtype=bool)
    incidence = np.zeros(qtna)
    for atom_id, atom_incidence in zip(surface_atoms, dots_per_atom):
        print(' found surface atom: ' + str(atom_id))
        is_surface[atom_id] = True
        incidence[atom_id] = atom_incidence / rssamples
        exposition[atom_id] = atom_incidence * atoms_area_per_dot[atom_id] 

    if writw_sp:
        write_points_xyz(sp_file, dots_surface)

    if print_surf_properties:
        centered_dots_surface = dots_surface - np.average(dots_surface, axis=0)
        origin = np.array([[0., 0., 0.]])
        dots_origin_dist = cdist(origin, centered_dots_surface).flatten()
        print("Surface points description: " + str(describe(dots_origin_dist)))
        print("reg_surface area: " + str(sum(exposition)))

    if not return_expositions:
        return is_surface
    if return_expositions:
        return is_surface, exposition

def compute_ecn_dav ( positions , print_convergence=True ):
    """Return the effective coordination number (ecn) and the average bound
    distance and the conective index matrix Pij.
    positions: np.array of (n,3) shape with atoms positions.
    pij_to_int: around the conective index matrix to intiger."""

    logging.debug("Initializing ECN analysis!")


    if (type(positions) != np.ndarray) or (np.shape(positions)[1] != 3):
        logging.error("positions must be a (n,3) shaped numpy.ndarray!"
                      + " Aborting...")
        sys.exit(1)

    logging.info("\n\nECN analysis:")
    dij = cdist(positions,positions) + 100*np.eye(len(positions))
    dav = np.max( cdist(positions,positions) , axis=0)
    dav_pre = np.zeros_like(dav)
    ecn_pre = np.zeros_like(dij)
    v=0
    logging.info("    \Delta sum_i(abs(dav_i))/N    \Delta sum_i(abs(ECN_i))/N")
    while np.sum(np.abs(dav_pre - dav ))/ len(dav) > 10E-8 or v < 2 :
        if v > 0 :
            dav_pre = dav * 1.
            ecn_pre = ecn * 1.
        Pij = np.exp(1-(2* dij / ( dav.reshape(-1,1) + dav.reshape(1,-1)   ))**6)
        ecn = np.sum(Pij, axis=1)
        dav = np.sum(Pij * dij,axis=1) / np.sum(Pij,axis=1)
        ecn = np.sum(Pij, axis=1)
        logging.info('   ' + str(np.sum(np.abs(dav_pre - dav ))/ len(dav))
                     + '  ' + str(np.sum(np.abs(ecn - ecn_pre))/len(dav)) )
        v+=1
    logging.info( '    bag_ecn: ' + '['+','.join(np.array( ecn , dtype=str))+']' )
    logging.info( '    bag_dav: ' + '['+','.join(np.array( dav , dtype=str))+']' )
    logging.info( '    ECN analysis finished...' )
    return ecn , dav , Pij



# reading info
cluster = ase.io.read(sys.argv[1])
if len(sys.argv) == 3:
   adatom_radius = float(sys.argv[2])
else: 
    adatom_radius = 1.1
positions = cluster.positions
chemical_symbols = np.array(cluster.get_chemical_symbols())

# compute dav
ecn, dav, pij = compute_ecn_dav(positions, print_convergence=False)

# finding surface atoms
surf_atoms = surfcore_classyfier_adatom(positions, dav/2., adatom_radius, ssamples=1000,
                                        writw_sp=False,
                                        return_expositions=False,
                                        print_surf_properties=True)

