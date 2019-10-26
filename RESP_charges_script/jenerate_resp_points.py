# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import matplotlib
matplotlib.use('Agg')
import sys
import itertools
import math
import numpy as np

#import prettyplotlib as ppl
import matplotlib.pyplot as plt

# <codecell>

# python script QMlogfile density 

def gaussian_reader(file_handle):
    """read file_handle as a gaussian output file
    retunrs a list of coordinates
    """
    confs = []
    for line in file_handle:
        if line.strip() == "Input orientation:":
            [file_handle.readline() for skip_lines in range(4)]
            confs.append([])
            for coord_line in file_handle:
                coord_line = coord_line.strip()
                if coord_line == "---------------------------------------------------------------------":
                    break
                coord_line = coord_line.split()
                confs[-1].append(dict(atomic_no=int(coord_line[1]), xyz=tuple(float(coord_line[i]) for i in range(3, 6))))
    return confs

def gamess_reader(file_handle):
    """reads file_handle as a gamess output file
    returns a list of coordinates
    """
    confs = []
    for line in file_handle:
        if line.strip() == "***** EQUILIBRIUM GEOMETRY LOCATED *****":
            while file_handle.readline().strip() != "------------------------------------------------------------": pass
            confs.append([])
            for coord_line in file_handle:
                coord_line = coord_line.strip()
                if not coord_line:
                    break
                coord_line = coord_line.split()
                confs[-1].append(dict(atomic_no=int(float(coord_line[1])), xyz=tuple(float(coord_line[i]) for i in range(2, 5))))

def el_xyz_reader(file_handle):
    """reads file_handle as a file of format <element> <x> <y> <z>
    returns a list of 1 set of coordinates
    """
    confs = [[]]
    for line in file_handle:
        line = line.strip().split()
        confs[-1].append(dict(element=line[0], xyz=tuple(float(line[i]) for i in range(1,4))))
    return confs

def orca_reader(file_handle):
    """reads file_handle as an orca output file
    returns a list of coordinates
    """
    confs = []
    for line in file_handle:
        if line.strip() == "CARTESIAN COORDINATES (ANGSTROEM)":
            file_handle.readline()
            confs.append([])
            for coord_line in file_handle:
                coord_line = coord_line.strip()
                if not coord_line:
                    break
                coord_line = coord_line.split()
                confs[-1].append(dict(element=coord_line[0], xyz=tuple(float(coord_line[i]) for i in range(1, 4))))
    return confs


def dist_sq(a, b):
    """returns the squared distance between points with coordinates a and b
    
    >>> dist_sq((1, 0, 0), (0, 0,0))
    1

    >>> dist_sq((1, 0, 0), (0, 1, 0))
    2

    >>> dist_sq([-1.76395651, -1.80511352,  0.39569726], [-0.96790506,  0.75137003, -1.72923407])
    11.684639209683272
    """
    return sum((x - y) * (x - y) for x, y in zip(a, b))

class Atom:
    """Represents an atom, right now only prints points at a certain radius"""
    # first 4 periods of periodic table
    elements = ["N/A", "H", "He",
                "Li", "Be", "B", "C", "N", "O", "F", "Ne",
                "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
                "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr"
                ]
    # some radii from Bondi 1962 of common/necessary organic atoms
    element_radii = {'H': 1.2, 'C': 1.7, 'N': 1.55, 'O': 1.52, 'S': 1.8, 'P':1.95, 'Se': 1.9}

    def __init__(self, xyz, atomic_no=None, element=None):
        """Initializes Atom with xyz coordinates, using atomic_no for atomic
           number or, failing that, element as a string for the element name
        """
        self.xyz = xyz
        try:
            self.element = Atom.elements[atomic_no]
            self.atomic_no = atomic_no
        except:
            self.atomic_no = Atom.elements.index(element)
            self.element = element
        self.radius = Atom.element_radii[self.element]
    
    def radius_points(self, scale=1, density=1, check_against=[]):
        """return points at the distance of radius multiplied by scale, with a
           density on the surface of density (points / square angstrom)
           if check_against is specified, it won't generate points within the
           scaled radius of the atoms in check_against
        """
        r = self.radius * scale
        # get other atoms to check against
        to_check = []
        for other in check_against:
            sum_r_sq = other.radius * scale + r
            sum_r_sq *= sum_r_sq

            if dist_sq(self.xyz, other.xyz) < sum_r_sq:
                to_check.append((other, (other.radius * scale)**2))
        # get number of points based on density and surface area
        # but this will be for each eight of the sphere
        num_points = math.pi * r * r * density / 2
        # inclination radius is always sphere's r
        # azimuth radius varies from 0 to sphere r
        # in fact, it's r * sin(theta), where theta is the inclination angle
        # the average of this is 2 / PI * r
        # so the average radii are r and 2 / PI * r
        # so the circumferences are 2 PI r and 2 PI (2 / PI r)
        # but the half-circumference of incl is PI r and the azimuthal circumference is still 2 PI (2 / PI r)
        # hence the surface area is 4 PI r^2, which is the textbook definition
        # but if we want to get number of points in each direction, we need to consider
        # x = PI r, and SA = x * (4 / PI) x
        # so x = sqrt(SA * PI / 4)
        # if instead of SA, we want to divide number of points, we can just substitute that, as below
        
        # partition number of points by inclination and azimuth d.o.f.
        inc_points = math.sqrt(num_points * math.pi / 4)

        # round up for both
        inc_points, azi_points_0 = math.ceil(inc_points), inc_points * 2 / math.pi
        avg_r = 2 / math.pi * r

        points = []
        counter = 0
        # pre-calculate some things
        theta_factor = math.pi / 2 / inc_points
        z_factorses = [(1, -1)] * inc_points + [(1, )]
        xy_factorses = [(1, )] + [(1, -1)] * inc_points
        # iterate over inclination points
        # + 1 to make sure we get top and bottom
        for inc, z_factors, xy_factors in zip(range(inc_points + 1), z_factorses, xy_factorses):
            theta = inc * theta_factor
            z_0 = r * math.cos(theta)
            inc_r = r * math.sin(theta)
            # get number of azimuthal points at this inclination
            azi_points = math.ceil(azi_points_0 * (inc_r + 0.001) / avg_r)
            for z in (z_0 * z_factor for z_factor in z_factors):
                counter += azi_points
                #print(azi_points)
                phi_factor = math.pi / 2 / azi_points

                for azi in range(azi_points):
                    phi = azi * phi_factor
                    x = inc_r * math.cos(phi)
                    y = inc_r * math.sin(phi)
                    for xy_rot in range(4 if x or y else 1):
                        final_offset = (x, y, z)
                        loc = tuple(a + b for a, b in zip(self.xyz, final_offset))
                        #print(math.sin(phi), inc_r)
                        #tuple(r_sin_theta * i for i in azi) + (z, )
                        for other, other_r_sq in to_check:
                            if other != self and dist_sq(loc, other.xyz) < other_r_sq:
                                break
                        else:
                            points.append(loc)
                        x, y = -y, x
        assert counter > num_points
        return points

# open and read file
with open(sys.argv[1]) as qm_file:
    line1 = qm_file.readline()
    xyzscale = 1
    if "GAMESS" in line1:
        print("GAMESS")
        confs = gamess_reader(qm_file)
    else:
        qm_file.readline()
        line1 = qm_file.readline()
        if "* O   R   C   A *" in line1:
            print("Orca")
            confs = orca_reader(qm_file)
            #convert to Bohr
            xyzscale = 1.889725989
        else:
            print("Gaussian")
            confs = gaussian_reader(qm_file)
    
    atoms = [Atom(**a) for a in confs[-1]]
    x, y, z = np.array([c['xyz'] for c in confs[-1]]).T

    # array of points for each scale of molecular surface
    #try different density, 6 is a good number
    pointses = []
    colors = []
    for scale in range(14, 21, 2):
        colors.append(scale)
        pointses.append([])
        for atom in atoms:
            pointses[-1].extend(atom.radius_points(density=int(sys.argv[2]), scale=scale/10, check_against=atoms))
    
    grand_colors = list(itertools.chain.from_iterable([color] * len(points) for color, points in zip(colors, pointses)))
    
    # make numpy arrays from pointses
    pointses = [np.array(points) for points in pointses]
    
    # go through all pointses, saving points to a file and plotting
    for color, points in zip(colors, pointses):
        np.savetxt('esppoints_{}'.format(color), points * xyzscale, header=str(len(points)), comments='')
        plt.figure(figsize=(10, 10))
        plt.scatter(points[:, 0], points[:, 1], c=points[:, 2])
        plt.axis('equal')
        plt.savefig('esppoints_{}.png'.format(color))
    
    #collect all points together for a giant scatter plot
    grand_points = np.vstack(pointses)
    
    
    for axes in (0, 1, 2), (1, 2, 0), (0, 2, 1):
        plt.figure(figsize=(10, 10))
        plt.scatter(grand_points[:, axes[0]], grand_points[:, axes[1]], c=grand_colors)
        plt.xlabel('xyz'[axes[0]])
        plt.ylabel('xyz'[axes[1]])
        plt.axis('equal')
        plt.savefig('esppoints_all_{}_{}_{}.png'.format(*axes))


