#####################################################################################
#####################################################################################

# This file is a companion to the paper "A statistical approach to knot confinement 
# via persistent homology" by D.Celoria and B.I. Mahler.
# What follows is a list of functions created to compute several quantities associated
# to PL embeddings of knots generated through Topoly (https://pypi.org/project/topoly/).
# The complete dataset generated used, together with the persistent homology groups computed
# with Ripser (https://github.com/Ripser/ripser) is available upon request to the authors.

#####################################################################################
#####################################################################################

# required imports
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import re
from ripser import ripser
import numpy as np
from ripser import Rips


# Takes as "points_list" input a Pl embedding of a knot (with unit length edges);
# The first part of the output consists of the integral of the first Betti curve 
# of the VR  persistent homology of the point cloud obtained by interpolating the
# endpoints of the PL knot with 10 points. 
# The second part of the output is the maximal value attained by the curve, or in 
# other words, the maximal number of bars that are alive at the same time.
# If the "draw_graph" is True, it also displays a picture of the Betti curve.

def betti_curve_features(points_list,draw_graph = True):
    all_points = [float(aa[i]) for aa in points_list for i in [0,1]]
    all_points = list(set(all_points))
    all_points.sort()
    vector_values = []
    epsilon = float(10**(-8))
    for tt in all_points[:-1]:
        contatore = 0
        for segmento in points_list:
            if segmento[0] < tt + epsilon < segmento[1]:
                contatore +=1
        vector_values.append(contatore)
    if draw_graph == True:
        x = [0]
        for yy in all_points:
            x.append(yy)
            x.append(yy)
        y = [0,0]
        for uu in vector_values:
            y.append(uu)
            y.append(uu)
        y.append(0)
        plt.plot(x,y)
        plt.title('betti curve')
        plt.show()
    integral = 0
    for rrr in range(len(all_points)-1):
        integral += (all_points[rrr+1] - all_points[rrr])*vector_values[rrr]
    return [integral, max(vector_values)]


# Takes as "points_list" input a Pl embedding of a knot (with unit length edges);
# outputs the curvature of the knot.

def compute_curvature(points_list):
    N = len(points_list)
    aux_angles = []
    for i in range(N):
        first = difference_arrays(points_list[(i+1)%N],points_list[i%N] )
        first = [el/np.linalg.norm(first) for el in first]
        second = difference_arrays(points_list[(i+2)%N],points_list[(i+1)%N] )
        second = [el/np.linalg.norm(second) for el in second]
        aux_angles.append(np.arccos(np.dot(first, second)))
    return float(sum(aux_angles)/(2*np.pi))


# Takes as "points_list" input a Pl embedding of a knot (with unit length edges);
# outputs the torsion of the knot.

def compute_torsion(points_list):
    N = len(points_list)
    aux_angles = []
    for i in range(N):
        first = difference_arrays(points_list[(i+1)%N],points_list[i%N] )
        first = [el/np.linalg.norm(first) for el in first]
        second = difference_arrays(points_list[(i+2)%N],points_list[(i+1)%N] )
        second = [el/np.linalg.norm(second) for el in second]
        third = difference_arrays(points_list[(i+3)%N],points_list[(i+2)%N] )    
        third = [el/np.linalg.norm(third) for el in third]
        first_normal = np.cross(first,second)
        first_normal = [el/np.linalg.norm(first_normal) for el in first_normal]
        second_normal = np.cross(second, third)
        second_normal = [el/np.linalg.norm(second_normal) for el in second_normal]
        aux_angles.append(np.arccos(np.dot(first_normal,second_normal)))
    return float(sum(aux_angles)/(2*np.pi))

# Takes as "points_list" input a Pl embedding of a knot (with unit length edges);
# outputs the maximal distance between two points in the knot (in other words the diameter
# of the circumscribing sphere).

def compute_max_distance(points_list):
    massimo = 0
    for ii in range(len(points_list)):
        for jj in range(len(points_list)):
            if ii > jj:
                dist = np.linalg.norm(difference_arrays(points_list[ii],points_list[jj]))
                if dist > massimo:
                    massimo = dist
    return massimo 
    
# Takes as "points_list" input a Pl embedding of a knot (with unit length edges);
# outputs the radius of gyration of the knot.

def radius_of_gyration(points_list):
    aux = 0
    for pp in points_list:
        for qq in points_list:
            aux += np.linalg.norm(difference_arrays(pp,qq))**2
    return aux/(2*(len(points_list)**2))
    
    
# This function takes as input two points in R^3 with distance 1, and interpolates them
# with "number_of_points" further points. This is used to create the point cloud from a PL
# knot embedding, by applying this function to each unit length edge.

def create_necklace_around_segment(startpoint, endpoint, number_of_points):
    new_points = []
    segment_difference = np.array([endpoint[i] - startpoint[i] for i in range(3)])
    segment_difference /= np.linalg.norm(segment_difference)
    first_orth = np.array([1,1,1], dtype = float)  #create a random vector
    first_orth -= first_orth.dot(segment_difference)*segment_difference
    first_orth /= np.linalg.norm(first_orth)
    second_orth = np.cross(segment_difference, first_orth)     #apply gram schmidt
    second_orth /= np.linalg.norm(second_orth)
    M = np.matrix([segment_difference,first_orth,second_orth]).transpose()  #matrix to rotate the points 
    for i in range(number_of_points):
        aux = np.array([(i+1)/number_of_points,0, 0])        
        aux2 =  M*np.array([aux[0],aux[2]*math.cos(aux[1]),aux[2]*math.sin(aux[1])]).reshape((3,1))
        new_points.append([aux2[j] + startpoint[j] for j in range(3)])
    return(new_points)

# This is just an auxiliary function that takes the difference of two 3d arrays
def difference_arrays(one, two):
    return [one[t] - two[t] for t in range(3)]

    
