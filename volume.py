import snappy
import numpy as np
import symmetries
from cypari2 import PariError

def get_volume(M, x, y):
    """
    Compute the volume of the (x, y)-Dehn filling of a given snappy manifold M.
    - The threshold value 0.942 is the volume of the smallest closed hyperbolic 3-manifold.
    - The computed volume has approximately 10 significant digits by default.
    """
    N = snappy.Manifold(M)
    N.dehn_fill((x, y))
    try:
        volume = float(N.volume())

        # Exclude non-hyperbolic fillings
        if volume < 0.942:
            return 0
        return volume

    except (ValueError, PariError) as e:
        print(f" Volume error for [{x}, {y}] ({e})", end='')
        return 0

def get_prec_volume(M, x, y):
    """
    Compute the volume of the (x, y)-Dehn filling of a given snappy manifold M with higher precision (64-bit float).
    """
    N = snappy.Manifold(M)
    N.dehn_fill((x,y))
    try:
        return float(N.volume(bits_prec=64))
    except (ValueError, PariError) as e:
        print(f" Prec_volume error for [{x}, {y}] ({e})", end='')
        return 0

def get_high_prec_volume(M, x, y):
    """
    Compute the volume of the (x, y)-Dehn filling of a given snappy manifold M in the snappy's high-precision setting.
    - Uses 212-bit precision (~63 decimal digits).
    """
    N = snappy.ManifoldHP(M)
    N.dehn_fill((x, y))
    try:
        return N.volume()
    except (ValueError, PariError) as e:
        print(f" High_prec_volume error for [{x}, {y}] ({e})", end='')
        return 0

def add_unique(lst, item):
    if item not in lst:
        lst.append(item)
    
def test_with_symmetry(M, symm, numb, tolerance, volume_matrix, prec_volume_matrix):
    """
    For a given snappy manifold M, this function identifies pairs of Dehn fillings of M
    that have the same volume but are not resulted from the known symmetries of M's volume formula.

    Inputs:
        M : A snappy manifold.
        symm : A list of known symmetries of M's volume formula.
        numb : The range of the filling coefficients.
        tolerance : The acceptable error in the computed volume.
        volume_matrix : A matrix storing the volumes of each filling to avoid redundant computations.
        prec_volume_matrix : A matrix storing the prec-volumes of each filling.
    """

    # Temporary variable for counting volume matchings.
    temp = 0
    found = False
    
    # This list stores filling coefficients that have already been checked by previous coefficients.
    multiple_indices = []
    
    # The space of Dehn filling coefficients exhibits central symmetry.
    rangex = range(-numb, numb + 1)
    rangey = range(0, numb + 1)
    
    for row_idx, x in enumerate(rangex):
        if found:
            break
        for col_idx, y in enumerate(rangey):
            if volume_matrix[row_idx, col_idx] > 0 and [x,y] not in multiple_indices:
                diff_matrix = volume_matrix - volume_matrix[row_idx, col_idx]
                cond_matrix = np.abs(diff_matrix) <= tolerance
                # This matrix stores 1 for entries whose volume is close enough to M(x, y)'s volume, 
                # and 0 otherwise.

                # This list stores filling coefficients that have the same volume as M(x, y) 
                # due to the symmetry of the volume formula.
                symm_indices = []
                for j in range(0,len(symm)):
                    n = symm[j][4]
                    x0 = symm[j][0] * x + symm[j][1]*y
                    y0 = symm[j][2] * x + symm[j][3]*y
                    if x0 % n == 0 and y0 % n == 0:
                        x1 = int(x0/n)
                        y1 = int(y0/n)
                    # For the structure of symmetries, see 'symmetries.py'.
                        
                        if abs(x1) <= numb and abs(y1) <= numb:
                            if 0 < y1 and [x,y] != [x1,y1]:
                                add_unique(symm_indices,[x1,y1])
                                cond_matrix[x1+numb, y1] = 0
                            elif y1 < 0 and [x,y] != [-x1,-y1]:
                                add_unique(symm_indices,[-x1,-y1])
                                cond_matrix[-x1+numb, -y1] = 0
                            elif y1 == 0 and [x,y] != [abs(x1),0]:
                                add_unique(symm_indices,[abs(x1),0])
                                cond_matrix[abs(x1)+numb, 0] = 0
                        # This excludes the matching of M(x, y) and M(ax + by, cx + dy), 
                        # where (a, b, c, d) is a known symmetry.
                        
                multiple_indices = multiple_indices + symm_indices
                
                matching_indices = np.argwhere(cond_matrix)
                if len(matching_indices) > 1:
                    matching_indices_list = matching_indices.tolist()
                    
                    # If the prec_volume is already computed, retrieve it from prec_volume_matrix.
                    # Otherwise, compute it.
                    if prec_volume_matrix[row_idx, col_idx] == 0:
                        volume = get_prec_volume(M, x, y)
                        prec_volume_matrix[row_idx, col_idx] = volume
                    else:
                        volume = prec_volume_matrix[row_idx, col_idx]
                        
                    HPvolume = 0
                        
                    # While iterating through the matching list, we remove coefficients if 
                    # their prec or high-precision volume differs from that of M(x, y).
                    # Due to a technical issue, we iterate through the list in reverse order.
                    for i in range(len(matching_indices_list) - 1, -1, -1):
                        matching_indices_list[i][0] -= numb
                        x1 = matching_indices_list[i][0]
                        y1 = matching_indices_list[i][1]
                        if [x,y] != [x1,y1]:
                            if prec_volume_matrix[x1+numb, y1] == 0:
                                volume1 = get_prec_volume(M, x1, y1)
                                prec_volume_matrix[x1+numb, y1] = volume1
                            else:
                                volume1 = prec_volume_matrix[x1+numb, y1]
                            if abs(volume - volume1) > 1e-15:
                                matching_indices_list.pop(i)
                            
                            # If the prec_volumes match, then compare the high-precision volumes.
                            else:
                                if HPvolume == 0:
                                    HPvolume = get_high_prec_volume(M, x, y)
                                HPvolume1 = get_high_prec_volume(M, x1, y1)
                                if abs(HPvolume - HPvolume1) > 1e-62:
                                    matching_indices_list.pop(i)
                                else:
                                    multiple_indices.append([x1,y1])
      
                    if len(matching_indices_list) > 1:
                        # Print and count the exceptional matchings.
                        print(' ' + str(matching_indices_list + symm_indices), end='')
                        temp += 1
                        
                        # If the number of the exceptional matchings exceeds five, exit the search.
                        if temp >= 5:
                            found = True
                            break
    
    return found

def test(manifolds, numb):
    """
    The main function for the experiments.
    
    Inputs:
        manifolds : A list of snappy manifolds
        numb : The range of the filling coefficients.
    """
    rangex = range(-numb, numb + 1)
    rangey = range(0, numb + 1)
    tolerance = 1e-9
    
    for M in manifolds:
        M = snappy.Manifold(M)
        print(M.name(), end='')
        
        volume_matrix = np.zeros((len(rangex), len(rangey)), dtype=np.double)
        prec_volume_matrix = np.zeros((len(rangex), len(rangey)), dtype=np.double)
        
        # Compute the volume for each Dehn filling.
        for row_idx, x in enumerate(rangex):
            for col_idx, y in enumerate(rangey):
            
                # M(x,0) and M(-x,0) represent the same filling.
                if y==0 and x <= 0:
                    volume = 0
                else:
                    volume = get_volume(M, x, y)
                volume_matrix[row_idx, col_idx] = volume
        
        sym = []
        # If the manifold has identified symmtries, retrieve them from 'symmetries'.
        for symm in symmetries.symm_list:
            if symm[0] == M.name():
                print(" Identified symmetries:", end='')
                for i in range(1,len(symm)):
                    print(' ' +str(symm[i]), end='')
                sym = symm[1:]

        found = test_with_symmetry(M, sym, numb, tolerance, volume_matrix, prec_volume_matrix)
        if found:
            print(" Unidentified symmetries?")
        else:
            print("")
