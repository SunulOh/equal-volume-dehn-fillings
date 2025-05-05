import ast
import snappy
import volume

"""
This program converts the text file "symmetries.txt" into a Python list.
Each element of the list represents a manifold's name followed by its volume formula symmetries.

Each symmetry is a list of five integers [a, b, c, d, n], which corresponds to the matrix:
    [a/n  b/n]
    [c/n  d/n]  
in PSL(2, Q), where a*d - b*c = n².

Example:
    ['m136', [0, 4, 1, 0, 2], [1, 0, 0, -1, 1], [0, -4, 1, 0, 2]]
"""

# Load symmetry data from "symmetries.txt"
with open("symmetries.txt", "r") as file:
    symm_list = [ast.literal_eval(line.strip()) for line in file]
    
def check(p, q):
    """
    Checks whether the volume symmetries hold for the given Dehn filling n(p, q).

    The function verifies if applying known symmetries preserves the computed volume.
    If a discrepancy is found, it reports an issue with the symmetry data.
    """
    found_error = False

    for symm in symm_list:
        manifold_name = symm[0]
        M = snappy.Manifold(manifold_name)

        for i in range(1, len(symm)):
            a, b, c, d, n = symm[i]

            # Compute the volumes of M(np, nq) and its symmetric counterpart
            volume_original = volume.get_high_prec_volume(M, p * n, q * n)
            volume_symmetric = volume.get_high_prec_volume(M, p * a + q * b, p * c + q * d)

            # If the volume does not match, report an error
            if volume_original != volume_symmetric:
                found_error = True
                print(f"Symmetry error detected for {manifold_name} at n({p}, {q}).")

    if not found_error:
        print(f"All symmetries verified for n({p}, {q}).")
