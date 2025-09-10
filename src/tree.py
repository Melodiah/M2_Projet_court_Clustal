from scipy.cluster.hierarchy import linkage, to_tree
import alignement
from Bio import Phylo
import matplotlib.pyplot as plt
from io import StringIO
import numpy as np

"""
Module servant à la création et la sauvegarde d'arbre guide par la méthode UPGMA
dans le cadre d'un alignement multiples de séquences.
"""
################################
# Création d'un arbre guide   #
################################

def build_tree(sequences, gap_open, gap_extend, score_fn, alphabet, protein=False):
    """
    Construit un arbre guide à partir des alignements pairwise.
    
    sequences : str -> séquences à aligner
    score_fn : fonction de score caractère-caractère
    gap_open, gap_extend : pénalités de gap
    protein ( boolean ) : Si True, utilise une table de conservation protéique.

    return : 
        Un arbre guide créer par méthode UPGMA
    """
    n = len(sequences)
    matrix = np.zeros((n, n))

    # Déterminer un score max/min attendu pour la normalisation
    if protein:
        # approx : score max en BLOSUM62 ≈ 11
        max_sub = 11
    else:
        max_sub = 1  # ADN : match = 1
    min_penalty = gap_open  # borne basse si tout gap

    for i in range(n):
        for j in range(i + 1, n):
            score = alignement.needleman_wunsch_affine(
                sequences[i], sequences[j],
                score_fn=score_fn,
                gap_open=gap_open,
                gap_extend=gap_extend
            )

            # Normalisation robuste
            L = min(len(sequences[i]), len(sequences[j]))
            max_score = L * max_sub
            min_score = -max(len(sequences[i]), len(sequences[j])) * min_penalty

            # Échelle en [0,1]
            identity = (score - min_score) / (max_score - min_score)
            identity = max(0, min(1, identity))  # borne finale

            dist = 1 - identity
            matrix[i, j] = dist
            matrix[j, i] = dist

    # Conversion pour scipy
    condensed = matrix[np.triu_indices(n, k=1)]
    Z = linkage(condensed, method="average")
    root_node, _ = to_tree(Z, rd=True)
    return root_node, Z

def linkage_to_newick(Z, labels):
    """
    Fonction formant un dendrogrammes à partir d'un arbre guide.
    """
    root_node, nodes = to_tree(Z, rd=True)

    def build_newick(node, newick, parentdist, leaf_names):
        if node.is_leaf():
            return "%s:%.5f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
        else:
            if newick != "":
                newick = "):%.5f%s" % (parentdist - node.dist, newick)
            else:
                newick = ");"
            newick = build_newick(node.get_left(), newick, node.dist, leaf_names)
            newick = build_newick(node.get_right(), ",%s" % newick, node.dist, leaf_names)
            newick = "(%s" % newick
            return newick

    return build_newick(root_node, "", root_node.dist, labels)

def export_tree_png_biopython(newick_str, filename):
    """
    Création d'un fichier png contenant l'image mise en paramètres
    """
    handle = StringIO(newick_str)
    tree = Phylo.read(handle, "newick")
    fig, ax = plt.subplots(figsize=(8, 6))
    Phylo.draw(tree, axes=ax, do_show=False)
    fig.savefig(filename, dpi=300)

    plt.close(fig)
