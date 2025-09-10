import argparse
import os
from Bio.Align import substitution_matrices
from Bio import SeqIO
import alignement
import score 
import tree

"""
Alignement multiple de séquences par alignement progressif.
Code créer à partir de l'algorithme décrit dans l'article Desmond, Higgings, Sharp (1989).

Prends en entrée un dossier contenant des séquences au format FASTA (protéique ou ADN),
Un dossier de sortie outdir et le mode de lancement ('DNA' ou 'protein').

Crée les fichiers suivants dans outdir :
    - msa.fasta : le fichier FASTA de l'alignement multiple
    - arbre.png : représentation graphique de l'arbre guide
"""

###############################
# Récupération des arguments  #
###############################


def parse_args():
    parser = argparse.ArgumentParser(description="Alignement multiple de séquences (format FASTA).")
    parser.add_argument('--input', type=str, required=True, help='Dossier d\'entrée contenant les fichiers FASTA')
    parser.add_argument('--output', type=str, required=True, help='Dossier contenant les fichiers de sorties')
    parser.add_argument('--mode', type=str, choices=['DNA', 'protein'], default='DNA', help='Type de séquences (DNA ou protein)')
    return parser.parse_args()

###############################
#          Main               #
###############################

if __name__ == "__main__":
    args = parse_args()
    input = args.input
    output = args.output

    #########################
    # Lecture des séquences #
    #########################
    Seq = []

    for filename in os.listdir(input):
        if filename.endswith(".fasta") or filename.endswith(".fa") or filename.endswith(".fna"):
            filepath = os.path.join(input, filename)
            for record in SeqIO.parse(filepath, "fasta"):
                Seq.append(record)

    if not Seq:
        print("Erreur : Aucun fichier FASTA valide trouvé dans le dossier d'entrée.")
        exit(1)

    #########################
    # Choix des paramètres  #
    #########################

    if args.mode == 'protein':
        blosum62 = substitution_matrices.load("BLOSUM62")
        score_fn = lambda a, b: score.protein_score(a, b, blosum62)
        gap_open = 10
        gap_extend = 1
        alphabet = "ACDEFGHIKLMNPQRSTVWY-"
        protein = True
    else:
        score_fn = score.score_DNA
        gap_open = 15
        gap_extend = 6
        alphabet = "ACGT-"
        protein = False

    ##########################
    # Alignement multiple    #
    ##########################
    print("Construction de l'arbre guide")
    guide_tree, Z = tree.build_tree(Seq, gap_open, gap_extend, score_fn=score_fn, alphabet=alphabet, protein=protein)  # <-- récupère racine et matrice linkage
    print("Arbre créer, début de l'alignement progressif")
    msa, ids = alignement.progressive_align_from_tree(guide_tree, Seq, score_fn, gap_open=gap_open, gap_extend=gap_extend, alphabet=alphabet)

    ##########################
    # écriture des résultats #
    ##########################


    with open(os.path.join(output ,"msa.fasta"), "w") as f:
        for sid, aligned_seq in zip([record.id for record in Seq], msa):
            f.write(f">{sid}\n{aligned_seq}\n")

    labels = [str(seq.id) for seq in Seq]

    newick_str = tree.linkage_to_newick(Z, labels)

    with open(os.path.join(output, "clustal.dnd"), "w") as f:
        f.write(newick_str + "\n")

    tree.export_tree_png_biopython(newick_str, os.path.join(output, "arbre.png"))
    alignement.write_clustal(msa, ids, os.path.join(output, "alignment.aln"), protein=protein)