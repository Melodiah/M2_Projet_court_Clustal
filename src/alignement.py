import numpy as np
import profil

"""
Module utilisé pour l'alignement de séquences, contenant différentes fonctions d'alignement global ou multiple ainsi qu'une écriture de fichier de sortie.
"""
##################################
# Fonction d'alignement pairwise #
##################################

def needleman_wunsch_affine(seq1, seq2, score_fn, gap_open=15, gap_extend=6):
    """
    Alignement global entre deux séquences par Needleman-Wunsch avec pénalités affines, afin de récupèrer le score optimal.
    
    seq1, seq2 : str -> séquences à aligner
    score_fn : fonction de score caractère-caractère
    gap_open, gap_extend : pénalités de gap

    return : score de l'alignement optimal
    """

    n, m = len(seq1), len(seq2)
    INF = float("inf")  # valeur très grande pour initialiser
    
    # --- 1. Initialisation des matrices ---
    M = np.full((n+1, m+1), INF)  # cas général (match/mismatch)
    X = np.full((n+1, m+1), INF)  # gap dans seq2
    Y = np.full((n+1, m+1), INF)  # gap dans seq1

    # condition de départ : alignement vide
    M[0,0] = 0
    
    # première colonne (alignement avec que des gaps dans seq2)
    for i in range(1, n+1):
        X[i,0] = gap_open + (i-1) * gap_extend
    
    # première ligne (alignement avec que des gaps dans seq1)
    for j in range(1, m+1):
        Y[0,j] = gap_open + (j-1) * gap_extend

    # --- 2. Remplissage des matrices ---
    for i in range(1, n+1):
        for j in range(1, m+1):
            # coût de matcher seq1[i-1] avec seq2[j-1]
            s = score_fn(seq1[i-1], seq2[j-1])
            
            # M = meilleur des cas où on vient de M, X ou Y
            M[i,j] = min(
                M[i-1,j-1] + s,  # prolongement d'un match
                X[i-1,j-1] + s,  # on ferme un gap dans seq2
                Y[i-1,j-1] + s   # on ferme un gap dans seq1
            )
            
            # X = on ouvre ou on prolonge un gap dans seq2
            X[i,j] = min(
                M[i-1,j] + gap_open,   # ouverture d'un gap
                X[i-1,j] + gap_extend  # extension d'un gap existant
            )
            
            # Y = on ouvre ou on prolonge un gap dans seq1
            Y[i,j] = min(
                M[i,j-1] + gap_open,   # ouverture d'un gap
                Y[i,j-1] + gap_extend  # extension d'un gap existant
            )

    # --- 3. Score final ---
    score = min(M[n,m], X[n,m], Y[n,m])
    
    return score

def needleman_wunsch_profiles(prof1, prof2, score_fn, gap_open=15, gap_extend=6):
    """
    Alignement global entre deux profils par Needleman-Wunsch avec pénalités affines.
    
    prof1, prof2 : list[dict] -> profils (colonnes de fréquences)
    score_fn : fonction de score caractère-caractère
    gap_open, gap_extend : pénalités de gap
    
    return : 
        alignement prof1, alignement prof2 (list de colonnes alignées)
    """
    n, m = len(prof1), len(prof2)
    INF = float("inf")

    # Trois matrices : M = match, X = gap dans prof2, Y = gap dans prof1
    M = np.full((n+1, m+1), -INF)
    X = np.full((n+1, m+1), -INF)
    Y = np.full((n+1, m+1), -INF)
    traceback = np.zeros((n+1, m+1), dtype=int)  # 0=M, 1=X, 2=Y

    M[0,0] = 0
    for i in range(1, n+1):
        X[i,0] = -(gap_open + (i-1)*gap_extend)
    for j in range(1, m+1):
        Y[0,j] = -(gap_open + (j-1)*gap_extend)

    # Remplissage
    print("Alignement sur ",n, " bases")
    for i in range(1, n+1):
        for j in range(1, m+1):
            s = profil.score_columns(prof1[i-1], prof2[j-1], score_fn)
            choices = [M[i-1,j-1]+s, X[i-1,j-1]+s, Y[i-1,j-1]+s]
            M[i,j] = max(choices)
            traceback[i,j] = np.argmax(choices)
            X[i,j] = max(M[i-1,j]-(gap_open+gap_extend), X[i-1,j]-gap_extend)
            Y[i,j] = max(M[i,j-1]-(gap_open+gap_extend), Y[i,j-1]-gap_extend)
    print("Fin de l'alignement")
    # écriture de l'alignement
    aln1, aln2 = [], []
    i, j = n, m
    state = np.argmax([M[i,j], X[i,j], Y[i,j]])
    while i>0 or j>0:
        if state == 0:  # M
            s = profil.score_columns(prof1[i-1], prof2[j-1], score_fn)
            if M[i,j] == M[i-1,j-1]+s: state = 0
            elif M[i,j] == X[i-1,j-1]+s: state = 1
            else: state = 2
            aln1.insert(0, prof1[i-1])
            aln2.insert(0, prof2[j-1])
            i -= 1; j -= 1
        elif state == 1:  # X = gap dans prof2
            aln1.insert(0, prof1[i-1])
            aln2.insert(0, {"-":1.0})  # colonne = 100% gaps
            i -= 1
            if X[i+1,j] == M[i,j]-(gap_open+gap_extend): state = 0
            else: state = 1
        else:  # Y = gap dans prof1
            aln1.insert(0, {"-":1.0})
            aln2.insert(0, prof2[j-1])
            j -= 1
            if Y[i,j+1] == M[i,j]-(gap_open+gap_extend): state = 0
            else: state = 2
    return aln1, aln2, M[n,m]  # retourner aussi le score final

###############################
# Alignement progressif MSA   #
###############################

def progressive_align_from_tree(node, sequences, score_fn, gap_open=15, gap_extend=6, alphabet="ACGT-"):
    """
    Alignement progressif basé sur un arbre guide.
    
    node : noeud de l'arbre (scipy to_tree)
    sequences : list[SeqRecord]
    score_fn : fonction de score caractère-caractère
    gap_open, gap_extend : pénalités de gap

    return : 
        msa (list of str) : Alignement multiple 
        ids (list of str) total des identifiants de séquences
    """
    if node.is_leaf():
        return [str(sequences[node.id].seq)], [sequences[node.id].id]
    left_msa, left_ids = progressive_align_from_tree(node.left, sequences, score_fn, gap_open, gap_extend, alphabet)
    right_msa, right_ids = progressive_align_from_tree(node.right, sequences, score_fn, gap_open, gap_extend, alphabet)
    # Construire profils gauche/droite
    prof_left = profil.build_profile(left_msa, alphabet)
    prof_right = profil.build_profile(right_msa, alphabet)
    # Aligner profil–profil
    aln_left, aln_right, score = needleman_wunsch_profiles(prof_left, prof_right, score_fn, gap_open=gap_open, gap_extend=gap_extend)
    # Fusionner en un seul MSA
    new_profile, new_msa = profil.merge_profiles(aln_left, left_msa, aln_right, right_msa)
    return new_msa, left_ids + right_ids


###############################
# Écriture au format Clustal   #
###############################
def write_clustal(msa, ids, filename, protein, width=60):
    """
    Écrit un alignement multiple au format Clustal .aln

    msa (list of str) : Alignement multiple 
    ids (list of str) Identifiants des séquences.
    Nom du fichier de sortie.
    protein ( boolean ) : Si True, utilise une table de conservation protéique.
    """
    n = len(msa)
    L = len(msa[0])

    # Table de conservation pour protéines (simplifiée)
    conserved_groups = [
        set("STA"), set("NEQK"), set("NHQK"), set("NDEQ"),
        set("QHRK"), set("MILV"), set("MILF"), set("HY"),
        set("FYW")
    ]

    def consensus_column(col):
        chars = [seq[col] for seq in msa]
        if all(c == chars[0] for c in chars):
            return "*"
        if protein:
            # Vérifier si toutes les aa appartiennent à un même groupe conservé
            for group in conserved_groups:
                if all(c in group or c == "-" for c in chars):
                    return ":"
        return " "

    with open(filename, "w") as f:
        f.write("CLUSTAL W multiple sequence alignment\n\n")

        # Parcours par blocs
        for start in range(0, L, width):
            end = min(L, start + width)

            # Écrire chaque séquence du bloc
            for i in range(n):
                f.write(f"{ids[i]:15} {msa[i][start:end]}\n")

            # Ligne de consensus
            consensus = "".join(consensus_column(c) for c in range(start, end))

            f.write(" " * 16 + consensus + "\n\n")

