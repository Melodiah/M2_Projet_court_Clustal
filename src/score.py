"""
Module utilisé pour la définition, l'édition de matrices de scores 
dans le cadre d'un alignement de séquences.
"""

################################
#  Variable Globale            #
################################

DNA_TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}


#################################
#   Fonctions de Scores         #
#################################

def score_DNA(a,b):
    """
    Renvoie le score de substitution entre deux nucléotides
    selon le modèle simple : +1 si identiques, -1 si transition, -2 si transversion.
    """
    if a == b : # Match entre nucléotides
        return 1
    elif(a, b) in DNA_TRANSITIONS: # Transition
        return -1
    else: # Substitution transversion
        return -2


def protein_score(a: str, b: str, matrix):
    """
    Renvoie le score de substitution entre deux acides aminés
    selon une matrice de substitution donnée (BLOSUM62 par défaut).
    """
    if a == "-" or b == "-":
        # Les gaps sont gérés par gap_open et gap_extend,
        # mais on met une grosse pénalité si jamais on tombe dessus.
        return -5
    
    pair = (a, b)
    if pair not in matrix:
        pair = (b, a)
    return matrix.get(pair, -4)  # -4 si la paire n'existe pas


