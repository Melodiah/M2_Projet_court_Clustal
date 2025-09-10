from collections import Counter, defaultdict

###############################
# Créations des profils MSA   #
###############################

def build_profile(msa, alphabet="ACGT-"):
    """
    Construit un profil à partir d'un alignement multiple (liste de séquences alignées).
    Chaque colonne est représentée par un dictionnaire {caractère: fréquence}.
    
    msa : list[str]  -> séquences alignées (même longueur)
    alphabet : str   -> caractères attendus (par défaut ADN + gap)
    
    return : list[dict]  -> profil (colonne par colonne)
    """
    profile = []
    L = len(msa[0])
    for col in range(L):
        counts = Counter(seq[col] for seq in msa)
        freqs = {a: counts.get(a, 0)/len(msa) for a in alphabet}
        profile.append(freqs)
    return profile

def score_columns(col1, col2, score_fn):
    """
    Calcule le score moyen entre deux colonnes de profils.
    
    col1, col2 : dict {caractère: fréquence}
    score_fn   : fonction donnant un score pour deux caractères
    
    return : float
    """
    total = 0.0
    for a, fa in col1.items():
        for b, fb in col2.items():
            total += fa * fb * score_fn(a, b)
    return total

def merge_profiles(aln1, msa1, aln2, msa2):
    """
    Fusionne deux profils alignés en un seul, en propageant les gaps
    dans toutes les séquences du MSA.
    
    aln1, aln2 : profils alignés (list de colonnes dict)
    msa1, msa2 : séquences correspondantes (list[str])
    
    return : (nouveau_profil, nouveau_MSA)
    """

    new_msa = []
    new_profile = []

    # --- 1. Reconstruire les séquences du premier groupe ---
    for s in msa1:
        new_seq = ""
        idx = 0
        for col in aln1:
            if "-" in col and col["-"] == 1.0 and len(col) == 1:
                # colonne de gap -> insérer un gap dans cette séquence
                new_seq += "-"
            else:
                # reprendre le prochain caractère de la séquence originale
                new_seq += s[idx]
                idx += 1
        new_msa.append(new_seq)

    # --- 2. Reconstruire les séquences du second groupe ---
    for s in msa2:
        new_seq = ""
        idx = 0
        for col in aln2:
            if "-" in col and col["-"] == 1.0 and len(col) == 1:
                new_seq += "-"
            else:
                new_seq += s[idx]
                idx += 1
        new_msa.append(new_seq)

    # --- 3. Construire le nouveau profil ---
    for c1, c2 in zip(aln1, aln2):
        merged = defaultdict(float)
        for k, v in c1.items():
            merged[k] += v
        for k, v in c2.items():
            merged[k] += v
        total = sum(merged.values())
        for k in merged:
            merged[k] /= total
        new_profile.append(dict(merged))

    # --- 4. Vérification ---
    L = len(new_msa[0])
    assert all(len(seq) == L for seq in new_msa), "Erreur : longueurs incohérentes"

    return new_profile, new_msa