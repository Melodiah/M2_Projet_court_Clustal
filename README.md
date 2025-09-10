# M2 Project - Program Reconstruction Using the Clustal Method

This is an **educational mini-project** based on the multiple sequence alignment method described by Higgins & Sharp (1988).
It is **not a full implementation of Clustal**, but a simplified version created **for educational purposes only**.  

The program can be run using two different modes:

- **DNA mode**: with a scoring scheme of +1 for a match, -1 for a transversion, and -2 for other substitutions.  
- **Protein mode**: using the BLOSUM62 scoring matrix.  

---

## Installation

Clone the repository:

```bash
git clone https://github.com/Melodiah/M2_projet_Clustal
```

Install dependencies with UV : 

```bash
uv sync
```

## Usage

You can run the program in either DNA or Protein mode.

### DNA mode
 ```bash
uv run src/main.py --input data_test/DNA/ --output outdir_DNA/ --mode DNA
```

### Protein mode

```bash
uv run src/main.py --input data_test/protein/ --output outdir_protein/ --mode protein
```
