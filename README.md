# BioInformatics
BioInformomatics Project
# Sequence Alignment, HMM Training, and Random Sequence Generation

This repository provides a comprehensive solution for sequence alignment, Hidden Markov Model (HMM) training, and random sequence generation, applied to biological sequence data. The project demonstrates multiple sequence alignment, scoring systems, and probabilistic modeling using HMMs. Additionally, it includes functionality for generating random nucleotide sequences to test the alignment algorithms.

## Overview

The repository features the following components:

1. **Sequence Alignment**: Performs pairwise and multiple sequence alignments using dynamic programming algorithms.
2. **Hidden Markov Model (HMM) Training**: Trains an HMM based on biological and randomly generated sequence datasets to calculate transition and emission probabilities.
3. **Random Sequence Generation**: Generates random nucleotide sequences with a specified length and nucleotide alphabet (A, C, G, T) for testing and comparison.
4. **Scoring Systems**: Calculates alignment scores using match, mismatch, and gap penalties.
5. **Normalization**: Normalizes transition and emission probabilities to ensure they sum to 1 for probabilistic accuracy.

## Features

- **Multiple Sequence Alignment (MSA)**: Identifies conserved positions across sequences and computes alignment scores.
- **HMM Training**: Builds a Hidden Markov Model (HMM) from aligned sequences, calculating the probabilities of sequence transitions and emissions.
- **Random Sequence Generation**: Generates random sequences for testing the alignment and HMM training process.
- **Dynamic Programming for Alignment Scoring**: Uses dynamic programming to calculate the optimal alignment score and paths for sequence pairs.
- **Normalization of Probabilities**: Normalizes the HMM's transition and emission probabilities to maintain a valid probabilistic model.

## Requirements

- Python 3.x
- Standard Python libraries (`collections`, `pprint`)

