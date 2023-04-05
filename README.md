# Random Motif Search Algorith

This program is a simple implementation of the random motif search algorithm

## Usage

```python
rnd_motif_search.py <number_of_iterations> 
```

## Output

- Plot of the best score found in each iteration.
- Text file with the best motifs found in each iteration and the consensus motif highlighted in   red.
- To open text file in terminal: cat <filename> | less -R

## Dependencies

- matplotlib
- seaborn

# Gibbs Sampler Algorithm

This program is a simple implementation of the Gibbs sampler algorithm for motif finding. It uses ANSI colors to print the best motif found in the sequence.

## Usage

```python
gibbs_sampler.py : with default parameters
gibbs_sampler.py <N> <n_iter> : with custom parameters
  - N: number of iterations inside the Gibbs sampler for sequence selection
  - n_iter: number of iterations for repeating the Gibbs sampler
```

## Output

- Plot of the best score found in each iteration.

- Text file with the best motifs found in each iteration highlighted in red.
- To open text file in terminal: cat <filename> | less -R

## Dependencies

- matplotlib
- seaborn