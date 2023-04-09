#!/usr/bin/env python
# encoding=utf-8
# Created by: xfulop

"""
    This program is a simple implementation of the random motif search algorithm.

    Usage: 
        rnd_motif_search.py <number_of_iterations> 
            -- number of iterations to run the algorithm
            -- default is 500

    Output: 
        -- plot of the best score found in each iteration
        -- txt file with the best motifs found in each iteration highlighted in red
        -- to open txt file in terminal: cat <filename> | less -R


    Dependencies:
        - matplotlib
        - seaborn
"""

import random
import sys
import matplotlib.pyplot as plt
import seaborn as sns


# Define the sequences
Dna = [
    "TTACCTTAAC",
    "GATGTCTGTC",
    "CCGGCGTTAG",
    "CACTAACGAG",
    "CGTCAGAGGT"
]

def RandomizedMotifSearch(Dna, k, t):
    # randomly select k-mers from each string in Dna
    motifs = []
    iterations = 0
    for i in range(t):
        rand_index = random.randint(0, len(Dna[i])-k)
        motif = Dna[i][rand_index:rand_index+k]
        motifs.append(motif)
    
    # initialize best motifs as the randomly selected motifs
    best_motifs = motifs[:]
    
    while True:
        iterations += 1
        profile = GenerateProfileMatrix(motifs)
        motifs = GenerateMotifs(profile, Dna)

        
        if Score(motifs) < Score(best_motifs):
            best_motifs = motifs[:]
        else:
            # print best motifs highlighted in red using ANSI
            # escape sequences using original sequences
            all_best_mot_out = []
            for i in range(len(best_motifs)):
                best_mot_out = Dna[i].replace(best_motifs[i], 
                                              f'\033[91m{best_motifs[i]}\033[0m')

                # Merge append into list declaration to avoid using global variable 
                all_best_mot_out.append(best_mot_out)
                consensus_motif = consensus(best_motifs)
            return best_motifs, Score(best_motifs), iterations, all_best_mot_out, consensus_motif

def GenerateProfileMatrix(motifs):
    t = len(motifs)
    k = len(motifs[0])
    profile = {'A': [], 'C': [], 'G': [], 'T': []}
    
    for j in range(k):
        count = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        for i in range(t):
            count[motifs[i][j]] += 1
        for nucleotide in ['A', 'C', 'G', 'T']:
            profile[nucleotide].append(count[nucleotide] / t)
    
    return profile

def GenerateMotifs(profile, Dna):
    k = len(profile['A'])
    motifs = []
    for i in range(len(Dna)):
        motif = ProfileMostProbableKmer(Dna[i], k, profile)
        motifs.append(motif)
    return motifs

def ProfileMostProbableKmer(text, k, profile):
    max_prob = -1
    most_probable_kmer = text[:k]
    
    for i in range(len(text) - k + 1):
        kmer = text[i:i+k]
        prob = 1
        for j in range(k):
            prob *= profile[kmer[j]][j]
        if prob > max_prob:
            max_prob = prob
            most_probable_kmer = kmer
    
    return most_probable_kmer

def Score(motifs):
    k = len(motifs[0])
    score = 0
    for j in range(k):
        count = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        for i in range(len(motifs)):
            count[motifs[i][j]] += 1
        max_count = max(count.values())
        score += len(motifs) - max_count
    return score

def consensus(motifs):
    k = len(motifs[0])
    count = {}
    consensus = ""
    for j in range(k):
        count = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        for i in range(len(motifs)):
            count[motifs[i][j]] += 1
        m = max(count, key=count.get)
        consensus += m
    return consensus

def visualize_until_convergence(score_lst):
    # plot the score for each iteration
    sns.set_style('darkgrid')
    plt.figure(figsize=(15, 5))
    plt.plot(score_lst)  
    plt.xlabel('Iterations')
    plt.ylabel('Score')
    plt.title('Randomized Motif Search Convergence')
    # find the lowest score and the corresponding iteration and highlight it
    min_score = min(score_lst)
    min_score_index = score_lst.index(min_score)
    plt.scatter(min_score_index, min_score, c='red', s=100)
    # legend
    if min_score_index == 0:
        plt.legend(['Score', 'Lowest score found after first iteration'])
    else:
        plt.legend(['Score', 'Lowest score found after {min_score_index} iterations'.format(min_score_index=min_score_index)])
    # start the axis from 1
    plt.xlim(1, len(score_lst)+1)
    # plt.text(min_score_index, min_score, f'Lowest score: {min_score} found after {min_score_index} iterations')
    plt.savefig(f'randomized_motif_search_{n_times}_iterations.png')
    plt.show()
    


# function for repeating the algorithm n times (user input)
# if there is no user input, the default is 500
n_times = int(sys.argv[1]) if len(sys.argv) > 1 else 500
motif_lst = []
score_lst = []
all_best_mot_out_lst = []
consensus_lst = []


def RepeatRandomizedMotifSearch(n_times):
    consensus = ''
    consensus_lst = []
    for _ in range(n_times):
        motifs, score, iterations, all_best_mot_out, consensus = RandomizedMotifSearch(Dna, 4, 5)
        motif_lst.append(motifs)
        score_lst.append(score)
        all_best_mot_out_lst.append(all_best_mot_out)
        consensus_lst.append(consensus)

    # sort the data based on the score
    data = list(zip(score_lst, all_best_mot_out_lst, consensus_lst))
    data_sorted = sorted(data, key=lambda x: x[0])
    
    # write the sorted data to a file
    with open('rnd_mtf_sch_output.txt', 'w') as f:
        f.write('Results of the randomized motif search algorithm:\n')
        f.write('sorted by score:\n')
        f.write('-----------------------------------------\n')
        for idx, (score, all_best_mot_out, consensus) in enumerate(data_sorted):
            f.write(f'N = {idx+1} of {n_times}\n')
            f.write(f'Best score = {score}\n')
            f.write(f'Consensus motif: {consensus}\n')
            f.write('Best motifs found after iterations:\n')
            for motif in all_best_mot_out:
                f.write(f'{motif}\n')
            f.write('-----------------------------------------\n')

def main():
    print(f'Running the randomized motif search algorithm {n_times} times')
    RepeatRandomizedMotifSearch(n_times)
    # visualize the score for each iteration
    print('Visualizing the score for each iteration')
    visualize_until_convergence(score_lst)
    print('Done!')

if __name__ == '__main__':
    main()