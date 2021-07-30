def scan_clumps(text, k, L, t):
    # Break down message in k-mers.
    # kmers = [ k-mer1, k-mer2, k-mer3, ...]
    kmers = [text[i:(i+k)] for i in range(len(text)-k+1)]
    kmers_per_window = L-k+1

    # Get all the k-mers from the first window of size L and calculate frequency
    f_kmers = {}
    for _kmer in kmers[0:kmers_per_window]:
        f_kmers[_kmer] = f_kmers.get(_kmer, 0) + 1

    # Keep valid clumps as initial solution
    solution = {x: v for x, v in f_kmers.items() if v >= t}

    # Sliding window of size L across the list of kmers
    # For each iteration:
    # 1. Decrement the frequency of the k-mer leaving the window (i)
    # 2. Increment the frequency of the k-mer entering the window (i+kmers_per_window)
    # 3. If the frequency of the k-mer entering the window is greater than t, save as candidate solution
    for i in range(0, len(kmers) - kmers_per_window):
        old_kmer = kmers[i]
        new_kmer = kmers[i + kmers_per_window]

        f_kmers[old_kmer] -= 1
        f_kmers[new_kmer] = f_kmers.get(new_kmer, 0) + 1
        if (f_kmers[new_kmer] >= t):
            solution[new_kmer] = max(solution.get(new_kmer, 0), f_kmers[new_kmer])

    return(solution)


with open("vibrio_cholerae_genome.txt", "r") as input_file:
    vc_genome = ''.join(line.rstrip() for line in input_file).lower()
    input_file.close()

clumps = scan_clumps(vc_genome, 9, 500, 3)
print(len(clumps))

#ANS : 1062
