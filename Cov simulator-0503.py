import random
from Bio import Align
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import nt_search

def load_sequence_from_fasta(file_path):
    # Load sequence from a local FASTA file
    try:
        with open(file_path, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                return str(record.seq)
    except FileNotFoundError:
        print("File not found.")
        return None

def load_sequences_from_fasta_list(file_paths):
    sequences = []
    for file_path in file_paths:
        try:
            with open(file_path, "r") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    sequences.append((str(record.seq), file_path))  # Store both sequence and file name
        except FileNotFoundError:
            print("File not found:", file_path)
    return sequences

def replicate_sequence(sequence, num_replications, substitution_rate, insertion_deletion_rate):
    replicated_sequences = []
    mutated_replicates = 0
    non_synonymous_replicates = 0
    
    for _ in range(num_replications):
        replicated_sequence = list(sequence)
        mutated = False
        i = 0
        while i < len(replicated_sequence):
            if random.random() < substitution_rate:
                # Substitution mutation
                replicated_sequence[i] = mutate_base(replicated_sequence[i])
                mutated = True
                i += 1
            elif random.random() < insertion_deletion_rate:
                if random.choice(['insert', 'delete']) == 'insert':
                    # Insertion mutation
                    new_base = random.choice(['A', 'T', 'C', 'G'])
                    replicated_sequence.insert(i, new_base)
                    mutated = True
                    i += 1  # Move past the inserted base
                else:
                    # Deletion mutation
                    deletion_length = random.randint(1, min(3, len(replicated_sequence) - i))
                    del replicated_sequence[i:i + deletion_length]
                    mutated = True
                    # No increment of i, will check next base in next iteration
            else:
                i += 1

        replicated_sequence = ''.join(replicated_sequence)
        
        # Check if the length of the sequence is a multiple of three
        if len(replicated_sequence) % 3 == 0:
            replicated_sequences.append(replicated_sequence)
            if mutated:
                mutated_replicates += 1
                original_protein_sequence = Seq(str(sequence)).translate(to_stop=True)
                mutated_protein_sequence = Seq(replicated_sequence).translate(to_stop=True)
                if original_protein_sequence != mutated_protein_sequence:
                    non_synonymous_replicates += 1
    
    return replicated_sequences, mutated_replicates, non_synonymous_replicates

def mutate_base(base):
    mutation_prob = random.random()
    if base == 'A':
        if mutation_prob < 0.3:  # Transition probability
            return 'G'
        elif mutation_prob < 0.6:  # Transversion probability
            return random.choice(['C', 'T'])
        else:
            return 'A'  # No mutation
    elif base == 'T':
        if mutation_prob < 0.3:  # Transition probability
            return 'C'
        elif mutation_prob < 0.6:  # Transversion probability
            return random.choice(['A', 'G'])
        else:
            return 'T'  # No mutation
    elif base == 'C':
        if mutation_prob < 0.3:  # Transition probability
            return 'T'
        elif mutation_prob < 0.6:  # Transversion probability
            return random.choice(['A', 'G'])
        else:
            return 'C'  # No mutation
    elif base == 'G':
        if mutation_prob < 0.3:  # Transition probability
            return 'A'
        elif mutation_prob < 0.6:  # Transversion probability
            return random.choice(['C', 'T'])
        else:
            return 'G'  # No mutation
    else:
        return base  # For any other base (e.g., N), return as is

def calculate_mutation_rate_per_cycle(mutation_rate_per_nucleotide, sequence_length):
    return mutation_rate_per_nucleotide * sequence_length

def calculate_sequence_similarity(seq1, seq2):
    aligner = Align.PairwiseAligner()
    alignments = aligner.align(seq1, seq2)
    best_alignment = alignments[0]  # Assuming only one alignment is returned
    max_score = best_alignment.score
    seq_length = max(len(seq1), len(seq2))
    similarity_percentage = (max_score / seq_length) * 100
    return similarity_percentage

def write_fasta_file_with_cycles(output_file_path, best_replicates_for_cycles):
    records = []
    for cycle_num, best_replicate in best_replicates_for_cycles.items():
        description = f"Cycle_{cycle_num}_Best_replicate"
        record = SeqIO.SeqRecord(Seq(best_replicate), id=description, description="")
        records.append(record)

    with open(output_file_path, "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")

def simulate_multiple_cycles(input_sequence, other_sequences, num_cycles, num_replications_per_cycle, substitution_rate, insertion_deletion_rate, output_file_path):
    best_similarity_overall = 0
    best_replicate_overall = None
    best_coronavirus_gene_overall = None
    current_input_sequence = input_sequence

    best_replicates_for_cycles = {}  # Dictionary to store best replicates for every 10 cycles
    similarities_per_cycle = []  # List to store similarities of best replicate for each cycle

    for cycle in range(num_cycles):
        replicated_sequences_cycle = []
        mutated_replicates_cycle = 0
        non_synonymous_replicates_cycle = 0

        for _ in range(num_replications_per_cycle):
            replicated_sequences, mutated_replicates, non_synonymous_replicates = replicate_sequence(current_input_sequence, 1, substitution_rate, insertion_deletion_rate)
            replicated_sequences_cycle.extend(replicated_sequences)
            mutated_replicates_cycle += mutated_replicates
            non_synonymous_replicates_cycle += non_synonymous_replicates

        percentage_mutated_cycle = (mutated_replicates_cycle / num_replications_per_cycle) * 100
        percentage_non_synonymous_cycle = (non_synonymous_replicates_cycle / mutated_replicates_cycle) * 100 if mutated_replicates_cycle != 0 else 0

        filtered_replicates = [rep for rep in replicated_sequences_cycle if Seq(rep).translate(to_stop=True) != Seq(current_input_sequence).translate(to_stop=True)]
        filtered_replicates.append(current_input_sequence)  # Add the current input sequence (best from previous cycle)

        best_similarity_cycle = 0
        best_replicate_cycle = None
        best_coronavirus_gene_cycle = None

        for replicate in filtered_replicates:
            for other_seq_tuple in other_sequences:
                other_seq, gene_file_name = other_seq_tuple
                similarity = calculate_sequence_similarity(Seq(replicate).translate(to_stop=True), Seq(other_seq).translate(to_stop=True))
                if similarity > best_similarity_cycle:
                    best_similarity_cycle = similarity
                    best_replicate_cycle = replicate
                    best_coronavirus_gene_cycle = gene_file_name

        if best_similarity_cycle > best_similarity_overall:
            best_similarity_overall = best_similarity_cycle
            best_replicate_overall = best_replicate_cycle
            best_coronavirus_gene_overall = best_coronavirus_gene_cycle

        current_input_sequence = best_replicate_cycle  # Set the best replicate as the input for the next cycle

        print(f"Cycle {cycle + 1}:")
        print("Percentage of mutated replicates:", percentage_mutated_cycle)
        print("Percentage of replicates with non-synonymous changes:", percentage_non_synonymous_cycle)
        print("Best replicate with similarity:", best_similarity_cycle)
        print()

        similarities_per_cycle.append(best_similarity_cycle)  # Record the similarity of the best replicate for this cycle

        # Store the best replicate for every 10 cycles
        if (cycle + 1) % 10 == 0:
            cycle_num = cycle + 1
            best_replicates_for_cycles[cycle_num] = best_replicate_cycle

    # Write the recorded similarities to a text file
    output_similarity_file_path = "similarities.txt"

    with open(output_similarity_file_path, "w") as similarity_file:
        similarity_file.write("Cycle\tSimilarity\n")
        for cycle_num, similarity in enumerate(similarities_per_cycle, start=1):
            similarity_file.write(f"{cycle_num}\t{similarity}\n")

    print("Overall best replicate with similarity:", best_similarity_overall)
    if best_replicate_overall:
        print("Nucleotide sequence of the overall best replicate:")
        print(best_replicate_overall)
        print("Best coronavirus spike gene:", best_coronavirus_gene_overall)
    else:
        print("No overall best replicate found.")

    # Write the best replicates for every 10 cycles into a single FASTA file
    write_fasta_file_with_cycles(output_file_path, best_replicates_for_cycles)

# Example usage
file_path = "F:\CoV-mutation-simulator-ongoing4\Pangolin.fasta"  # Input
other_sequences_file_paths = ["F:\CoV-mutation-simulator-ongoing4\SARS2.fasta"]   # Similarity reference
input_sequence = load_sequence_from_fasta(file_path)
other_sequences = load_sequences_from_fasta_list(other_sequences_file_paths)

if input_sequence is None:
    print("Failed to load input sequence.")
elif not other_sequences:
    print("Failed to load other sequences.")
else:
    substitution_frequency = 24  # Substitution frequency relative to insertion/deletion
    substitution_rate = 3.76e-5
    insertion_deletion_rate = substitution_rate / substitution_frequency
    num_cycles = 1000
    num_replications_per_cycle = 10000
    output_file_path = "output.fasta"  # Specify the output file path
    simulate_multiple_cycles(input_sequence, other_sequences, num_cycles, num_replications_per_cycle,
                             substitution_rate, insertion_deletion_rate, output_file_path)
