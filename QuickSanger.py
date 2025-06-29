#This is my Aliview pipeline for the most efficient cutting and merging of Aliview files 
#The way that this works is a simple pipeline 

#Lets break down the basic pipeline for this 

#ALIGNMENT 

#First---The input reads from files that are in FASTA format coupled with corresponding ab1 files
#Runs can be in single or batch for alignment and same for merging
#merging and alignment are not coupled in the pipeline, but it is simple to batch merge the created folder

#The two original FASTA files that are pretrimmed are aligned by reverse complement of the reverse primer 
#Then all gaps and mismatches are identified 

#The software then computes an index aligned between the ab1 and the FASTA sequence 
#This is done by "finding" the first 20 bp's of the original FASTA file in the ab1 and 
#that index becomes the start position

#Traversing both indices, gaps and mismatches are identified
#As a part of the "correct_and_clean_alignment" method--there is a confidence score threshold you can adjust 
#By comparing confidence scores of reads -- Alanview either makes a substiution to fill gaps or correct mismatches or
# a deletion to fix misreads

#ANY extremely bad chromatography reads should be taken out of the batch by manual sort

#Finally the corrections are made and the aligned file(s) are exported into a locally created aligned_files folder 

#MERGE 

#Also runs in either single or batch
#Just merges the two sequences after alignment and stores in new folder

#GENERAL 

#File input formatting for batch should be as follows 
#Folder 1            Folder 2 
#Sample_1.FWD.fas    Sample_1.FWD.ab1
#Sample_1.REV.fas    Sample_2.REV.ab1
#Sample_2.FWD.fas    Sample_2.FWD.ab1
#...                 ....

#Anything that is unresolvable by the machine will be exported to a "needs manual check" folder 
#for manual analysis 

#**An example of a fasta formated version of this alignment that should be saved into the aligned folder would be:

"""
> HB20_COX-LCO1490
GGTGACCAAAAAATCAAAATAAATGTTGATATAAAATTGGATCTCCACCTCCTATTGGATCAAAAAATGAAGTATTAAAATTTCGATCAAATAATAATATTGTAATTGCTCCAGCTAGAACTGGAAGAGATATAATTAATAAGATTGCTGTAATAAATACTGATCATGGAAATAAAGAAATTTGGTCATAATTTAATGAAAAATTTTTTATTATTATAATTGTAACTATTAAATTTAATGATCCTATAATTGATGAGATTCCAGATATATGTAAGGAGAAAATTGCAAAATCAACTGAAGGGGATGAATGATATAAATATGCAGATAATGGAGGATAAACTGTTCATCCTGTTCCTGGTCTTGGATAAAATAAATTTCTTAATAATAATATAAATAATGAAGGAGGAAGCAATCAGAATCTAATATTATTTATTCGAGGAAATGCTATATCTGGAGATCCTAATATTAAAGGAATTAATCAATTTCCAAAACCTCCAATTAGGAATGGTATAACTATAAAAAAAATTATTAGAAATGCATGTCTTGTTACAACTGTATTATAAATTTGATCATTATTAATTCATGAACCGGGGGATCTTAATTCTATGCGAACAATCAACCTTATTGATGAGCCTAATATTC------------------------------------------------------
> HB20_COX-HCO2198
------------------------------------------CTCCACCTCCTATTGGATCAAAAAATGAAGTATTAAAATTTCGATCAAATAATAATATTGTAATTGCTCCAGCTAGAACTGGAAGAGATATAATTAATAAGATTGCTGTAATAAATACTGATCATGGAAATAAAGAAATTTGGTCATAATTTAATGAAAAATTTTTTATTATTATAATTGTAACTATTAAATTTAATGATCCTATAATTGATGAGATTCCAGATATATGTAAGGAGAAAATTGCAAAATCAACTGAAGGGGATGAATGATATAAATATGCAGATAATGGAGGATAAACTGTTCATCCTGTTCCTGGTCTTGGATAAAATAAATTTCTTAATAATAATATAAATAATGAAGGAGGAAGCAATCAGAATCTAATATTATTTATTCGAGGAAATGCTATATCTGGAGATCCTAATATTAAAGGAATTAATCAATTTCCAAAACCTCCAATTAGGAATGGTATAACTATAAAAAAAATTATTAGAAATGCATGTCTTGTTACAACTGTATTATAAATTTGATCATTATTAATTCATGAACCGGGGGATCTTAATTCTATGCGAACAATCAACCTTATTGATGAGCCTAATATTCCTGATCATAAAGCTAATAAAATATATAAAATTCCAATATCTTTATGATTTTGTT

"""


#Imports
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.pairwise2 import format_alignment
import os
from pathlib import Path
from tkinter import Tk, filedialog

#Loads the sequence 
def load_fasta_sequence(file_path):
    record = SeqIO.read(file_path, "fasta")
    return str(record.seq), record.id

#Used for both REV.fas and REV.ab1
def reverse_complement_seq(seq):
    return str(Seq(seq).reverse_complement())

#adjust alignment weights if needed but this works well. 
def align_sequences(seq1, seq2):
    alignments = pairwise2.align.localms(seq1, seq2, 2, -1, -0.5, -0.1)
    best_alignment = alignments[0]
    return best_alignment

#Shows alignment score in terminal
def display_alignment(alignment):
    seqA, seqB, score, start, end = alignment
    print(f"\n Alignment Score: {score}\n")
    

#Creates local folder for aligned files
def save_alignment_as_fasta(seqA, seqB, id1, id2, output_folder="aligned_files"):
    os.makedirs(output_folder, exist_ok=True)
    output_path = os.path.join(output_folder, f"{id1}_{id2}_aligned.fasta")
    with open(output_path, "w") as f:
        f.write(f">{id1}\n{seqA}\n")
        f.write(f">{id2}\n{seqB}\n")
    print(f"Aligned FASTA saved to: {output_path}")

#Read and differentiate a dataset of confidence scores correlating to ab1 indices
def read_ab1_confidence_scores(ab1_path, is_reverse=False):
    record = SeqIO.read(ab1_path, "abi")
    sequence = str(record.seq)
    confidence_scores = record.letter_annotations['phred_quality']

    if is_reverse:
        sequence = str(Seq(sequence).reverse_complement())
        confidence_scores = confidence_scores[::-1]

    return sequence, confidence_scores

#Uses first 20 FASTA to find offset reigon of the ab1
def find_chromatogram_start(full_seq, trimmed_seq):
    start_index = full_seq.find(trimmed_seq[:20])
    if start_index == -1:
        raise ValueError("Trimmed sequence not found in chromatogram sequence.")
    return start_index

#Builds map fo indices for direct ab1 to FASTA comparison 
def build_index_map(aligned_seq):
    index_map = []
    chrom_index = 0

    for base in aligned_seq:
        if base == '-':
            index_map.append(None)
        else:
            index_map.append(chrom_index)
            chrom_index += 1

    return index_map


#Large alignment method with corrections
def correct_and_clean_alignment(alignment, conf_seq1, conf_seq2, offset_A, offset_B, confidence_threshold=20):
    seqA, seqB, score, start, end = alignment
    seqA = list(seqA)
    seqB = list(seqB)

    index_map_A = build_index_map(seqA)
    index_map_B = build_index_map(seqB)

    # Find overlap boundaries (first and last aligned positions) so only checking for gaps within region and not overhangs 
    first_aligned_pos = next(i for i, (a, b) in enumerate(zip(seqA, seqB)) if a != '-' and b != '-')
    last_aligned_pos = len(seqA) - 1 - next(i for i, (a, b) in enumerate(zip(reversed(seqA), reversed(seqB))) if a != '-' and b != '-')

    for index in range(last_aligned_pos, first_aligned_pos - 1, -1):
        baseA = seqA[index]
        baseB = seqB[index]

        if baseA == '-' and baseB == '-':
            continue

        chrom_index_A = index_map_A[index]
        chrom_index_B = index_map_B[index]

        chrom_index_A = chrom_index_A + offset_A if chrom_index_A is not None else None
        chrom_index_B = chrom_index_B + offset_B if chrom_index_B is not None else None

        confA = conf_seq1[chrom_index_A] if chrom_index_A is not None and chrom_index_A < len(conf_seq1) else 0
        confB = conf_seq2[chrom_index_B] if chrom_index_B is not None and chrom_index_B < len(conf_seq2) else 0

        # mismatches
        if baseA != baseB and baseA != '-' and baseB != '-':
            if confA > confB:
                seqB[index] = baseA
            else:
                seqA[index] = baseB

        # gaps where the other sequence is confident 
        elif baseA == '-' and baseB != '-':
            if confB >= confidence_threshold:
                seqA[index] = baseB
            else:
                # Low confidence: delete both
                del seqA[index]
                del seqB[index]

        elif baseB == '-' and baseA != '-':
            if confA >= confidence_threshold:
                seqB[index] = baseA
            else:
                # Low confidence: delete both
                del seqA[index]
                del seqB[index]

    corrected_seqA = ''.join(seqA)
    corrected_seqB = ''.join(seqB)

    #Recalculate for adjusted positions (updated code)
    first_aligned_pos_2 = next(i for i, (a, b) in enumerate(zip(corrected_seqA, corrected_seqB)) if a != '-' and b != '-')
    last_aligned_pos_2 = len(seqA) - 1 - next(i for i, (a, b) in enumerate(zip(reversed(corrected_seqA), reversed(corrected_seqB))) if a != '-' and b != '-')

    # Only check for final gaps within the overlap region if manual check needed
    overlapA = corrected_seqA[first_aligned_pos_2:last_aligned_pos_2+1]
    overlapB = corrected_seqB[first_aligned_pos_2:last_aligned_pos_2+1]
    if '-' in overlapA or '-' in overlapB:
        os.makedirs('manual_check_required', exist_ok=True)
        file_index = len(os.listdir('manual_check_required')) // 2 + 1  # count pairs

        # Find gap indices in overlap region (relative to overlap)
        gap_indices_A = [i for i, base in enumerate(overlapA) if base == '-']
        gap_indices_B = [i for i, base in enumerate(overlapB) if base == '-']

        # Convert to indices relative to full sequence
        gap_indices_A_full = [first_aligned_pos + i for i in gap_indices_A]
        gap_indices_B_full = [first_aligned_pos + i for i in gap_indices_B]

        with open(f'manual_check_required/sample_{file_index}_seq.fasta', 'w') as f:
           f.write(f">Corrected_SeqA_Sample_{file_index}\n{corrected_seqA}\n")
           f.write(f">Corrected_SeqB_Sample_{file_index}\n{corrected_seqB}\n")

        print(f"\nWarning: Gaps remain in overlap. Sequences saved to manual_check_required/sample_{file_index}_seqA.fasta and seqB.fasta")
        if gap_indices_A_full:
            print(f"Gap(s) in seqA at indices: {gap_indices_A_full}")
        if gap_indices_B_full:
            print(f"Gap(s) in seqB at indices: {gap_indices_B_full}")

    return corrected_seqA, corrected_seqB

def merge_aligned_sequences(seqA, seqB):
    # Find overlap boundaries
    def find_overlap_bounds(a, b):
        left = next(i for i, (x, y) in enumerate(zip(a, b)) if x != '-' and y != '-')
        right = len(a) - 1 - next(i for i, (x, y) in enumerate(zip(reversed(a), reversed(b))) if x != '-' and y != '-')
        return left, right

    left, right = find_overlap_bounds(seqA, seqB)

    # Overhangs
    prefix = ''
    suffix = ''
    if left > 0:
        # Use non-gap overhangs from either sequence
        prefix = ''.join([a if a != '-' else b for a, b in zip(seqA[:left], seqB[:left])])
    if right < len(seqA) - 1:
        suffix = ''.join([a if a != '-' else b for a, b in zip(seqA[right+1:], seqB[right+1:])])

    # Merge overlap
    overlapA = seqA[left:right+1]
    overlapB = seqB[left:right+1]
    merged_overlap = []
    manual_check_needed = False
    for a, b in zip(overlapA, overlapB):
        if a == b:
            merged_overlap.append(a.replace('-', ''))  # remove gap if both are gap, will be ''
        elif a == '-':
            merged_overlap.append(b)
        elif b == '-':
            merged_overlap.append(a)
        else:
            # mismatch in overlap
            manual_check_needed = True
            merged_overlap.append(a)  # arbitrary, but will trigger manual check

    merged_seq = prefix + ''.join(merged_overlap) + suffix

    # Only trigger manual check if there are gaps or mismatches in the overlap
    if manual_check_needed or '-' in ''.join(merged_overlap):
        print("Warning: Overlap region contains mismatches or gaps. Saving for manual check.")
        os.makedirs('manual_check_required', exist_ok=True)
        file_index = len(os.listdir('manual_check_required')) // 2 + 1

        with open(f'manual_check_required/merge_issue_seqA_{file_index}.fasta', 'w') as f:
            f.write(f">Merge_Issue_SeqA_{file_index}\n{seqA}\n")

        with open(f'manual_check_required/merge_issue_seqB_{file_index}.fasta', 'w') as f:
            f.write(f">Merge_Issue_SeqB_{file_index}\n{seqB}\n")

        return None  # Skip merging for now

    return merged_seq


def save_merged_sequence(merged_seq, sample_id, output_folder="merged"):
    os.makedirs(output_folder, exist_ok=True)
    output_path = os.path.join(output_folder, f"{sample_id}_merged.fasta")
    with open(output_path, "w") as f:
        f.write(f">{sample_id}_merged\n{merged_seq}\n")
    print(f"Merged sequence saved to: {output_path}")

def select_file(prompt="Select a file"):
    root = Tk()
    root.withdraw()
    file_path = filedialog.askopenfilename(title=prompt)
    if not file_path:
        raise FileNotFoundError("No file selected.")
    return file_path


def select_folder(prompt="Select a folder"):
    root = Tk()
    root.withdraw()
    folder_path = filedialog.askdirectory(title=prompt)
    if not folder_path:
        raise FileNotFoundError("No folder selected.")
    return folder_path


def process_sample(file1, file2, ab1_file1, ab1_file2):
    seq1, id1 = load_fasta_sequence(file1)
    seq2, id2 = load_fasta_sequence(file2)

    seq2_rc = reverse_complement_seq(seq2)

    alignment = align_sequences(seq1, seq2_rc)

    print(f"\nProcessing: {os.path.basename(file1)} and {os.path.basename(file2)}")
    print("\nInitial Alignment:")
    display_alignment(alignment)

    chrom_seq1, conf_seq1 = read_ab1_confidence_scores(ab1_file1, is_reverse=False)
    chrom_seq2, conf_seq2 = read_ab1_confidence_scores(ab1_file2, is_reverse=True)

    offset_A = find_chromatogram_start(chrom_seq1, seq1)
    offset_B = find_chromatogram_start(chrom_seq2, reverse_complement_seq(seq2))

    corrected_seqA, corrected_seqB = correct_and_clean_alignment(alignment, conf_seq1, conf_seq2, offset_A, offset_B)
    """
    print("\nCorrected Alignment:")
    print(corrected_seqA)
    print(corrected_seqB)
    """
    save_alignment_as_fasta(corrected_seqA, corrected_seqB, id1, id2)


def main():
    mode = input("Select mode: [1] Single Run or [2] Batch Run or [3] Merge: ").strip()

    if mode == '1':
        print("Select the FORWARD FASTA file:")
        file1 = select_file()

        print("Select the REVERSE FASTA file:")
        file2 = select_file()

        print("Select the FORWARD AB1 chromatogram file:")
        ab1_file1 = select_file()

        print("Select the REVERSE AB1 chromatogram file:")
        ab1_file2 = select_file()

        process_sample(file1, file2, ab1_file1, ab1_file2)
    
    elif mode == '2':
        print("Select folder with FASTA files:")
        fasta_folder = select_folder()

        print("Select folder with AB1 files:")
        ab1_folder = select_folder()

        fasta_files = sorted([os.path.join(fasta_folder, f) for f in os.listdir(fasta_folder)])
        ab1_files = sorted([os.path.join(ab1_folder, f) for f in os.listdir(ab1_folder)])

        if len(fasta_files) % 2 != 0 or len(ab1_files) % 2 != 0 or len(fasta_files) != len(ab1_files):
            print("Error: Files are not in expected counts or pairs. Ensure folders contain (FWD, REV, FWD, REV...) format.")
            return

        num_samples = len(fasta_files) // 2
        for i in range(num_samples):
            file1 = fasta_files[2 * i]
            file2 = fasta_files[2 * i + 1]
            ab1_file1 = ab1_files[2 * i]
            ab1_file2 = ab1_files[2 * i + 1]

            process_sample(file1, file2, ab1_file1, ab1_file2)

        print("\nBatch processing complete.")

    elif mode == '3':

        merge_mode = input("\nSelect merge mode:\n1: Single Merge\n2: Batch Merge\nEnter choice (1 or 2): ").strip()

        if merge_mode == '1':
            file_path = select_file("Select the aligned FASTA file to merge:")
            records = list(SeqIO.parse(file_path, "fasta"))

            if len(records) != 2:
                print("Error: Aligned FASTA must contain exactly two sequences.")
                return

            seqA = str(records[0].seq)
            seqB = str(records[1].seq)
            merged_seq = merge_aligned_sequences(seqA, seqB)

            if merged_seq:
                save_merged_sequence(merged_seq, records[0].id)

        elif merge_mode == '2':
            print("Select folder with aligned FASTA files:")
            aligned_folder = select_folder()

            aligned_files = sorted([os.path.join(aligned_folder, f) for f in os.listdir(aligned_folder)])

            for file_path in aligned_files:
                records = list(SeqIO.parse(file_path, "fasta"))

                if len(records) != 2:
                    print(f"Skipping {file_path}: does not contain exactly two sequences.")
                    continue

                seqA = str(records[0].seq)
                seqB = str(records[1].seq)
                merged_seq = merge_aligned_sequences(seqA, seqB)

                if merged_seq:
                    save_merged_sequence(merged_seq, records[0].id)

            print("\nBatch merging complete.")

        else:
            print("Invalid selection. Please restart and enter 1 or 2.")


    else:
        print("Invalid mode selected.")


if __name__ == "__main__":
    main()
