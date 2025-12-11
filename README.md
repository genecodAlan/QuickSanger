# QuickSanger!
Raw Sanger sequence alignment and merging. Auto gap and mismatch corrections based on AB1 chromatogram file confidence scores. Features batch and single alignment, as well as batch or single merge. Bioinformatics tool for rapid Sanger analysis---automatic alternative to Aliview---AlanView :)

## Features 🛠️
- Single or batch mode for alignment and merging

- Auto reverse complementing of reverse reads

 Local alignment using Biopython

Auto-correction of:
- Mismatches
- Gaps (based on chromatogram confidence)
- AB1 confidence score-based corrections using customizable thresholds

⚠️ Unresolvable cases sent to manual_check_required/ for manual curation

# Pipeline Overview
## 1. Input
Input must include:

Two trimmed FASTA files (forward and reverse)

Their corresponding AB1 files

Modes:

Single: Manually select files one by one

Batch: Provide two folders (FASTA and AB1), files are assumed to be ordered as:


Folder 1 (FASTA)         | Folder 2 (AB1)
-------------------------|--------------------------
Sample_1.FWD.fas         | Sample_1.FWD.ab1
Sample_1.REV.fas         | Sample_1.REV.ab1
Sample_2.FWD.fas         | Sample_2.FWD.ab1
Sample_2.REV.fas         | Sample_2.REV.ab1
...                      | ...

**use the sample folders to try for yourself**

## 2. Alignment
Reverse complements the reverse primer read

Performs local alignment using pairwise2

Identifies:

Mismatches

Gaps

Calculates the start index of the FASTA read within AB1 using the first 20 bp for mapping

Scores each base with Phred quality scores from AB1:
![image](https://github.com/user-attachments/assets/1a45423d-c6d2-40d6-80c8-a099ab80ce38)


Uses higher-scoring base to correct mismatch

Deletes gaps from low-quality reads

Skips unresolved regions or flags for manual checking

## 3. Output

aligned_files/: Corrected and aligned FASTA files

merged/: Merged consensus sequences

manual_check_required/: Problematic sequences needing manual inspection

Example Output FASTA (from aligned_files/)

![image](https://github.com/user-attachments/assets/6579ab3f-e5f8-482a-bc05-de4d1555fd23)


> HB20_COX-LCO1490
GGTGACCAAAAAATCAAAATAAATGTTGATATAAAATTGGATCTCCACCTCCTATTGGATCAAAAAATGAAGTATT...
> HB20_COX-HCO2198
------------------------------------------CTCCACCTCCTATTGGATCAAAAAATGAAGTATT...

# Notes
Trim poor-quality AB1 reads before running

Confidence threshold in corrections is customizable (confidence_threshold=20 by default)

Alignment and merging are decoupled for flexibility

# Requirements
Python 3.7+
Dependencies (included in requirements.txt):
biopython
tkinter (preinstalled with most Python versions)


