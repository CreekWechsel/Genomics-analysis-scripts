### Requirements

- Python 3.x
- biopython
- numpy
- pandas

### 1_assembly-stat_v2.py

This script is designed to implement the same functionality as `assembly-stat`, but with more detailed output information. 

#### Usage

```shell
python assembly-stat_v2.py input.fasta
```

`input.fasta`: Any genome file

## 2_rename_fa.py

This script is designed to modify the sequence names (IDs) in a published genome's FASTA file. To use the script, you will need to prepare a table file with two columns, separated by tabs. The first column should contain the old IDs as they appear in the FASTA file, and the second column should contain the new IDs you want to assign. The script will generate a new FASTA file with the updated IDs, appending .renamed to the original file name.

#### Usage

```shell
python rename_fa.py input.fasta rename.table
```

`input.fasta`: The FASTA file whose sequence IDs you want to rename.

`rename.table`: A table file with two columns: the first column contains the old IDs, and the second column contains the new IDs.
