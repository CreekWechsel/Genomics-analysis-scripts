### Requirements

- Python 3.x
- biopython
- numpy
- pandas

### 1_assembly-stat

This script is designed to implement the same functionality as `assembly-stat`, but with more detailed output information. 

#### Usage

```shell
python assembly-stat_v2.py input.fasta
```

**Arguments**

`input.fasta`: Any genome file

**———————————————————————————————————————————————————**

### 2_rename_fa

This script is designed to modify the sequence names (IDs) in a published genome's FASTA file. To use the script, you will need to prepare a table file with two columns, separated by tabs. The first column should contain the old IDs as they appear in the FASTA file, and the second column should contain the new IDs you want to assign. The script will generate a new FASTA file with the updated IDs, appending .renamed to the original file name.

#### Usage

```shell
python rename_fa.py input.fasta rename.table
```

**Arguments**

`input.fasta`: The FASTA file whose sequence IDs you want to rename.

`rename.table`: A table file with two columns: the first column contains the old IDs, and the second column contains the new IDs.

**———————————————————————————————————————————————————**

### 3_repeat_to_gff

The main function of this script is to convert different formats of repeat sequence files (such as `.out` files generated by RepeatMasker, `.annot` files generated by RepeatProteinMask, or `.dat` files generated by TRF) into GFF3 format. Specifically:

- For `.dat` files (TRF output), it extracts sequence information and formats it into GFF3.
- For `.out` files (RepeatMasker output), it filters out simple repeat sequences and converts them into GFF3.
- For `.annot` files (RepeatProteinMask output), it extracts protein repeat information and formats it into GFF3.

#### Usage

```
python repeat_to_gff.py --in ./trf.dat --out ./trf.gff
python repeat_to_gff.py --in ./RepeatMasker.out --out ./RepeatMasker.gff
python repeat_to_gff.py --in ./Proteinmask.annot --out ./Proteinmask.gff
```

**Arguments**

`--in`:Specify the input file path. The input file must be in one of the following formats: `.dat` (TRF output), `.out` (RepeatMasker output), or `.annot` (RepeatProteinMask output).

`--out`: Specify the output GFF3 file name where the converted data will be saved.

`--pre`: Set a prefix before the repeat element ID. By default, this is an empty string. If a prefix is specified, it will be added before each repeat element ID. For example, if `--pre=MyPrefix_` is set, the ID may appear as `MyPrefix_TR100001`.

`--verbose`: Enable the output of running progress information to the screen. This is a flag, meaning that if it is specified, the program will display progress details during execution.

##### Example Usage

```bash
python script_name.py --in input.dat --out output.gff3 --pre MyPrefix_ --verbose
```
