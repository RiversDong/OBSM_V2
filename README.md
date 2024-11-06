# OBSM_V2
Version 2 of the Optimal Branch-Specific Model (OBSM)

## Basic Information
- **Version**: 2.0
- **Authors**: Chuan Dong, Ruiyuan Li, and Chengjun Zhang
- **Email**: [zhangcj@zafu.edu.cn](mailto:zhangcj@zafu.edu.cn)

Please feel free to contact us with any comments or questions.

## Usage
The program accepts the following command-line arguments:

- `-o` : Specifies the output directory for the evolutionary force results.
- `-aln` : The file path to the codon sequence alignment.
- `-tree` : The file path to the tree representing the alignment sequences.
- `-fo` : The branch of focus for the analysis.

### Example Usage

```bash
perl obsm_v2.pl -aln /data/chuand/fusion_gene/obsm/OBSM/OBSM_program/example/dna_seq_for_paml.txt \
               -tree /data/chuand/fusion_gene/obsm/OBSM/OBSM_program/example/gene_tree.trees \
               -o /data/chuand/fusion_gene/obsm/OBSM/OBSM_program/test3 \
               -fo Ppardus

## Notes

- **Output Path**: Ensure that the output path (`-o`) does not contain any other files or directories. If there are files present in the output path, it will cause an error. It is recommended to specify an **empty directory** for the output.

- **Absolute Path**: Always use an **absolute path** for the output directory to avoid any issues with relative paths or path resolution.

- **Beta Version**: This tool is currently in **beta version** and has not been officially released. It is for internal testing and feedback purposes only.

