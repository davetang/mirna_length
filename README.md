README
-----------------------------

### Generate the transition matrices and first base frequencies

`R --no-save < dinucleotide.R`

### Count nucleotide frequencies in the fasta files

`count_base.pl hg38.fa`

A: 898285419
C: 623727342
G: 626335137
T: 900967885
Other: 159970322

`count_base.pl ce10.fa`

A: 32371753
C: 17782012
G: 17759040
T: 32373265
Other: 0

`count_base.pl danRer7.fa`

A: 446195643
C: 258711649
G: 258833640
T: 446029177
Other: 2694734

`count_base.pl mm10.fa`

A: 773280124
C: 552647817
G: 552690118
T: 774165441
Other: 78088274

### Save nucleotide frequencies from different genomes as Robjects

`R --no-save < genome_freq.R`

### Generate random sequence using genome nucleotide frequencies

`seq_by_gen_freq.R 22 human_nuc_freq.Robject 1000000 human`

### Generate random sequence using dinucleotide frequencies

`seq_by_markov_chain.R 22 human_trans_mat.Robject human_init_prob.Robject 1000000 human`

### Generate random sequence using equal nucleotide frequencies

`seq_by_equal.R 22 1000000 human`

### Generate the transition matrices

`image/transition.R human_trans_mat.Robject`
