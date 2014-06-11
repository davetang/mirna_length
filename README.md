README
-----------------------------

### Generate the transition matrix and first base frequency

`R --no-save < dinucleotide.R`

### Count bases in the fasta file

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

### Save these genome frequencies as Robjects

`R --no-save < genome_freq.R`

### Generate random sequence using genome nucleotide frequencies

`seq_by_gen_freq.R 22 human_nuc_freq.Robject human`

### Generate random sequence using dinucleotide frequencies

`seq_by_markov_chain.R 22 human_trans_mat.Robject human_init_prob.Robject human`

### Generate random sequence using equal nucleotide frequencies

`seq_by_equal.R 22 human`
