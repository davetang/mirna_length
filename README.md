README
-----------------------------

Run pipeline.sh with the number of random sequences to generate as the first parameter. It is possible to run multiple instances of pipeline.sh, to speed things up (though this hasn't been extensively tested). The script downloads necessary files, generates the random sequences, and maps them back to the respective genomes.

This work was inspired by an observation I made back in 2012 (http://davetang.org/muse/2012/03/08/why-mirna-are-22-or-23-nucleotides-long/).

### Clone this repository

`git clone https://github.com/davetang/mirna_length.git`

### Perform the entire analysis

`pipeline.sh 1000000`

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

`seq_by_gen_freq.R 22 hg38_nuc_freq.Robject 1000000 hg38`

### Generate random sequence using dinucleotide frequencies

`seq_by_markov_chain.R 22 hg38_trans_mat.Robject hg38_init_prob.Robject 1000000 hg38`

### Generate random sequence using equal nucleotide frequencies

`seq_by_equal.R 22 1000000 hg38`

### Generate the transition matrices

`Rscript image/transition.R ce10_trans_mat.Robject`<br />
`Rscript image/transition.R danRer7_trans_mat.Robject`<br />
`Rscript image/transition.R hg38_trans_mat.Robject`<br />
`Rscript image/transition.R mm10_trans_mat.Robject`<br />
`epstopdf image/ce10.eps`<br />
`epstopdf image/danRer7.eps`<br />
`epstopdf image/hg38.eps`<br />
`epstopdf image/mm10.eps`
