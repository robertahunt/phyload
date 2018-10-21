# Flu HA data

## [`flu_HA.fasta`](flu_HA.fasta)
This codon alignment has 49 human seasonal H1N1 sequences from 1918 to 2008.
There is at most 1 sequence per year and each sequence has 1,692 nucleotides.

## [`flu_HA.tre`](flu_HA.tre)
This is the tree inferred under GTR+GAMMA using `RAxML` version 8.2.12.
Inferred model parameters
 - alpha shape parameter = 0.440894
 - relative exchange rates (ac ag at cg ct gt) =  1.882161 7.009179 0.914813 0.495852 7.666181 1.000000
 - base frequencies = 0.340152 0.190828 0.225045 0.243974


## [`flu_HA.png`](flu_HA.png)
Here is the tree above:  
![](flu_HA.png)
Here is a fun flu fact (FFF).
If you look at the tips they march through time (pretty well) except there is a gap in the years from the 50s to the 70s.
This is because human seasonal H1N1 had gone extinct in the 50s but was re-introduced into the population in the 70s.
So 20 years had passed but because the virus had been frozen 0 _evolutionary_ years had passed.

## [`raxml`](raxml)
Inference of tree from sequence data and all associated files.
