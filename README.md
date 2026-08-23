# Relevant SP Scripts


### Swave command

```
find ~ -type f -path "*/venv-cactus-v3.2.1/bin/activate"

# activate from output above

python ~/clleva/tools/Swave/Swave.py call \
    --input_path assemblies.tsv \
    --ref_path  /home/rocm-user/clleva/real_arm/graphs/pansn/Nipponbare_chr1.fna \
    --gfa_source pggb \
    --gfa_path graphs/pggb/chr1_pangenome.combined.fa.gz.bf3285f.11fba48.8088a73.smooth.final.gfa \
    --decomposed_vcf graphs/pggb/

python  ~/clleva/tools/Swave/Swave.py  convert_seq \
    --vcf_path pggb-swave/swave.sample_level.vcf \  
    --gfa_path chr1_pangenome.combined.fa.gz.bf3285f.11fba48.ddb9d60.smooth.final.gfa \
   --ref_path  /home/rocm-user/clleva/real_arm/graphs/pansn/Nipponbare_chr1.fna 
# -> pggb-swave/swave.sample_level.converted.vcf

python ~/clleva/tools/Swave/Swave.py call \
    --input_path assemblies.tsv \
    --ref_path  /home/rocm-user/clleva/real_arm/graphs/pansn/Nipponbare_chr1.fna \
    --gfa_source minigraph \
    --gfa_path graphs/minigraph/chr1_pangenome.minigraph.gfa \
    --output_path minigraph-swave

```