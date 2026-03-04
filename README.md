# Peptide_2_Mol
Repo for extracting pharmacophores from peptides and finding molecules form there
This is the creation of the repo with this thread on gemini: https://gemini.google.com/share/0ed57d383f41


Steps:

python /mnt/d/GitHub_repos/Peptide_2_Mol/boltz_alanine_scan.py -i /mnt/e/TFEB_hopping/F8.pdb -o /mnt/e/boltz_inputs -b Rags -c A -r "1-24"

for i in boltz_inputs/*.yaml; do boltz predict $i --use_msa_server --out_dir boltz_outputs/; done

python /mnt/d/GitHub_repos/Peptide_2_Mol/parse_boltz_results.py -d ./boltz_outputs -w "_WT" -o hotspot_analysis.csv
