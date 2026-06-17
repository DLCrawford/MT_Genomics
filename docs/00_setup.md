# 00 — Setup (Triton 2 + GitHub + conda)

This is the one-time setup. Once everything below is in place, every future session only needs the [Per-session](#per-session-startup) block at the bottom.

## Hosts

| Host | SSH | Notes |
|---|---|---|
| Triton 2 (current HPC) | `ssh dcrawford@t2.idsc.miami.edu` | LSF / BSUB scheduler. Working dir: `/projectnb/dcrawford/MT_Genomics2/`. |
| Triton (legacy) | `ssh dcrawford@triton.ccs.miami.edu` | Retired. Original BSUB scripts ran here. |
| GitHub | `git@github.com:DLCrawford/MT_Genomics.git` | Private repo. SSH auth via `mito_gen_key`. |

## 1. SSH key for GitHub (one-time, on Triton 2)

```bash
ssh-keygen -t ed25519 -C dcrawford@miami.edu -f ~/.ssh/mito_gen_key

cat > ~/.ssh/config << 'EOF'
Host github.com
    IdentityFile ~/.ssh/mito_gen_key
    User git
EOF

chmod 600 ~/.ssh/mito_gen_key ~/.ssh/config
chmod 700 ~/.ssh
```
# The key fingerprint is:
	SHA256:HAsiLgnYPa0ZF4JyD9QclFCgbtofoqWhWxGXPimAvoo dcrawford@miami.edu
Add `~/.ssh/mito_gen_key.pub` to GitHub → *Settings → SSH and GPG Keys → New SSH key* (title: `mito_gen_key`).

Verify:
```bash
ssh -T git@github.com
# Expected: "Hi DLCrawford! You've successfully authenticated..."
```

## 2. Initialize the working directory

```bash
cd /projectnb/dcrawford
mkdir -p MT_Genomics2 && cd MT_Genomics2
git init
git config user.email "dcrawford@miami.edu"
git config user.name  "DLCrawford"
git remote add origin git@github.com:DLCrawford/MT_Genomics.git
git config pull.rebase false
git pull origin main
```

## 3. Conda env

Triton's default `git` is too old, so we use the conda env for git too. Same env carries the bioinformatics tools, at least until conflicts force a split.

```bash
module load anaconda3
conda create -n mito_genomics -c bioconda -c conda-forge -y \
    git fastqc multiqc trimmomatic bwa samtools bcftools
conda activate mito_genomics
```

If a tool conflicts during install, fall back to a stage-specific env (e.g. `fastqc_env`, `BioInfo_env`) and update the relevant `jobs/*.sh` to activate it. Document the split here.

## 4. ~/.bashrc — make conda available in every interactive shell

```bash
# Add to ~/.bashrc
if [[ $- == *i* ]]; then
    module load anaconda3
    source "$(conda info --base)/etc/profile.d/conda.sh"
fi
```

The interactive-shell guard (`[[ $- == *i* ]]`) keeps `module load` from running in non-interactive contexts (job scripts, scp, rsync) where it can cause trouble.

## 5. Reference genome (one-time)

The MT reference + GFF live at `/projectnb/dcrawford/SSM_Mito/Fh_MT_ref/`.

```bash
conda activate mito_genomics

# (optional) install NCBI EDirect if you need to pull NC_012312.1 fresh
sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
esearch -db nucleotide -query "NC_012312.1" \
  | efetch -format fasta > Fhet_MT.fasta

# MT annotations (subset the full Fhet GFF to chrMT)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/125/445/GCF_011125445.2_MU-UCD_Fhet_4.1/GCF_011125445.2_MU-UCD_Fhet_4.1_genomic.gff.gz
gunzip GCF_011125445.2_MU-UCD_Fhet_4.1_genomic.gff.gz
grep "NC_012312.1" GCF_011125445.2_MU-UCD_Fhet_4.1_genomic.gff > Fhet_MT.gff

# Index the reference for bwa + samtools
REF=/projectnb/dcrawford/SSM_Mito/Fh_MT_ref/Fhet_MT.fasta
samtools faidx "$REF"
bwa index "$REF"
```

## 6. Logs directory

```bash
mkdir -p /projectnb/dcrawford/MT_Genomics2/logs
```

All `jobs/*.sh` write `.out` / `.err` here.

## Per-session startup

Every new shell on Triton 2:

```bash
module load anaconda3
conda activate mito_genomics
cd /projectnb/dcrawford/MT_Genomics2
git pull origin main
```

## Standard git workflow

```bash
git add .
git commit -m "brief description of what changed"
git push origin main
```

## Sync from Mac (push local changes to Triton 2)

From Mac terminal:

```bash
rsync -avzhm --progress \
  ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2/ \
  dcrawford@t2.idsc.miami.edu:/projectnb/dcrawford/MT_Genomics2/
```

Reverse direction (HPC → Mac, e.g. fresh `.md` files):

```bash
rsync -avzhm --progress \
  --include='*/' --include='*.md' --include='CHANGE*' --exclude='*' \
  dcrawford@t2.idsc.miami.edu:/projectnb/dcrawford/MT_Genomics2/ \
  ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2/
```
