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
    IdentitiesOnly yes
    User git
EOF

chmod 600 ~/.ssh/mito_gen_key ~/.ssh/config
chmod 700 ~/.ssh
```

`IdentitiesOnly yes` forces ssh to offer ONLY `mito_gen_key` to github.com. Without it, ssh will offer every key in your agent first and GitHub may rate-limit the connection before it ever reaches the right key — manifests as `Permission denied (publickey)` even though the key is registered.

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

## 7. Mac-side setup (one-time, parallel to §1–§4)

The Mac mirror at `~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2/` is its own git checkout (not just an rsync target). Both Mac and Triton 2 push/pull against `git@github.com:DLCrawford/MT_Genomics.git`.

### 7a. SSH key for GitHub (Mac)

Check whether you already have a key:
```bash
ls -la ~/.ssh/
cat ~/.ssh/id_ed25519.pub 2>/dev/null
```

If a key exists and you'd rather use it: copy its `.pub` line into github.com → *Settings → SSH and GPG Keys → New SSH key* (label, e.g., `MAC_mito_gen_key`). If not, generate one:
```bash
ssh-keygen -t ed25519 -C dcrawford@miami.edu -f ~/.ssh/id_ed25519
chmod 700 ~/.ssh
chmod 600 ~/.ssh/id_ed25519
cat ~/.ssh/id_ed25519.pub          # paste this into GitHub
```

Verify:
```bash
ssh -T git@github.com
# Expected: "Hi DLCrawford! You've successfully authenticated..."
```

### 7b. zsh comment behavior (only relevant on Mac)

zsh does not treat `#` as a comment in interactive shells by default — pasting commented recipes will produce `zsh: command not found: #`. Fix once per session, or permanently in `~/.zshrc`:

```bash
setopt interactive_comments
```

### 7c. Bring the Mac directory under version control

If the Mac project directory was created without `git init` (e.g., during a restructure), and the GitHub repo already has history that you want to preserve, clone the repo into a temp dir and graft its `.git/` into the Mac project. This avoids force-push and keeps history intact.

```bash
cd /tmp
git clone git@github.com:DLCrawford/MT_Genomics.git mt_clone
mv mt_clone/.git ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2/.git
rm -rf /tmp/mt_clone

cd ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2
git config user.email "dcrawford@miami.edu"
git config user.name  "DLCrawford"
git config pull.rebase false

git status                          # review: modified vs untracked vs deleted
git add -A
git commit -m 'Restructure: numbered BSUB scripts, jobs/config.sh, docs/'
git push origin main
```

After the Mac push, on Triton 2:
```bash
cd /projectnb/dcrawford/MT_Genomics2
git pull origin main
```

If Triton 2 has untracked or modified files that block the merge and Mac is authoritative:
```bash
git reset --hard HEAD                # drops tracked modifications
git clean -fd -- docs jobs           # drops blocking untracked files in those folders
git pull origin main
```

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
