# Pre-run script pasted into the Seqera Cloud pipeline settings ("Pre-run script").
# 1) Bactopia's sample-table reader needs a LOCAL file, so copy the table from cloud storage.
# 2) Fill the dataset cache once: the in-pipeline download was too slow on Fusion (hit the 8 h limit).
# Replace <your-bucket> with the bucket of your compute environment.
aws s3 cp s3://<your-bucket>/inputs/all120/samples_all120.txt /tmp/samples_all120.txt
set -o pipefail
D=s3://<your-bucket>/bactopia-datasets
B=https://datasets.bactopia.com/datasets
for f in v4.1.0/amrfinderplus.tar.gz mlst.tar.gz gtdb-rs207.genomic-reps.dna.k31.lca.json.gz mash-refseq88.k21.msh.xz; do
  if ! aws s3 ls "$D/$f" > /dev/null 2>&1; then
    echo "Fetching $f"
    if ! curl -fsSL --retry 5 "$B/$f" | aws s3 cp - "$D/$f"; then
      echo "Download of $f failed, removing partial file"
      aws s3 rm "$D/$f" > /dev/null 2>&1
      exit 1
    fi
  fi
done
