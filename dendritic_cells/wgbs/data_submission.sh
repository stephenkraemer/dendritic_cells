# # Transfer raw data via aspera

# ## Private ssh key

# private ssh key was not found / maybe could not even be copied to home on odcf-transfer
# copy it to cluster filesystem instead

mkdir /icgc/dkfzlsdf/analysis/hs_ontogeny/notebook-data/NlIosmWpSIbw7MA8/OV3eavXF6PkzOpYz/ssh
scp sra-2.ssh.priv kraemers@odcf-transfer.dkfz.de:/icgc/dkfzlsdf/analysis/hs_ontogeny/notebook-data/NlIosmWpSIbw7MA8/OV3eavXF6PkzOpYz/ssh


# ## Upload data

# - failed with "disk write failed", a known problem. Solved by submitting all files individually in for loop

#cd /icgc/dkfzlsdf/analysis/hs_ontogeny/notebook-data/NlIosmWpSIbw7MA8/OV3eavXF6PkzOpYz/dc-methylomes
cd /icgc/dkfzlsdf/analysis/hs_ontogeny/notebook-data/NlIosmWpSIbw7MA8/kgBLAl9d1HSGDSCI/dc-methylomes
for fp in *.fastq.gz; do
    ascp \
        -i /icgc/dkfzlsdf/analysis/hs_ontogeny/notebook-data/NlIosmWpSIbw7MA8/OV3eavXF6PkzOpYz/ssh/sra-2.ssh.priv \
        -QT \
        -l400m \
        --symbolic-links=follow \
        -k1 \
        $fp \
        asp-sra@upload.ncbi.nlm.nih.gov:incoming
done

# ### on the command

# –Q (for adaptive flow control) – needed for disk throttling!
# –T to disable encryption
# –k1 enable resume of failed transfers
# –l (maximum bandwidth of request, try 100m ; you may need to go up or down from there for optimum transfer)
# –r recursive copy
# –i <path/key file> (note in this example, the user has placed the key in ~/.ssh/)
# --symbolic-links=follow - follow symbolic links and upload the target files, this should be the default

# Details on ncbi upload
# <directory> is either 'incoming' or 'test'.  *Please direct your uploads to the 'test' directory* until you are confident your transmission command will work as intended.
# Files deposited in the 'incoming' directory will automatically be moved into the archive. Any files in the incoming area will typically be copied into the archive within 1-4 hours and will be removed from the incoming area once the copy is complete.
# For file transfer tips please read our Aspera transfer guide:
# http://www.ncbi.nlm.nih.gov/books/NBK242625/
# - Experiment with transfers starting at 100 Mbps and working up to 400 Mbps. Select the bandwidth setting that gives good performance with unattended operation.


# GEO submission
# ==============

# Data were uploaded to
# ftp://geoftp:rebUzyi1@ftp-private.ncbi.nlm.nih.gov/uploads/stephen_kraemer_ClixiOAs/dendritic-cells-rosenbauer

# Metadata sheet name was
# metadata_dendritic-cells-rosenbauer.xlsx

#mput /icgc/dkfzlsdf/analysis/hs_ontogeny/notebook-data/NlIosmWpSIbw7MA8/OV3eavXF6PkzOpYz/dc-methylomes/mcalls_*
#mput /icgc/dkfzlsdf/analysis/hs_ontogeny/notebook-data/NlIosmWpSIbw7MA8/OV3eavXF6PkzOpYz/dc-methylomes/metadata_dendritic-cells-rosenbauer.xlsx

mput /icgc/dkfzlsdf/analysis/hs_ontogeny/notebook-data/NlIosmWpSIbw7MA8/kgBLAl9d1HSGDSCI/dc-methylomes/mcalls_*
mput /icgc/dkfzlsdf/analysis/hs_ontogeny/notebook-data/NlIosmWpSIbw7MA8/kgBLAl9d1HSGDSCI/dc-methylomes/metadata_dendritic-cells-rosenbauer.xlsx
