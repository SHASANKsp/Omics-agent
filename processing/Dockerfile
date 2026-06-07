# ─────────────────────────────────────────────────────────────
#  Omics Processing Image
#  Tools: SRA Toolkit · fastp · STAR (+ STARsolo) · HISAT2
#         featureCounts (subread) · MACS3 · samtools · Bismark
#  Base:  Ubuntu 22.04 LTS
#
#  URL notes (learned the hard way):
#   - HISAT2 binaries: cloud.biohpc.swmed.edu, not GitHub
#   - Subread binaries: SourceForge, not GitHub (GitHub only has v2.0.2)
#   - fastp binary:    opengene.org/fastp/fastp (always latest)
#   - SRA Toolkit:     sratoolkit.current-ubuntu64.tar.gz (always latest)
#   - STAR, samtools, Bismark: standard GitHub releases (work fine)
# ─────────────────────────────────────────────────────────────
FROM ubuntu:22.04

LABEL maintainer="omics-pipeline"
LABEL description="Omics processing: STAR/STARsolo, HISAT2, featureCounts, MACS3, Bismark, SRA Toolkit"

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC

# ── System packages ───────────────────────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget curl git ca-certificates \
    build-essential cmake \
    python3 python3-pip python3-dev \
    zlib1g-dev libbz2-dev liblzma-dev \
    libncurses5-dev libcurl4-openssl-dev \
    libssl-dev libxml2-dev \
    perl \
    pigz \
    unzip \
    && rm -rf /var/lib/apt/lists/*

# ── samtools 1.21 ─────────────────────────────────────────────
# Standard GitHub release — URL format is stable
RUN wget -q https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 \
    && tar -xjf samtools-1.21.tar.bz2 \
    && cd samtools-1.21 \
    && ./configure --prefix=/usr/local \
    && make -j8 \
    && make install \
    && cd .. && rm -rf samtools-1.21 samtools-1.21.tar.bz2 \
    && samtools --version | head -1

# ── STAR 2.7.11b ─────────────────────────────────────────────
# Static binary from GitHub releases — works reliably
RUN wget -q https://github.com/alexdobin/STAR/releases/download/2.7.11b/STAR_2.7.11b.zip \
    && unzip -q STAR_2.7.11b.zip \
    && mv STAR_2.7.11b/Linux_x86_64_static/STAR /usr/local/bin/STAR \
    && rm -rf STAR_2.7.11b STAR_2.7.11b.zip \
    && STAR --version

# ── HISAT2 2.2.1 ─────────────────────────────────────────────
# IMPORTANT: binary is on cloud.biohpc.swmed.edu, NOT GitHub releases.
# GitHub release page only has source + a link pointing here.
RUN wget -q \
        "https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download" \
        -O hisat2.zip \
    && unzip -q hisat2.zip \
    && HDIR=$(find . -maxdepth 1 -name "hisat2-*" -type d | head -1) \
    && cp ${HDIR}/hisat2         /usr/local/bin/ \
    && cp ${HDIR}/hisat2-align-s /usr/local/bin/ \
    && cp ${HDIR}/hisat2-align-l /usr/local/bin/ \
    && cp ${HDIR}/hisat2-build   /usr/local/bin/ \
    && cp ${HDIR}/hisat2-build-s /usr/local/bin/ \
    && cp ${HDIR}/hisat2-build-l /usr/local/bin/ \
    && cp ${HDIR}/hisat2-inspect /usr/local/bin/ \
    && cp ${HDIR}/*.py           /usr/local/bin/ \
    && rm -rf "${HDIR}" hisat2.zip \
    && hisat2 --version | head -1

# ── Subread 2.1.0 (featureCounts) ────────────────────────────
# IMPORTANT: GitHub only has v2.0.2. v2.0.3 onwards are SourceForge only.
# Latest stable: 2.1.0 (March 2025). Use -O to name the file explicitly
# because SourceForge download URLs have no filename in the path.
RUN wget -q \
        "https://sourceforge.net/projects/subread/files/subread-2.1.0/subread-2.1.0-Linux-x86_64.tar.gz/download" \
        -O subread.tar.gz \
    && tar -xzf subread.tar.gz \
    && SDIR=$(find . -maxdepth 1 -name "subread-*-Linux-x86_64" -type d | head -1) \
    && cp -r ${SDIR}/bin/* /usr/local/bin/ \
    && rm -rf "${SDIR}" subread.tar.gz \
    && featureCounts -v 2>&1 | head -1

# ── fastp (latest stable) ─────────────────────────────────────
# opengene.org/fastp/fastp always serves the latest pre-built Linux binary.
# More reliable for Docker than pinned GitHub release URLs.
RUN wget -q http://opengene.org/fastp/fastp \
        -O /usr/local/bin/fastp \
    && chmod +x /usr/local/bin/fastp \
    && fastp --version 2>&1 | head -1

# ── SRA Toolkit (current stable) ─────────────────────────────
# NCBI recommends sratoolkit.current-* URL — avoids broken pinned version paths.
# Resolves to 3.4.1 as of June 2026. Directory name varies by version so we
# use find to locate it rather than hardcoding.
RUN wget -q \
        "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz" \
        -O sratoolkit.tar.gz \
    && tar -xzf sratoolkit.tar.gz \
    && SRADIR=$(find . -maxdepth 1 -name "sratoolkit.*-ubuntu64" -type d | head -1) \
    && cp ${SRADIR}/bin/prefetch      /usr/local/bin/ \
    && cp ${SRADIR}/bin/fasterq-dump  /usr/local/bin/ \
    && cp ${SRADIR}/bin/vdb-config    /usr/local/bin/ \
    && rm -rf sratoolkit.tar.gz "${SRADIR}" \
    && prefetch --version 2>&1 | head -1

# ── MACS3 (peak caller for ATAC-seq / ChIP-seq) ──────────────
# pip install is the canonical method — no binary hosting issues.
RUN pip3 install --no-cache-dir macs3 \
    && macs3 --version

# ── Bowtie2 2.5.4 ────────────────────────────────────────────
# Required by Bismark for genome indexing and alignment.
# Must be installed BEFORE Bismark.
RUN wget -q \
        https://github.com/BenLangmead/bowtie2/releases/download/v2.5.4/bowtie2-2.5.4-linux-x86_64.zip \
        -O bowtie2.zip \
    && unzip -q bowtie2.zip \
    && cp bowtie2-2.5.4-linux-x86_64/bowtie2*  /usr/local/bin/ \
    && rm -rf bowtie2-2.5.4-linux-x86_64 bowtie2.zip \
    && bowtie2 --version | head -1

# ── Bismark v0.25.1 (methylation) ────────────────────────────
# Tags use 'v' prefix from v0.24.2 onwards (e.g. v0.25.1, not 0.25.1).
# The extracted directory name drops the 'v': Bismark-0.25.1
RUN wget -q \
        https://github.com/FelixKrueger/Bismark/archive/refs/tags/v0.25.1.tar.gz \
        -O bismark.tar.gz \
    && tar -xzf bismark.tar.gz \
    && mv Bismark-0.25.1 /opt/bismark \
    && ln -s /opt/bismark/bismark                       /usr/local/bin/bismark \
    && ln -s /opt/bismark/bismark_genome_preparation    /usr/local/bin/bismark_genome_preparation \
    && ln -s /opt/bismark/bismark2bedGraph              /usr/local/bin/bismark2bedGraph \
    && ln -s /opt/bismark/coverage2cytosine             /usr/local/bin/coverage2cytosine \
    && ln -s /opt/bismark/bismark_methylation_extractor /usr/local/bin/bismark_methylation_extractor \
    && rm -f bismark.tar.gz \
    && bismark --version 2>&1 | head -1

# ── Python packages ───────────────────────────────────────────
RUN pip3 install --no-cache-dir \
    requests \
    rich \
    pyyaml

# ── SRA Toolkit: disable interactive config prompts ──────────
RUN mkdir -p /root/.ncbi \
    && printf '/LIBS/GUID = "00000000-0000-0000-0000-000000000000"\n/config/default = "false"\n' \
       > /root/.ncbi/user-settings.mkfg

# ── Working directories ───────────────────────────────────────
RUN mkdir -p /data /references /output /tmp/sra_cache
ENV SRA_CACHE=/tmp/sra_cache

WORKDIR /data
CMD ["/bin/bash"]
