"""Default build specific values (i.e. URLs) for each genome build"""

build_37_urls = {
    "reference": {
        "id": "hg19",
        "name": "Human (GRCh37/hg19)",
        "fastaURL": "https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg19/hg19.fasta",
        "indexURL": "https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg19/hg19.fasta.fai",
        "cytobandURL": "https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg19/cytoBand.txt",
        "aliasURL": (
            "https://s3.amazonaws.com/igv.org.genomes/hg19/hg19_alias.tab"
        ),
    },
    "tracks": [
        {
            "name": "Refseq Genes",
            "format": "refgene",
            "id": "hg19_genes",
            "url": "https://s3.amazonaws.com/igv.org.genomes/hg19/ncbiRefSeq.sorted.txt.gz",
            "indexURL": "https://s3.amazonaws.com/igv.org.genomes/hg19/ncbiRefSeq.sorted.txt.gz.tbi",
            "visibilityWindow": -1,
            "supportsWholeGenome": False,
            "removable": False,
            "order": 1000000,
            "infoURL": "https://www.ncbi.nlm.nih.gov/gene/?term=$$",
            "type": "annotation",
        }
    ],
}

build_38_urls = {
    "reference": {
        "id": "hg38",
        "name": "Human (GRCh38/hg38)",
        "fastaURL": "https://igv.org/genomes/data/hg38/hg38.fa",
        "indexURL": "https://igv.org/genomes/data/hg38/hg38.fa.fai",
        "cytobandURL": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBandIdeo.txt.gz",
        "aliasURL": "https://igv.org/genomes/data/hg38/hg38_alias.tab",
        "twoBitURL": (
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit"
        ),
        "chromSizesURL": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes",
    },
    "tracks": [
        {
            "name": "Refseq All",
            "format": "refgene",
            "url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz",
            "order": 0,
            "type": "annotation",
            "height": 150,
            "color": "rgb(255, 41, 135)",
        },
        {
            "name": "MANE Transcripts",
            "format": "bigbed",
            "url": "https://hgdownload.soe.ucsc.edu/gbdb/hg38/mane/mane.bb",
            "visibilityWindow": -1,
            "order": 0,
            "type": "annotation",
            "height": 150,
        },
        {
            "name": "Refseq Select",
            "format": "refgene",
            "url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeqSelect.txt.gz",
            "indexed": False,
            "order": 1000001,
            "infoURL": "https://www.ncbi.nlm.nih.gov/gene/?term=$$",
            "type": "annotation",
            "height": 150,
        },
    ],
}
