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
