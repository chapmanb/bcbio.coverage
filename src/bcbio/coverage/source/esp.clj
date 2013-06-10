(ns bcbio.coverage.source.esp
  "Coverage from the Exome Sequencing Project (http://evs.gs.washington.edu/EVS/).
   Covert custom files into UCSC compatible BigWig for query."
  (:require [clojure.string :as string]
            [clojure.java.io :as io]
            [clojure.java.shell :as shell]
            [me.raynes.fs :as fs]
            [bcbio.align.gref :as gref]
            [bcbio.run.itx :as itx]
            [bcbio.variation.variantcontext :as gvc]))

;; ## ESP coverage

(defn- wig->bigwig
  "Convert wig file into binary BigWig with UCSC's wigToBigWig"
  [wig-file ref-file]
  (let [bigwig-file (str (itx/file-root wig-file) ".bw")
        seq-dict (gref/get-seq-dict ref-file)]
    (when (itx/needs-run? bigwig-file)
      (itx/with-tx-file [tx-bigwig-file bigwig-file]
        (itx/with-temp-dir [tmp-dir (fs/parent wig-file)]
          (let [sizes-file (str (io/file tmp-dir "chrom-sizes.txt"))]
            (with-open [wtr (io/writer sizes-file)]
              (doseq [cur-chrom (.getSequences seq-dict)]
                (.write wtr (format "%s\t%s\n" (.getSequenceName cur-chrom)
                                    (.getSequenceLength cur-chrom)))))
            (shell/sh "wigToBigWig" wig-file sizes-file tx-bigwig-file)))))
    bigwig-file))

(defn- coverage->wig
  "Covert directory of ESP coverage into wiggle file"
  [esp-dir]
  (let [coverage-files (fs/glob (fs/file esp-dir) "*.coverage.all_sites.txt")
        wig-file (str (io/file esp-dir
                               (-> coverage-files first fs/base-name (string/split #"\.") first
                                   (str "-coverage.wig"))))]
    (when (itx/needs-run? wig-file)
      (itx/with-tx-file [tx-wig-file wig-file]
        (with-open [wtr (io/writer tx-wig-file)]
          (doseq [cov-file coverage-files]
            (let [chrom (-> cov-file str (string/split #"\.") second (string/replace "chr" ""))]
              (with-open [rdr (io/reader cov-file)]
                (.write wtr (format "variableStep chrom=%s\n" chrom))
                (doseq [line (remove #(.startsWith % "#") (line-seq rdr))]
                  (let [parts (string/split line #" ")]
                    (.write wtr (format "%s %s\n" (nth parts 1) (nth parts 3)))))))))))
    wig-file))

(defn coverage->bigwig
  "Covert directory of ESP coverage into BigWig format.
   Works from an unpacked directory from ESP coverage downloads:
   http://evs.gs.washington.edu/evs_bulk_data/ESP6500SI-V2.coverage.all_sites.txt.tar.gz"
  [esp-dir ref-file]
  (wig->bigwig (coverage->wig esp-dir) ref-file))

;; ## ESP variants

(defn- above-maf?
  "Check if a variant is above a minimum threshold allele frequency.
   In ESP VCF files, the global (All) frequency is the last attribute
   in the MAF INFO field."
  [min-maf vc]
  (>= (Float/parseFloat (last (get-in vc [:attributes "MAF"])))
      min-maf))

(defn variants-in-region
  "Retrieve variants in a genomic region with a minimum allele frequency."
  [vcf-file ref-file region min-maf]
  (with-open [vcf-get (gvc/get-vcf-retriever ref-file vcf-file)]
    (vec (filter (partial above-maf? min-maf)
                 (gvc/variants-in-region vcf-get region)))))
