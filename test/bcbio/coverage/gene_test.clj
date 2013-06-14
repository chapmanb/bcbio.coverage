(ns bcbio.coverage.gene_test
  "Tests for assessing coverage of reads on genes of interest."
  (:require [clojure.java.io :as io]
            [midje.sweet :refer :all]
            [bcbio.coverage.gene :as gene]
            [bcbio.coverage.io.bam :as bam]
            [bcbio.coverage.io.bed :as bed]
            [bcbio.coverage.source.esp :as esp]
            [bcbio.coverage.workflow.wgsexome :as wgsexome]
            [bcbio.run.itx :as itx]))

(background
 (around :facts
         (let [data-dir (str (io/file "test" "data"))
               bam-file (str (io/file data-dir "aligned-reads.bam"))
               ref-file (str (io/file data-dir "GRCh37.fa"))
               bw-file (str (io/file data-dir "esp-coverage.bw"))
               gene-bed (str (io/file data-dir "gene-regions.bed"))]
           ?form)))

(facts "Identify gene regions with minimal sequencing coverage"
  (let [params {:coverage 10
                :downsample 50
                :block {:min 20 :distance 5}}
        bases [1 3 4 5 8 10 11 14]]
    (gene/split-into-blocks 3 bases) => [[1 3 4 5 8 10 11 14]]
    (gene/split-into-blocks 2 bases) => [[1 3 4 5] [8 10 11] [14]]
    (gene/split-into-blocks 1 bases) => [[1] [3 4 5] [8] [10 11] [14]]
    (gene/problem-coverage bam-file gene-bed ref-file params) => nil
    (gene/problem-coverage bw-file gene-bed ref-file params) => nil))

(facts "Convert gene names into coordinates"
  (let [name-file (str (io/file data-dir "genenames.txt"))
        name-exon-file (str (io/file data-dir "genenames-exons.bed"))]
    (itx/remove-path name-exon-file)
    (gene/get-coord-bed gene-bed {}) => gene-bed
    ;; Test calls out to Ensembl, slow
    ;;(gene/get-coord-bed name-file {:organism "human"}) => name-exon-file
    ))

(facts "Manipulations for intervals and BED files"
  (bed/merge-intervals []) => []
  (bed/merge-intervals [{:chr "A" :start 1 :end 3}
                        {:chr "A" :start 2 :end 5}
                        {:chr "A" :start 7 :end 9}
                        {:chr "B" :start 1 :end 3}]) =>
                        [{:chr "A" :start 1 :end 5}
                         {:chr "A" :start 7 :end 9}
                         {:chr "B" :start 1 :end 3}])

(facts "Extract variants from ESP by region and allele frequency"
  (let [esp-vcf-file (str (io/file data-dir "esp-snps_indels.vcf"))
        vcs (esp/variants-in-region esp-vcf-file ref-file
                                    {:chr "22" :start 6900 :end 7300} 1.0)]
    (map (juxt :chr :start) vcs) => [["22" 6920] ["22" 7226]]))

(facts "Compare whole genome and exome coverage"
  (let [wgs-config (str (io/file data-dir "wgsexome-compare.yaml"))
        out-file (str (itx/file-root wgs-config) ".csv")]
    (itx/remove-path out-file)
    (wgsexome/compare-from-config wgs-config out-file) => out-file))

(facts "GATK based iterator over piled up read regions."
  (with-open [iter (bam/prep-bam-region-iter [bam-file] ref-file {:chr "MT" :start 250 :end 300})]
    (-> (bam/get-align-contexts iter) first .size) => 5))
