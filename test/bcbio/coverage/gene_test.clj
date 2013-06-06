(ns bcbio.coverage.gene_test
  "Tests for assessing coverage of reads on genes of interest."
  (:require [clojure.java.io :as io]
            [midje.sweet :refer :all]
            [bcbio.coverage.gene :as gene]
            [bcbio.coverage.io.bed :as bed]
            [bcbio.run.itx :as itx]))

(background
 (around :facts
         (let [data-dir (str (io/file "test" "data"))
               bam-file (str (io/file data-dir "aligned-reads.bam"))
               bw-file (str (io/file data-dir "esp-coverage.bw"))
               gene-bed (str (io/file data-dir "gene-regions.bed"))]
           ?form)))

(facts "Identify gene regions with minimal sequencing coverage"
  (let [params {:coverage 10
                :block {:min 20 :distance 5}}
        bases [1 3 4 5 8 10 11 14]]
    (gene/split-into-blocks 3 bases) => [[1 3 4 5 8 10 11 14]]
    (gene/split-into-blocks 2 bases) => [[1 3 4 5] [8 10 11] [14]]
    (gene/split-into-blocks 1 bases) => [[1] [3 4 5] [8] [10 11] [14]]
    (gene/problem-coverage bam-file gene-bed params) => nil
    (gene/problem-coverage bw-file gene-bed params) => nil))

(facts "Convert gene names into coordinates"
  (let [name-file (str (io/file data-dir "genenames.txt"))
        name-exon-file (str (io/file data-dir "genenames-exons.bed"))]
    (itx/remove-path name-exon-file)
    (gene/get-coord-bed gene-bed {}) => gene-bed
    ;; Test calls out to Ensembl, slow
    ;; (gene/get-coord-bed name-file {:organism "human"}) => name-exon-file
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
