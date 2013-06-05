(ns bcbio.coverage.gene_test
  "Tests for assessing coverage of reads on genes of interest."
  (:require [clojure.java.io :as io]
            [midje.sweet :refer :all]
            [bcbio.coverage.gene :as gene]))

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
