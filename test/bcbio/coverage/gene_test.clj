(ns bcbio.coverage.gene_test
  "Tests for assessing coverage of reads on genes of interest."
  (:require [clojure.java.io :as io]
            [midje.sweet :refer :all]
            [bcbio.coverage.gene :as gene]))

(background
 (around :facts
         (let [data-dir (str (io/file "test" "data"))
               bam-file (str (io/file data-dir "aligned-reads.bam"))
               gene-bed (str (io/file data-dir "gene-regions.bed"))]
           ?form)))

(facts "Identify gene regions with minimal sequencing coverage"
  (let [params {:coverage 10
                :block 10}]
    (gene/problem-coverage bam-file gene-bed params)) => nil)
