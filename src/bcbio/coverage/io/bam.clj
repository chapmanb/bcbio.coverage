(ns bcbio.coverage.io.bam
  "Reading from sequence reads in BAM format, leveraging Picard."
  (:import [net.sf.samtools SAMFileReader SAMFileReader$ValidationStringency]
           [net.sf.picard.sam BuildBamIndex])
  (:require [clojure.java.io :as io]
            [clojure.string :as string]
            [bcbio.run.itx :as itx]))

(defn index-bam
  "Generate BAM index, skipping if already present."
  [in-bam]
  (let [index-choices [(str in-bam ".bai") (string/replace in-bam ".bam" ".bai")]]
    (if-let [index-file (->> index-choices
                             (remove itx/needs-run?)
                             (filter #(itx/up-to-date? % in-bam))
                             first)]
      index-file
      (do
        (SAMFileReader/setDefaultValidationStringency SAMFileReader$ValidationStringency/LENIENT)
        (BuildBamIndex/createIndex (SAMFileReader. (io/file in-bam)) (io/file (first index-choices)))
        (first index-choices)))))

(defn get-bam-source
  "Retrieve Picard SAMFileReader for input BAM, with querying."
  [bam-file]
  (let [index-file (index-bam bam-file)]
    (SAMFileReader. (io/file bam-file) (io/file index-file))))
