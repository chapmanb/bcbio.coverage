(ns bcbio.coverage.io.bigwig
  "Parse and query BigWig files: http://genome.ucsc.edu/goldenPath/help/bigWig.html"
  (:import [org.broad.igv.bbfile BBFileReader]))

(defn get-source
  [bw-file]
  (proxy [BBFileReader java.io.Closeable] [bw-file]
    (close [])))
