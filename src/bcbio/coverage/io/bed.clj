(ns bcbio.coverage.io.bed
  "Parse BED format files, handling sorting and other issues."
  (:import [org.broad.tribble.bed BEDCodec]
           [org.broad.tribble.index IndexFactory]
           [org.broad.tribble AbstractFeatureReader])
  (:require [clojure.java.io :as io]
            [clojure.string :as string]))

;; ## Tribble parsing

;; XXX Handle full sorting and parsing with Tribble
(declare sort-bed-file)

(defn get-bed-source
  "Provide tribble feature source for a BED formatted file."
  [bed-file ref-file]
  (let [batch-size 500
        work-bed (sort-bed-file bed-file ref-file)
        idx (IndexFactory/createIntervalIndex (io/file work-bed) (BEDCodec.) batch-size)]
    (AbstractFeatureReader/getFeatureReader work-bed (BEDCodec.) idx)))

;; ## Quick and dirty parsing

(defn- parse-bed-line
  [line]
  (let [[chrom start end name] (take 4 (string/split line #"\t"))]
    {:chr (string/trim chrom)
     :start (Integer/parseInt start)
     :end (Integer/parseInt end)
     :name (string/trim name)}))

(defn get-iterator
  ([bed-file ref-file]
     (.iterator (get-bed-source bed-file ref-file)))
  ([rdr]
     (map parse-bed-line (line-seq rdr))))
