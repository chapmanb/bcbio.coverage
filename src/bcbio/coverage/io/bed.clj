(ns bcbio.coverage.io.bed
  "Parse BED format files, handling sorting and other issues."
  (:import [org.broad.tribble.bed BEDCodec]
           [org.broad.tribble.index IndexFactory]
           [org.broad.tribble AbstractFeatureReader])
  (:require [clojure.java.io :as io]
            [clojure.string :as string]))

;; ## Utilities

(defn is-bed?
  "Check if an input file appears to be BED formatted."
  [in-file]
  (letfn [(is-int? [x]
            (try
              (Integer/parseInt x)
              true
              (catch java.lang.NumberFormatException e
                false)))]
    (with-open [rdr (io/reader in-file)]
      (let [test-line (first (remove #(.startsWith % "#") (line-seq rdr)))
            parts (string/split test-line #"\t")]
        (and (>= (count parts) 3)
             (is-int? (nth parts 1))
             (is-int? (nth parts 2)))))))

(defn- merge-intervals-chr
  "Merge sorted intervals present on the same chromosome"
  [intervals]
  (loop [xs intervals
         final []]
    (if-let [x (first xs)]
      (let [prev (last final)]
        (if (and prev (< (:start x) (:end prev)))
          (recur (rest xs) (conj (vec (butlast final))
                                 (-> prev
                                     (assoc :start (min (:start x) (:start prev)))
                                     (assoc :end (max (:end x) (:end prev))))))
          (recur (rest xs) (conj final x))))
      final)))

(defn merge-intervals
  "Merge BED intervals specified as maps with :chr :start :end"
  [intervals]
  (->> intervals
       (group-by :chr)
       (map #(sort-by (juxt :start :end) %))
       (map second)
       (mapcat merge-intervals-chr)))

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
