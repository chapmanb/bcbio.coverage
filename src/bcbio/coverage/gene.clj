(ns bcbio.coverage.gene
  "Identify gene regions with low and no coverage blocks."
  (:require [clojure.java.io :as io]
            [clojure.string :as string]
            [clojure.tools.cli :refer [cli]]
            [me.raynes.fs :as fs]
            [bcbio.coverage.io.bam :as bam]
            [bcbio.coverage.io.bed :as bed]))

;; ## Coverage retrieval

(defmulti get-coverage*
  "Retrieve coverage, dispatching on type of input file."
  (fn [ftype & args] ftype))

(defmethod get-coverage* :bed
  ^{:doc "Retrieve coverage at a position from a BED file of by-position coverage"}
  [_ source contig pos])

(defmethod get-coverage* :bam
  ^{:doc "Retrieve position coverage directly from a BAM file"}
  [_ source contig pos]
  (with-open [iter (.queryOverlapping source contig (inc pos) (inc pos))]
    (count (iterator-seq iter))))

(defprotocol RetrieveCoverage
  (get-coverage [this contig pos]))

(defrecord CoverageRetriever [ftype source]
  RetrieveCoverage
  (get-coverage [_ contig pos]
    (get-coverage* ftype source contig pos))
  java.io.Closeable
  (close [_]
    (.close source)))

(defn get-coverage-retriever
  "Get a retriever that calculates coverage per position, handling multiple input types."
  [in-file]
  (let [ftype (-> in-file fs/extension string/lower-case (subs 1) keyword)
        source (case ftype
                 :bam (bam/get-bam-source in-file))]
    (CoverageRetriever. ftype source)))

(defn region-problem-coverage
  "Calculate stats for problematic coverage in a chromosome region."
  [retriever contig start end params]
  (let [cov-by-pos (reduce (fn [coll i]
                             (let [n (get-coverage retriever contig i)]
                               (if (< n (:coverage params))
                                 (assoc coll i n)
                                 coll)))
                           {} (range start end))]
    {:count (count cov-by-pos)
     :size (- end start)}))

(defn gene-problem-coverage
  "Calculate stats for problematic coverage across a set of gene regions."
  [retriever coords params]
  (let [regions (map #(region-problem-coverage retriever (:chr %) (:start %) (:end %) params)
                     coords)]
    {:name (-> coords first :name)
     :percent-nocoverage (* 100.0 (/ (apply + (map :count regions))
                                     (apply + (map :size regions))))}))

(defn problem-coverage
  "Identify problematic coverage for all supplied gene regions."
  [coverage-file gene-file params]
  (let [retriever (get-coverage-retriever coverage-file)]
    (with-open [rdr (io/reader gene-file)]
      (doseq [coords (map second (group-by :name (bed/get-iterator rdr)))]
        (println (gene-problem-coverage retriever coords params))))))

(defn -main [& args]
  (let [[options args banner] (cli args)]
    (when (not= (count args) 2)
      (println "Usage: gene <BAM or BED coverage file> <BED file of gene regions>")
      (System/exit 1))
    (apply problem-coverage args)))
