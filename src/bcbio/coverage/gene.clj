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

(defn split-into-blocks
  "Split a sequence of positions into blocks of consecutive items within n bases."
  [n in-xs]
  (loop [xs (sort in-xs)
         cur-block []
         final []]
    (cond
     (empty? xs)
     (conj final cur-block)
     (or (empty? cur-block)
         (<= (- (first xs) (last cur-block)) n))
     (recur (rest xs) (conj cur-block (first xs)) final)
     :else
     (recur (rest xs) [(first xs)] (conj final cur-block)))))

(defn- identify-nocoverage-blocks
  "Extract blocks of nocoverage bases from raw positions"
  [cov-by-pos params]
  (->> (keys cov-by-pos)
       (split-into-blocks (get-in params [:block :distance]))
       (remove empty?)
       (map (juxt first last))
       (filter (fn [[s e]] (> (- e s) (get-in params [:block :min]))))))

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
     :blocks (identify-nocoverage-blocks cov-by-pos params)
     :size (- end start)}))

(defn gene-problem-coverage
  "Calculate stats for problematic coverage across a set of gene regions."
  [retriever coords params]
  (let [regions (map #(region-problem-coverage retriever (:chr %) (:start %) (:end %) params)
                     coords)]
    {:name (-> coords first :name)
     :coords (map #(select-keys % [:chr :start :end]) coords)
     :blocks (mapcat :blocks regions)
     :size (apply + (map :size regions))
     :percent-nocoverage (* 100.0 (/ (apply + (map :count regions))
                                     (apply + (map :size regions))))}))

(defn problem-coverage
  "Identify problematic coverage for all supplied gene regions.
   params is a map of attributes to influence what we identify
    :coverage - Minimum coverage required for a position to be considered covered.
    :block - Map of parameters for calculating blocks of nocoverage
      :min - Minimum size of a block to report
      :distance - Allowed distance between nocoverage bases to be considered in a block"
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
    (apply problem-coverage (concat args [{:coverage 10
                                           :block {:min 50 :distance 5}}]))))
