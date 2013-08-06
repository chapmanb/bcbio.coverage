(ns bcbio.coverage.workflow.multicompare
  "Compare coverage across multiple experiments and technologies."
  (:require [clojure.data.csv :as csv]
            [clojure.java.io :as io]
            [clojure.tools.cli :refer [cli]]
            [clj-yaml.core :as yaml]
            [bcbio.coverage.gene :as gene]
            [bcbio.coverage.io.bed :as bed]
            [bcbio.coverage.source.ensembl :as ensembl]
            [bcbio.run.itx :as itx]))

(defn- do-compare
  "Compare multiple coverage approaches across set of regions."
  [exps region-file ref-file params]
  (let [gene-coord-file (ensembl/get-coord-bed region-file params)]
    (with-open [rdr (io/reader gene-coord-file)]
      (let [coords (vec (bed/get-iterator rdr))]
        (->> exps
             (mapcat (fn [exp]
                       (map #(assoc % :exp (:name exp))
                            (gene/coverage-report-bygene (map :coverage (:samples exp))
                                                         coords ref-file params))))
             (group-by :name)
             (sort-by first)
             (map second))))))

(defn- cmp->csv
  "Convert a comparison for a gene into CSV output"
  [cmp exps]
  (let [cmps-by-exp (group-by :exp cmp)]
    (concat [(-> cmp first :name)
             (-> cmp first :size)]
            (mapcat (fn [e]
                      [(get-in cmps-by-exp [e 0 :avg-reads])
                       (get-in cmps-by-exp [e 0 :percent-nocoverage])])
                    exps))))

(defn- write-to-csv
  [cmps exps out-file]
  (let [header (concat ["gene" "size"]
                       (mapcat (fn [x] [(str x "_avg_reads") (str x "_pct_nocov")]) exps))]
    (with-open [wtr (io/writer out-file)]
      (csv/write-csv wtr (cons header (map #(cmp->csv % exps) cmps))))))

(defn compare-from-config
  "Perform comparison, reading input information from YAML configuration file."
  [config-file out-file cores]
     (let [config (-> config-file slurp yaml/parse-string
                      (assoc-in [:params :cores] cores))
           cmps (do-compare (:experiments config) (:regions config) (:ref-file config)
                            (:params config))]
       (write-to-csv cmps (map :name (:experiments config)) out-file))
     out-file)

(defn -main [& args]
  (let [[options args banner]
        (cli args
             ["-c" "--cores" "Number of cores to use" :default 1
              :parse-fn #(Integer. %)])]
    (when (not= (count args) 2)
      (println banner)
      (println "Usage: multicompare <config file> <output CSV file>")
      (System/exit 1))
    (compare-from-config (first args) (second args) (:cores options))))
