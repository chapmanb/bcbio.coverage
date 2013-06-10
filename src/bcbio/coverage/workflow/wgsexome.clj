(ns bcbio.coverage.workflow.wgsexome
  "Compare whole genome sequencing data to exome data, assessing coverage and variant calls."
  (:require [clojure.java.io :as io]
            [clj-yaml.core :as yaml]
            [bcbio.coverage.gene :as gene]
            [bcbio.coverage.io.bed :as bed]))

(defn- compare-region
  "Compare coverage based metrics and variations in a set of gene exons."
  [wgs-config exome-config coords])

(defn sqn-compare
  "Compare a WGS and exome experiment in specific gene regions of interest."
  [wgs-config exome-config region-file params]
  (let [gene-coord-file (gene/get-coord-bed region-file params)]
    (with-open [rdr (io/reader gene-coord-file)]
      (doseq [coords (map second (group-by :name (bed/get-iterator rdr)))]
        (compare-region wgs-config exome-config coords)))))
