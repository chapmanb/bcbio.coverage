(ns bcbio.coverage.main
  "Main entry point for command line programs."
  (:require [clojure.string :as string]
            [bcbio.coverage.gene]
            [bcbio.coverage.workflow.multicompare :as multicompare]
            [bcbio.coverage.workflow.wgsexome :as wgsexome])
  (:gen-class))

(def ^{:doc "Mapping of command line arguments to sub commands"
       :private true}
  main-map
  {:gene bcbio.coverage.gene/-main
   :wgsexome wgsexome/-main
   :multicompare multicompare/-main})

(defn -main [& args]
  (if-let [main-fn (get main-map (keyword (first args)))]
    (apply main-fn (rest args))
    (do
      (println (str "Unexpected command: '" (first args) "' "
                    "Available commands: " (string/join ", " (map name (keys main-map)))))))
  (shutdown-agents)
  (System/exit 0))
