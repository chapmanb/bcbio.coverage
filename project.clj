(defproject bcbio.coverage "0.0.1-SNAPSHOT"
  :description "Investigate coverage metrics for variant calling experiments"
  :url "https://github.com/chapmanb/bcbio.coverage"
  :license {:name "MIT" :url "http://www.opensource.org/licenses/mit-license.html"}
  :repositories [["jensembl" {:url "http://jensembl.sourceforge.net/m2-repo" :checksum :ignore :snapshots false }]]
  :dependencies [[org.clojure/clojure "1.5.1"]
                 [org.clojure/data.csv "0.1.2"]
                 [org.clojure/tools.cli "0.2.2"]
                 [org.clojars.chapmanb/bigwig "r39"]
                 [org.clojars.chapmanb/picard "1.90"]
                 [org.clojars.chapmanb/sam "1.90"]
                 [org.clojars.chapmanb/tribble "1.90"]
                 [org.clojars.chapmanb/variant "1.90"]
                 [org.clojars.chapmanb/gatk-lite "2.6.5"]
                 [org.bioclojure/bio.ensembl "0.2.0"
                  :exclusions [uk.ac.roslin/ensembl-data-access uk.ac.roslin/ensembl-config]]
                 [uk.ac.roslin/ensembl-data-access "1.18"]
                 [uk.ac.roslin/ensembl-config "1.72"]
                 [clj-yaml "0.4.0"]
                 [incanter/incanter-core "1.5.0"]
                 [incanter/incanter-charts "1.5.0"]
                 [log4j "1.2.17"]
                 [org.slf4j/slf4j-log4j12 "1.7.5"]
                 [de.kotka/lazymap "3.1.1"]
                 [me.raynes/fs "1.4.3"]
                 [ordered "1.3.2" :exclusions [org.clojure/clojure]]
                 [lonocloud/synthread "1.0.4"]]
  :plugins [[lein-midje "3.0.1"]]
   ;; Remove X11 requirement for charts
   ;; http://www.jfree.org/phpBB2/viewtopic.php?t=1012
  :jvm-opts ["-Djava.awt.headless=true"]
  :profiles {:dev {:dependencies
                   [[midje "1.5.1"]]}}
  :main bcbio.coverage.main)
