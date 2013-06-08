(defproject bcbio.coverage "0.0.1-SNAPSHOT"
  :description "Investigate coverage metrics for variant calling experiments"
  :url "https://github.com/chapmanb/bcbio.coverage"
  :license {:name "MIT" :url "http://www.opensource.org/licenses/mit-license.html"}
  :dependencies [[org.clojure/clojure "1.5.1"]
                 [org.clojure/tools.cli "0.2.2"]
                 [org.clojars.chapmanb/bigwig "r39"]
                 [org.clojars.chapmanb/picard "1.90"]
                 [org.clojars.chapmanb/sam "1.90"]
                 [org.clojars.chapmanb/tribble "1.90"]
                 [org.clojars.chapmanb/variant "1.90"]
                 [org.clojars.intronic/bio.ensembl "0.1.1-SNAPSHOT"]
                 [log4j "1.2.17"]
                 [de.kotka/lazymap "3.1.1"]
                 [me.raynes/fs "1.4.3"]
                 [ordered "1.3.2" :exclusions [org.clojure/clojure]]
                 [lonocloud/synthread "1.0.4"]]
  :plugins [[lein-midje "3.0.1"]]
  :profiles {:dev {:dependencies
                   [[midje "1.5.1"]]}}
  :main bcbio.coverage.main)
