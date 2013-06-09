(ns rosalind.core
  (:use [clojure.string :only (split)]))


;; Introduction to the Bioinformatics Armory

(defn freqs [strand]
  (for [sym [\A \C \G \T]]
    ((frequencies strand) sym)))

(freqs "AGACTCACCCGGGACTGCCTAAAAACGAGTCTGTACTCCTTGAGGAGAATTCACAGGAGGCCGTAGCGTTGATTGCCGCAACCCCGTCTGATCTTTTGGATCACTTGATGTTGCGTGGTATTTGGGATGGTAACACGTAGTTCTGTAGCGCGTGTGTGCTAGGGTAAGACTTATGGCTGTCAAGCAGCGTGGGCAAGCCCAGACCCGATCATCATGATGGCCTGACTGGCACTACGTATACAGGGATAGTCGTTCTTTATAAAGATTTAGAGCGTTCCCTGACGATGTCGGGCCTTGCTCTTAACCCCCAATCTCTAACCAGGGCGCGTCTCAAACTAGTGCGGCAAAGCGTCTTAACGTATCGTTGTTTGCTGAGCGCAACGGATCTACATCTCAGTTGACACAACCCTAAACGAGCAAAGTCTCTACGAATCATAACGGTCGCAACACGGGTGAACCTCAGGTTGTTAGCGTTGCTAACGAGTGTCGATTTTCATACTCTCTCCCCTGAGTATGCTTAATACAAGGGGCGAACCAGTTGATATGTGCGGTCTCGTAGGACAAGAATAGCCTGAATCGAGCATGACCCTTGTTTGTTTGGCAGATAGACACACGGAACTTACGCCCTGTGGAGTTCTTCTTCTACGTAGAAGGCACCACTACACTCGCTGCCCCGACATGAACACCCACACCAAGTGAATGAGATCCGCATCCATTTAGTGGGTCGCTAAACCCCATGCCCCGCGTCACATCCGTCCGGAGAAAAACATTTAGCTCCACTAGTTCCTCGGGTATTCGACCAAACACTTGGTTGATATATTGTGTACAATGCTAGGAGGTTAAGTAGGTTTCAGGGCCTTGGACATTGAGAGGGGGCCGGAAAGA")

;; Introduction to Protein Databases

(defn processes-from-db [protein]
  (->> (str "http://www.uniprot.org/uniprot/" protein ".txt")
       slurp
       (re-seq #"(?m)DR\s+GO;\s+GO:\d+; P:(.+?)\;")
       (map second)))

(doseq [s (processes-from-db "B1IMN7")]
  (println s))
