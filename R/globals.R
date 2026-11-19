# column names used inside ggplot2 aes() and data-frame helpers
utils::globalVariables(c(".F", ".G", ".ID", ".L", ".LG", ".PANEL", ".PNL", ".U", ".X", ".XC", ".Y",
                         ".fc", ".fr", ".fw", ".x", ".y", "A", "B", "C", "CONC_USED", "GROUP", "L",
                         "PROFILE", "TIME", "U", "r"))

# a small cache of values that are expensive to build and never change
.hp_cache <- new.env(parent = emptyenv())
