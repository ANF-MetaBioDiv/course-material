data(mite, package = "vegan")

inb_ind <- c("kl", "ch", "hartigan", "scott", "cindex", "db", "silhouette",
             "duda", "pseudot2", "beale", "ratkowsky", "ball", "ptbiserial",
             "gap", "frey", "mcclain", "gamma", "gplus", "tau", "dunn",
             "hubert", "sdindex", "dindex", "sdbw")

# Find best number of clusters
set.seed(1000)

do_nbclust <- function(data, index) {
    print(paste0("Trying ", index, " index..."))
    tmp <- try(NbClust::NbClust(data = data,
                    distance = "euclidean",
                    diss = NULL,
                    min.nc = 2,
                    max.nc = 5,
                    method = "average",
                    index = index), silent = TRUE)
    if ("Best.nc" %in% names(tmp)) tmp[["Best.nc"]][1]
}

num_clust <- lapply(inb_ind, function(x) do_nbclust(mite, x))
names(num_clust) <- inb_ind
num_clust <- unlist(num_clust)

num_clust

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

cat("Based on a number of criteria, we will select",
       getmode(num_clust),
       "clusters.\n")
