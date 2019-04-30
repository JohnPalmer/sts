.calculate_distances <- function(i, this_D, nrow_this_D, data, x, y, t, ID, group, distance_function="NegExp") {
  these_indexes = (i+1):nrow_this_D
  this_pair_table = data.table(this_D[rep(i, length(these_indexes)), ], this_D[these_indexes, ], check.names = TRUE)
  if(distance_function=="NegExp"){
    this_pair_table[, dist:= exp(-1*geosphere::distCosine(cbind(get(x), get(y)), cbind(get(paste0(x, ".1")),get(paste0(y, ".1")))))]
  } else{
  this_pair_table[, dist:= geosphere::distCosine(cbind(get(x), get(y)), cbind(get(paste0(x, ".1")),get(paste0(y, ".1"))))]
  }
  this_pair_table[, group_combo:=paste(get(group), get(paste0(group, ".1")), sep="_")]
  this_pair_table[get(group) == get(paste0(group, ".1")), group_combo:=get(group)]
  this_pair_table[get(group) != get(paste0(group, ".1")), group_combo:="mix"]
  return(this_pair_table[, .(N=.N, distsum=sum(dist)), by=group_combo])
}

.cmp_calculate_distances = compiler::cmpfun(.calculate_distances)

.stp = function(data = D, x = "long", y = "lat", t = "time", ID = "ID", group = "race", distance_function="NegExp", cores=2, small=FALSE){
  if(!is(data, "data.table")){
    data = data.table(data)
  }
  data_list = split(data, by=t)
  dist_table = rbindlist(
    lapply(
      data_list, function(this_D) {
        nrow_this_D = nrow(this_D)
        if(nrow_this_D>1){
          if(small & nrow_this_D<=1000){
            these_indexes = combn(1:nrow_this_D, 2)
            this_pair_table = data.table(this_D[these_indexes[1,]], this_D[these_indexes[2,]], check.names = TRUE)
            if(distance_function=="NegExp"){
            this_pair_table[, dist:= exp(-1*geosphere::distCosine(cbind(get(x), get(y)), cbind(get(paste0(x, ".1")),get(paste0(y, ".1")))))]
            } else{
              this_pair_table[, dist:= geosphere::distCosine(cbind(get(x), get(y)), cbind(get(paste0(x, ".1")),get(paste0(y, ".1"))))]
            }
            this_pair_table[, group_combo:=paste(get(group), get(paste0(group, ".1")), sep="_")]
            this_pair_table[get(group) == get(paste0(group, ".1")), group_combo:=get(group)]
            this_pair_table[get(group) != get(paste0(group, ".1")), group_combo:="mix"]
            return(this_pair_table[, .(N=.N, distsum=sum(dist)), by=group_combo])

          } else{
          rbindlist(
          if(cores>1){
            parallel::mclapply(1:(nrow_this_D-1), function(i) {
          return(.cmp_calculate_distances(i=i, this_D=this_D, nrow_this_D=nrow_this_D, data=data, x=x, y=y, t=t, ID=ID, group=group, distance_function=distance_function))}, mc.cores=cores)}
          else{
            lapply(1:(nrow_this_D-1), function(i) {
              print(paste0(i, " of ", (nrow_this_D-1)))
              flush.console()
              return(.cmp_calculate_distances(i=i, this_D=this_D, nrow_this_D=nrow_this_D, data=data, x=x, y=y, t=t, ID=ID, group=group, distance_function=distance_function))})
            })
          }
        }
  }))
  dist_table2 = dist_table[, .(N=sum(N), distsum=sum(distsum)), by=group_combo]
  Pww = dist_table2[group_combo=="w", distsum/N]
  Pbb = dist_table2[group_combo=="b", distsum/N]
  Ptt = dist_table2[, .(distsum=sum(distsum), N=sum(N))][, distsum/N]
  Nw = data[time==0 & race=="w", .N]
  Nb = data[time==0 & race=="b", .N]
  Nt = Nw+Nb
  STP = (Nw*Pww + Nb*Pbb)/(Nt*Ptt)
  return(STP)
}

#' Spatio-Temporal Proximity Index
#'
#' \code{stp} calculates the spatio-temporal proximity index for a population in
#' which each individual is labelled as belonging to one of two groups.
#'
#' This function requires a set of pairwise distances between all individuals in the population and a variable that contains the group labels. It returns the spatio-temporal proximity index and, optionally, a set of boostrap replicates from which uncertainty may be estimated.
#' @export
stp = compiler::cmpfun(.stp)

#' Spatio-Temporal Proximity Index Bootstrapped
#'
#' \code{stpboot} calculates the spatio-temporal proximity index for a population in
#' which each individual is labelled as belonging to one of two groups.
#'
#' This function requires a set of pairwise distances between all individuals in the population and a variable that contains the group labels. It returns the spatio-temporal proximity index and, optionally, a set of boostrap replicates from which uncertainty may be estimated.
#' @export
stpboot = function(data = D, x = "long", y = "lat", t = "time", ID = "ID", group = "race", distance_function="NegExp", cores=2, small=FALSE, nboots=300){
  if(cores>1){
    unlist(parallel::mclapply(1:nboots, function(bootn) {
      nrow_data = nrow(data)
      rdata = data[sample(1:nrow_data, nrow_data, replace=TRUE)]
      return(stp(data=rdata, x=x, y=y, t=t, ID=ID, group=group, distance_function=distance_function, cores=1, small=small))
    }, mc.cores=cores))
      }else{
  sapply(1:nboots, function(bootn) {
    print(paste0("replicate ", bootn))
    flush.console()
    nrow_data = nrow(data)
    rdata = data[sample(1:nrow_data, nrow_data, replace=TRUE)]
    return(stp(data=rdata, x=x, y=y, t=t, ID=ID, group=group, distance_function=distance_function, cores=1, small=small))
  })
  }
}

#' Spatio-Temporal Proximity Index Empirical Sampling Distribution
#'
#' \code{stpesd} calculates the spatio-temporal proximity index for a population in
#' which each individual is labelled as belonging to one of two groups.
#'
#' This function requires a set of pairwise distances between all individuals in the population and a variable that contains the group labels. It returns the spatio-temporal proximity index and, optionally, a set of boostrap replicates from which uncertainty may be estimated.
#' @export
stpesd = function(data = D, x = "long", y = "lat", t = "time", ID = "ID", group = "race", distance_function="NegExp", cores=2, small=FALSE, sample_size=100, nreps=300){
  if(cores>1){
    unlist(parallel::mclapply(1:nreps, function(repn) {
      nrow_data = nrow(data)
      rdata = data[sample(1:nrow_data, sample_size, replace=TRUE)]
      return(stp(data=rdata, x=x, y=y, t=t, ID=ID, group=group, distance_function=distance_function, cores=1, small=small))
    }, mc.cores=cores))
  }else{

  sapply(1:nreps, function(repn) {
    print(paste0("replicate ", repn))
    flush.console()
    nrow_data = nrow(data)
    rdata = data[sample(1:nrow_data, sample_size, replace=TRUE)]
    return(stp(data=rdata, x=x, y=y, t=t, ID=ID, group=group, distance_function=distance_function, cores=cores, small=small))
  })
  }
}

#' Spatio-Temporal Proximity Index Bias Estimation
#'
#' \code{stpbias} calculates the spatio-temporal proximity index for a population in
#' which each individual is labelled as belonging to one of two groups.
#'
#' This function requires a set of pairwise distances between all individuals in the population and a variable that contains the group labels. It returns the spatio-temporal proximity index and, optionally, a set of boostrap replicates from which uncertainty may be estimated.
#' @export
stpbias = function(data = D, x = "long", y = "lat", t = "time", ID = "ID", group = "race", distance_function="NegExp", cores=2, small=FALSE, sample_size=100, nreps=300){
  pop_stp = stp(data=data, x=x, y=y, t=t, ID=ID, group=group, distance_function=distance_function, cores=cores, small=FALSE)
sample_sizes_result <- lapply(sample_size, function(this_size){
  esd = stpesd(data=data, x=x, y=y, t=t, ID=ID, group=group, distance_function=distance_function, cores=cores, small=small, sample_size=sample_size, nreps=nreps)
  pop_table = data.frame(data[, .N, by=get(group)])
  names(pop_table) = c("group", "N")
  return(list("bias"=mean(esd)-pop_stp, "esd"=esd, "mean_esd"=mean(esd), "pop_stp"=pop_stp, "sample_size"=sample_size, "nreps"=nreps, "pop_groups"=pop_table))
})
names(sample_sizes_result) <- paste0("sample_size_", sample_size)
return(sample_sizes_result)
  }
