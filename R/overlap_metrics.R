overlap_metrics <- function(compare_matrix, env_data, maha_distM, in_Ellipsoid) {
  get_metrics <- lapply(1:dim(compare_matrix)[2],
                        function(x){
                          compare <- in_Ellipsoid[, compare_matrix[,x]]
                          zeros1 <- which(rowSums(compare) %in% 0)
                          compare <- cbind(1:dim(compare)[1],compare)
                          compare <- compare[-zeros1,]
                          union_npoints <- dim(compare)[1]
                          union_prop <- colSums(compare[,-1])/union_npoints
                          union_sum <-rowSums(compare[,-1])
                          intersect_id <- which(union_sum==2)
                          np_intesection_global<- length(intersect_id)
                          intersection <- np_intesection_global/union_npoints
                          Ellip1_in_Ellpi2 <- union_prop[1]/union_prop[2]
                          Ellip2_in_Ellpi1 <- union_prop[2]/union_prop[1]


                          ellipsoid_over <- list(np_union = union_npoints,
                                                 np_intersection=np_intesection_global,
                                                 intersection_prop =   intersection,
                                                 Ellip1_in_Ellpi2 = Ellip1_in_Ellpi2,
                                                 Ellip2_in_Ellpi1 = Ellip2_in_Ellpi1)
                          ellipsoid_over <-as.data.frame(ellipsoid_over)


                          env_data_in<- env_data[compare[intersect_id,1],]
                          cov_center_inter_centroid <- hsi::cov_center(env_data_in,
                                                                       mve = T,
                                                                       level = 0.99999,
                                                                       vars =1:ncol(env_data_in))

                          mh_inter_centroid <- mahalanobis(x =env_data_in,
                                                           center =  cov_center_inter_centroid$centroid,
                                                           cov = cov_center_inter_centroid$covariance)


                          mh_inter_centroid <- exp(-0.5 *  mh_inter_centroid )


                          env_data_in <- data.frame(cbind(env_data_in,maha_distM[compare[intersect_id,1],
                                                                                 compare_matrix[,x]]),mh_inter_centroid)


                          return(list(ellipsoid_over,env_data_in))

                        })
  return(get_metrics)
}
