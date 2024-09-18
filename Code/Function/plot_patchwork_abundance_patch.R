#' Plot patchwork abundance when given it.
#'
#' @param result_list1 From the desolve simulator function
#' @param result_list2 From the desolve simulator function
#'
#' @return
#' @export
#'
#' @examples
plot_patchwork_abundance <- function(result_list1, result_list2 = NA) {

    plot_titles <- c(
    "LOW COV/HIGH FREQ",
    "HIGH COV/HIGH FREQ",
    "LOW COV/LOW FREQ",
    "HIGH COV/LOW FREQ"
  )

  if (length(result_list2) == 1) {
    plots <- lapply(1:4, function(i) {
      plot_abundance_patch_time(result_list1[[i]][[1]],
        aggregate = "yes"
      ) +
        ggtitle(plot_titles[i])
    })
  } else {
    plots <- lapply(1:4, function(i) {
      plot_abundance_patch_time(result_list1[[i]][[1]],
        aggregate = "yes"
      ) +
        plot_abundance_patch_time(result_list2[[i]][[1]],
          aggregate = "yes"
        ) +
        ggtitle(plot_titles[i])
    })
  }


  wrap_plots(plots) + plot_layout(guides = "collect")
}
