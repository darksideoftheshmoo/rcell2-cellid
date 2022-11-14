test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("load_tsv_masks: load from examples", {
  # Example: Choose one set of example images:
  data.dir <- system.file("extdata/sample_datasets/sample_time_series/Position4/out_all_masks.tsv.gz",
                          package = "rcell2.examples")
  # Load
  d <- read_coords_tsv(masks_tsv_path = data.dir)
  
  # Load original
  ref.tsv <- system.file("extdata/test_data/out_all_masks.RDS",
                         package = "rcell2.examples")
  
  # Test
  expect_equal(d, ref.tsv)
})
